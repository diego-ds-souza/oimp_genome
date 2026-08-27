#!/usr/bin/env bash
#
# Sort every orthogroup into the nine categories reported in the paper:
# single-copy universal, multicopy universal, species-specific, Cerambycidae-,
# Chrysomelidae-, Curculionoidea- and Phytophaga-specific, other shared
# orthologs, and unassigned genes.
#
# Family membership comes from species_taxonomy.tsv rather than from name
# prefixes, so adding or renaming a species does not silently change the
# categories.
#
# Run from the repository root:
#   conda activate oimp_04_orthology
#   bash 04_orthology/08_orthology_categories.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
PROJECT="${PROJECT:-results/04_orthology}"          # section output root
OF_BASE="${OF_BASE:-${PROJECT}/orthofinder_out}"    # OrthoFinder output, step 03
FA_DIR="${FA_DIR:-${PROJECT}/clean_proteomes}"      # collapsed proteomes, step 02
TREE_NEX="${TREE_NEX:-${PROJECT}/concat/partitions_aa.txt_rooted_v2.nex}"
                                                    # rooted tree, from step 07
TAXONOMY_TSV="${TAXONOMY_TSV:-${PROJECT}/species_taxonomy.tsv}"
                                                    # family per species, step 07
OUT_DIR="${OUT_DIR:-${PROJECT}/plots}"              # output directory
ORDER_TXT="${ORDER_TXT:-${OUT_DIR}/species_order_manual.txt}"
# -----------------------------------------------------------------------------

mkdir -p "$OUT_DIR"

command -v python3 >/dev/null || {
    echo "ERROR: python3 not found"
    exit 1
}

[[ -f "$TAXONOMY_TSV" ]] || {

    echo "ERROR: taxonomy file not found:"
    echo "$TAXONOMY_TSV"

    exit 1
}

# 1) Locate the OrthoFinder results, preferring the path step 03 recorded.
if [[ -f "${OF_BASE}/LATEST_RESULTS_PATH.txt" ]]; then

    RESULTS_DIR=$(cat "${OF_BASE}/LATEST_RESULTS_PATH.txt")

else

    RESULTS_DIR=$(
        ls -d ${OF_BASE}/run_*/Results_* 2>/dev/null |
        sort |
        tail -n1 ||
        true
    )

    if [[ -z "${RESULTS_DIR}" ]]; then

        RESULTS_DIR=$(
            ls -d ${OF_BASE}/Results_* 2>/dev/null |
            sort |
            tail -n1 ||
            true
        )

    fi

fi

if [[ ! -d "${RESULTS_DIR}" ]]; then

    echo "ERROR: Could not find OrthoFinder Results directory"
    exit 1

fi

GENECOUNT="${RESULTS_DIR}/Orthogroups/Orthogroups.GeneCount.tsv"

if [[ ! -f "$GENECOUNT" ]]; then

    echo "ERROR: Gene count file not found:"
    echo "$GENECOUNT"

    exit 1

fi

echo "Using OrthoFinder results:"
echo "$RESULTS_DIR"
echo

# 2) Classify every orthogroup and write the per-species summary tables.
python3 - \
"$GENECOUNT" \
"$FA_DIR" \
"$TREE_NEX" \
"$TAXONOMY_TSV" \
"$OUT_DIR" \
"$ORDER_TXT" <<'PY'

import sys
import os
import re
import glob
import pandas as pd

gene_count_file = sys.argv[1]
fa_dir          = sys.argv[2]
tree_nexus      = sys.argv[3]
taxonomy_file   = sys.argv[4]
out_dir         = sys.argv[5]
order_file      = sys.argv[6]

# Helpers: read FASTA headers, and normalise species names.
def normalize_name(x):

    """
    Normalize names between OrthoFinder and FASTA files.
    """

    x = str(x).strip().strip("'")

    x = re.sub(
        r'\.(fa|faa|fasta)$',
        '',
        x,
        flags=re.I
    )

    x = re.sub(
        r'_protein$',
        '',
        x,
        flags=re.I
    )

    return x

def find_fasta(species, raw_name):

    """
    Find matching proteome file.

    Accepts:

        species.fa
        species.faa
        species.fasta

    and wildcard versions.
    """

    candidates = [

        f"{raw_name}.fa",
        f"{raw_name}.faa",
        f"{raw_name}.fasta",

        f"{raw_name}*.fa",
        f"{raw_name}*.faa",
        f"{raw_name}*.fasta",

        f"{species}.fa",
        f"{species}.faa",
        f"{species}.fasta",

        f"{species}*.fa",
        f"{species}*.faa",
        f"{species}*.fasta",

    ]

    for pattern in candidates:

        files = sorted(
            glob.glob(
                os.path.join(
                    fa_dir,
                    pattern
                )
            )
        )

        if files:

            return files[0]

    return None

def count_fasta(path):

    return sum(
        1
        for line in open(path)
        if line.startswith(">")
    )

# 1) Read the OrthoFinder gene-count matrix.
og = pd.read_csv(
    gene_count_file,
    sep="\t"
)

orf_columns = [

    c
    for c in og.columns
    if c not in [
        "Orthogroup",
        "Total"
    ]

]

species_columns = {}

for col in orf_columns:

    species = normalize_name(col)

    species_columns[species] = col

print("OrthoFinder species detected:")
for s in species_columns:

    print(" ", s)

print()

# 2) Match each species column to its proteome file.
fasta_files = {}

for species, raw in species_columns.items():

    fasta_files[species] = find_fasta(
        species,
        raw
    )

print("FASTA matches:")

for species, fasta in fasta_files.items():

    print(
        f" {species} -> {fasta}"
    )

print()

species = [

    s
    for s in species_columns
    if fasta_files[s] is not None

]

if len(species) == 0:

    raise SystemExit(
        "No OrthoFinder species columns matched FASTAs."
    )

print(
    "Matched species:",
    len(species)
)

print()

# 3) Count the proteins contributed by each species.
protein_counts = {

    sp: count_fasta(
        fasta_files[sp]
    )

    for sp in species

}

# 4) Read the family assignment of each species.
taxonomy = pd.read_csv(
    taxonomy_file,
    sep="\t",
    dtype=str
).fillna("")

required_columns = {
    "Species",
    "Family",
    "Superfamily"
}

missing = required_columns - set(taxonomy.columns)

if missing:

    raise SystemExit(
        "Taxonomy table missing columns: "
        + ", ".join(sorted(missing))
    )

taxonomy["Species"] = taxonomy["Species"].map(
    normalize_name
)

taxonomy = taxonomy.set_index(
    "Species"
)

missing_taxonomy = [

    sp
    for sp in species
    if sp not in taxonomy.index

]

if missing_taxonomy:

    raise SystemExit(
        "Species missing from taxonomy table:\n"
        +
        "\n".join(missing_taxonomy)
    )

# Taxonomic groups: Curculionoidea is a superfamily, Cerambycidae and
# Chrysomelidae are families.

def get_group(sp):

    family = taxonomy.loc[
        sp,
        "Family"
    ]

    superfamily = taxonomy.loc[
        sp,
        "Superfamily"
    ]

    if superfamily == "Curculionoidea":

        return "Curculionoidea"

    if family == "Cerambycidae":

        return "Cerambycidae"

    if family == "Chrysomelidae":

        return "Chrysomelidae"

    return "Outgroup"

groups = {

    sp: get_group(sp)

    for sp in species

}

curculionoidea = [

    sp
    for sp in species
    if groups[sp] == "Curculionoidea"

]

cerambycidae = [

    sp
    for sp in species
    if groups[sp] == "Cerambycidae"

]

chrysomelidae = [

    sp
    for sp in species
    if groups[sp] == "Chrysomelidae"

]

outgroup = [

    sp
    for sp in species
    if groups[sp] == "Outgroup"

]

print("Groups:")

print(
    " Curculionoidea:",
    len(curculionoidea)
)

print(
    " Cerambycidae:",
    len(cerambycidae)
)

print(
    " Chrysomelidae:",
    len(chrysomelidae)
)

print(
    " Outgroup:",
    len(outgroup)
)

print()

# 5) Define the orthology categories.
categories = [

    "SingleCopy_Universal",

    "Multicopy_Universal",

    "Species_Specific",

    "Curculionoidea_Specific",

    "Cerambycidae_Specific",

    "Chrysomelidae_Specific",

    "Phytophaga_Specific",

    "Other_Orthologs",

    "Unassigned_Genes"

]

summary = pd.DataFrame(

    0,

    index=species,

    columns=categories

)

# Helpers for the category tests.
def get_counts(row):

    return {

        sp: int(row[species_columns[sp]])

        for sp in species

    }

def present_species(counts):

    return {

        sp
        for sp, value in counts.items()
        if value > 0

    }

# 6) Assign every orthogroup to a category.
for _, row in og.iterrows():

    if row["Orthogroup"] == "Total":

        continue

    counts = get_counts(row)

    present = present_species(counts)

    if len(present) == 0:

        continue

    # ------------------------------------------------------------
    # Universal orthogroups
    # ------------------------------------------------------------

    if len(present) == len(species):

        if all(
            counts[sp] == 1
            for sp in species
        ):

            summary[
                "SingleCopy_Universal"
            ] += 1

        else:

            for sp in species:

                summary.loc[
                    sp,
                    "Multicopy_Universal"
                ] += counts[sp]

        continue

    # ------------------------------------------------------------
    # Species-specific
    # ------------------------------------------------------------

    if len(present) == 1:

        sp = next(iter(present))

        summary.loc[
            sp,
            "Species_Specific"
        ] += counts[sp]

        continue

    # ------------------------------------------------------------
    # Group membership
    # ------------------------------------------------------------

    in_curc = any(
        sp in present
        for sp in curculionoidea
    )

    in_cer = any(
        sp in present
        for sp in cerambycidae
    )

    in_chry = any(
        sp in present
        for sp in chrysomelidae
    )

    in_out = any(
        sp in present
        for sp in outgroup
    )

    # ------------------------------------------------------------
    # Cerambycidae-specific
    # ------------------------------------------------------------

    if (
        in_cer
        and not in_curc
        and not in_chry
        and not in_out
    ):

        for sp in cerambycidae:

            summary.loc[
                sp,
                "Cerambycidae_Specific"
            ] += counts.get(sp,0)

        continue

    # ------------------------------------------------------------
    # Chrysomelidae-specific
    # ------------------------------------------------------------

    if (
        in_chry
        and not in_curc
        and not in_cer
        and not in_out
    ):

        for sp in chrysomelidae:

            summary.loc[
                sp,
                "Chrysomelidae_Specific"
            ] += counts.get(sp,0)

        continue

    # ------------------------------------------------------------
    # Curculionoidea-specific
    # ------------------------------------------------------------

    if (
        in_curc
        and not in_cer
        and not in_chry
        and not in_out
    ):

        for sp in curculionoidea:

            summary.loc[
                sp,
                "Curculionoidea_Specific"
            ] += counts.get(sp,0)

        continue

    # ------------------------------------------------------------
    # Phytophaga-specific
    # ------------------------------------------------------------

    if (
        in_curc
        and (in_cer or in_chry)
        and not in_out
    ):

        for sp in present:

            summary.loc[
                sp,
                "Phytophaga_Specific"
            ] += counts[sp]

        continue

    # ------------------------------------------------------------
    # Other orthologs
    # ------------------------------------------------------------

    for sp in present:

        summary.loc[
            sp,
            "Other_Orthologs"
        ] += counts[sp]

# 7) Count the genes that fall in no orthogroup.
assigned = pd.Series(
    0,
    index=species,
    dtype=int
)

for _, row in og.iterrows():

    if row["Orthogroup"] == "Total":

        continue

    for sp in species:

        assigned[sp] += int(
            row[species_columns[sp]]
        )

for sp in species:

    summary.loc[
        sp,
        "Unassigned_Genes"
    ] = max(
        0,
        protein_counts[sp] - assigned[sp]
    )

# 8) Order the species to match the tips of the species tree.
def extract_nexus_order(path):

    if not os.path.isfile(path):

        return []

    text = open(path).read()

    labels = []

    for token in re.split(
        r"[\(\),:;=\s]+",
        text
    ):

        token = normalize_name(token)

        if token in species:

            if token not in labels:

                labels.append(token)

    return labels

if os.path.isfile(order_file):

    order = [

        normalize_name(x.strip())

        for x in open(order_file)

        if x.strip()

    ]

elif os.path.isfile(tree_nexus):

    order = extract_nexus_order(
        tree_nexus
    )

else:

    order = species

order = [

    sp
    for sp in order
    if sp in species

]

for sp in species:

    if sp not in order:

        order.append(sp)

summary = summary.loc[order]

# 9) Write the long and wide count tables.
long_table = (

    summary
    .reset_index()
    .melt(
        id_vars="index",
        var_name="Category",
        value_name="Count"
    )
    .rename(
        columns={
            "index":"Species"
        }
    )

)

long_table.to_csv(

    os.path.join(
        out_dir,
        "orthology_histogram_counts_long.csv"
    ),

    index=False

)

summary.reset_index().rename(

    columns={
        "index":"Species"
    }

).to_csv(

    os.path.join(
        out_dir,
        "orthology_histogram_counts_wide.csv"
    ),

    index=False

)

with open(

    os.path.join(
        out_dir,
        "orthology_tip_order.txt"
    ),

    "w"

) as f:

    for sp in summary.index:

        f.write(
            sp + "\n"
        )

print()
print("Finished successfully.")
print("Outputs written to:")
print(out_dir)

PY

echo "Done."
