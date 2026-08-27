#!/usr/bin/env bash
#
# Download the annotated proteomes of the 24 comparison beetle genomes from
# NCBI. Accessions are listed in Table S2 and in accessions.txt beside this
# script; the O. impluviata proteome comes from section 03 and is added to the
# same directory by hand.
#
# Run from the repository root:
#   conda activate oimp_04_orthology
#   bash 04_orthology/00_download_proteomes.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
ACC="${ACC:-04_orthology/accessions.txt}"      # one assembly accession per line
OUTDIR="${OUTDIR:-data/proteomes}"             # where the proteomes land
# -----------------------------------------------------------------------------

for t in datasets unzip; do
  command -v "$t" >/dev/null || { echo "ERROR: $t not found" >&2; exit 1; }
done
[[ -s "$ACC" ]] || { echo "ERROR: no such file: $ACC" >&2; exit 1; }

mkdir -p "$OUTDIR"

# 1) Fetch protein FASTA and genomic GFF3 for each accession. The GFF3 is
#    needed in step 02, which maps each protein back to its parent gene.
while read -r line; do
  acc="${line%%[[:space:]]*}"                  # strip the trailing "# Species" note
  [[ -z "$acc" || "$acc" == \#* ]] && continue
  if [[ -s "${OUTDIR}/${acc}/protein.faa" ]]; then
    echo "have ${acc}, skipping"; continue
  fi
  echo "[fetch] ${acc}"
  datasets download genome accession "$acc" \
           --include protein,gff3 --filename "${OUTDIR}/${acc}.zip"
  unzip -qo "${OUTDIR}/${acc}.zip" -d "${OUTDIR}/${acc}_tmp"
  mkdir -p "${OUTDIR}/${acc}"
  find "${OUTDIR}/${acc}_tmp" -name 'protein.faa'   -exec mv {} "${OUTDIR}/${acc}/" \;
  find "${OUTDIR}/${acc}_tmp" -name 'genomic.gff'   -exec mv {} "${OUTDIR}/${acc}/" \;
  rm -rf "${OUTDIR}/${acc}_tmp" "${OUTDIR}/${acc}.zip"
done < "$ACC"

echo "wrote ${OUTDIR}/<accession>/{protein.faa,genomic.gff}"
