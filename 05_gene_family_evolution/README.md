# 05 — Gene family evolution

Methods section 6, second part: testing with CAFE whether orthogroup size
differences exceed what a neutral birth–death process would produce.

## Scripts

| Script | Software | Purpose |
|---|---|---|
| `01_prepare_cafe_input.py` | — | Build the CAFE count table and apply its two filters |
| `02_make_ultrametric.R` | ape | Convert the ML tree to a chronogram calibrated at 227 Ma |
| `03_fix_tree.py` | — | Lift numerically-zero branches that CAFE5 cannot evaluate |
| `04_check_tree.R` | ape | Diagnose a tree CAFE5 refuses to use |
| `05_run_cafe.sh` | CAFE | Fit the four models |
| `06_plot_cafe_column.py` | matplotlib | Expansion/contraction column for Fig. 2 |
| `07_summarise_cafe.py` | — | Tie the result to the annotated gene families |

## Running

```bash
conda env create -f 05_gene_family_evolution/environment.yml
conda activate oimp_05_gene_family_evolution

python  05_gene_family_evolution/01_prepare_cafe_input.py
Rscript 05_gene_family_evolution/02_make_ultrametric.R
python  05_gene_family_evolution/03_fix_tree.py
Rscript 05_gene_family_evolution/04_check_tree.R
bash    05_gene_family_evolution/05_run_cafe.sh
python  05_gene_family_evolution/06_plot_cafe_column.py
python  05_gene_family_evolution/07_summarise_cafe.py
```

`cafe` on bioconda is linux-64 only.

## Notes

- Two filters are applied when the count table is built, both from the CAFE5
  documentation: orthogroups absent from one of the two clades descending from
  the root are removed, since CAFE5 assumes every family was present in the
  common ancestor; and orthogroups with **≥ 30 genes in a single species** are
  set aside, because their size range makes every likelihood infinite. The
  `large` run in step 05 attempts the latter separately and is expected to fail.
- CAFE5 needs a rooted, binary, ultrametric tree; IQ-TREE gives a phylogram.
  Step 02 roots on *T. castaneum*, drops node labels, applies penalized-
  likelihood rate smoothing with `ape::chronos()` and scales the root to 227 Ma.
- Step 03 lifts any branch below 0.25 Ma to that floor. Rooting with
  `resolve.root = TRUE` can leave a branch of ~1e-14, and the birth–death
  transition probabilities go to 0/0 as t → 0, making every family likelihood
  infinite; because the branch is positive, a naive check misses it. Step 04
  prints the diagnostics that reveal the problem.
- Four models are fitted, all with `-p` so the root size follows a Poisson
  rather than CAFE5's default uniform distribution: `base` (one rate λ),
  `base_error` (the same plus an error model estimated from the data),
  `gamma_k3` (three gamma rate categories, as a consistency check) and `large`.
  **`base_error` is the run reported**, because the error model keeps
  annotation noise from being read as gene family change.
- A family-wide *p* < 0.05 marks an orthogroup evolving inconsistently with the
  genome-wide background; branch-specific *p*-values assign a change to the
  *O. impluviata* branch in particular.
