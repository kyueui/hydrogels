# CLAUDE.md

Guidance for AI assistants (and humans) working in this repository.

## What this repo is

Research code accompanying the paper **"Data-Driven De Novo Design of
Super-Adhesive Hydrogels."** It is **not** a software product or library — it is
a collection of reproducibility scripts, datasets, and analysis notebooks. The
goal of the work is to use machine learning + Bayesian optimization (SMBO) to
propose new hydrogel monomer formulas with maximized adhesive strength
("F_a", recorded in the data as `Glass (kPa)`), validated through iterative
rounds of wet-lab experiments.

There are two largely independent parts:

1. **Root directory** — the ML / sequential model-based optimization (SMBO)
   pipeline that proposes hydrogel formulas.
2. **`data_mining/`** — bioinformatics analysis of natural adhesive proteins
   (sequence alignment, consensus, amino-acid composition) used to derive the
   initial 180 hydrogel formulas that seed the ML pipeline.

## Repository layout

```
.
├── *.py                       # SMBO optimization scripts (run directly)
├── hu_utils.py                # Core shared library (models, plotting, CV, log parsing)
├── plt_configs.py             # Global matplotlib rcParams (Arial, tick styling)
├── requirement.txt            # Pinned dependencies (note: NOT requirements.txt)
├── *.xlsx                     # Raw experimental datasets (input + intermediate results)
├── round{1,2,3,4}_*_eval.ipynb # Benchmarking / correlation analysis notebooks (large)
├── data/                      # Processed CSV datasets + pickled RMSE dicts for rounds
│   ├── df_180.csv  df_289.csv  df_316.csv  df_341.csv
│   └── rmse_dict_*.pkl
├── reproduce/                 # 1000-repeat SMBO logs + scripts to tally chosen formulas
│   ├── count_logs_by_EI.py
│   ├── count_logs_by_PRED.py
│   └── *_saved_models.log     # Large (15–25 MB) precomputed run logs
├── data_mining/               # Protein bioinformatics pipeline (see its README.md)
└── SI-consensus-sequence.zip  # Supplementary consensus sequences
```

## Environment & dependencies

- Originally tested on **Python 3.8.3 / Ubuntu 18.04** (this container runs a
  newer Python; expect version-pin friction if reinstalling).
- Install with `pip install -r requirement.txt` (filename is
  `requirement.txt`, singular — easy to mistype).
- Key pinned versions matter for reproducibility: `scikit-learn==1.0.2`,
  `numpy==1.22`, `pandas==1.3.5`, `scipy==1.4.1`, `xgboost==1.6.2`.
- The optimization stack is **scikit-optimize** (`skopt`) on top of
  scikit-learn estimators. `GPy`/`GPyOpt` are also listed. `data_mining` uses
  **biopython** plus external `clustalo` / `clustalw` binaries (vendored in
  `data_mining/`).

## The ML / SMBO pipeline (root scripts)

Every optimization script follows the **same template**:

1. Load a dataset (round 1 reads `Original Data_ML_20220829.xlsx`; later rounds
   read `data/df_289.csv`, `df_316.csv`, etc.).
2. Define the 6 monomer-ratio features and the target:
   ```python
   x_columns = ["Nucleophilic-HEA", "Hydrophobic-BA", "Acidic-CBEA",
                "Cationic-ATAC", "Aromatic-PEA", "Amide-AAm"]
   y_column  = "Glass (kPa)"          # round 1
   # later rounds use "Glass (kPa)_max"
   ```
3. Build a model via `hu_utils.setup_gridsearch_model(<name>)` and `.fit(X, y)`.
4. Negate `y` (the optimizers **minimize**, so adhesion is maximized by
   minimizing `-y`).
5. Run `skopt.gp_minimize` with `acq_func="EI"`, seeding `x0`/`y0` with the real
   data, then print the candidate formulas sorted by acquisition value.

The scripts differ only in **which model proposes hypothetical values and which
maximizes EI**, and in the EI/penalty strategy:

| Script | Strategy (per README "Usage") |
|---|---|
| `CLmax_gp.py` | GP with constant-liar = `y_max` |
| `CLmin_gp.py` | GP with constant-liar = `y_min` |
| `LP_gp.py` | GP with a local-penalization term in EI |
| `gp_gp.py` | GP predictions as hypothetical values (kriging-believer) |
| `gp_rfr.py` | GP hypothetical values + RFR EI maximizer |
| `rfr_rfr.py` | RFR hypothetical values + RFR EI maximizer |
| `rfr_gp.py` | RFR hypothetical values + GP EI maximizer |
| `rfr_gp_star.py` | RFR-GP with warm start (10 RFR-generated points added) |
| `ENU_gp.py` / `ENU_rfr.py` | Brute-force **enumeration** (sample ~1e7 random L1-normalized formulas, predict, take top-1000) instead of SMBO |
| `*_round2.py` / `*_round3.py` | Same methods rerun on the round-2 (`df_289`) and round-3 (`df_316`) datasets |

Run any of them directly, e.g.:
```bash
python rfr_gp.py
```
Output is a NumPy argsort index list followed by lines of
`<id> [formula ratios] ML predicted: <value>`. Runtime is typically <10 min.

### Reproducing the paper's chosen formulas

`skopt`'s SMBO is numerically nondeterministic across runs, so the paper reran
each script **1000×** and picked the most frequent formulas. Those aggregate
logs live in `reproduce/*.log`. Tally them with the `fire`-based CLIs:

```bash
cd reproduce
python count_logs_by_EI.py   --round_num=1 --cutoff=50 --show_num=50 rfr_gp_saved_models.log
python count_logs_by_PRED.py --round_num=1 --cutoff=50 --show_num=50 rfr_gp_saved_models.log
```
`round_num` maps to a dataset-size threshold (1→180, 2→289, 3→316); only rows
with `id >= threshold` (i.e. newly proposed points) are counted. Formulas with
frequency >100 are flagged as "chosen for experiments."

## `hu_utils.py` — the shared library

A single ~1600-line module imported by every script and notebook (always via
`importlib.reload(hu_utils)` so edits take effect in long Jupyter sessions).
Key entry points:

- `setup_gridsearch_model(model, rs=1126)` — returns a `GridSearchCV` for a
  given model string: `"RFR"`/`"RFRsk"` (skopt vs sklearn RF), `"ETR"`, `"XGB"`,
  `"GP"`, `"LASSO"`, `"RIDGE"`, `"KRR"`, `"SVR"`, `"KNN"`, `"MLP"`,
  `"Dummy"`, etc. This is the canonical place to add/tune a model.
- `BestEstimatorCV` — wraps a fitted best estimator with CV stats and hold-out
  plotting (`output_stats`, `plot_hold_out`, `plot_importance`).
- `normalize(ratio_arr)` — L1-normalize a formula so its 6 ratios sum to 1.
- `pred_ei_x_to_df*` / `load_logs*` / `return_sorted_df` — parse SMBO result
  objects and saved logs into DataFrames.
- `plot_compared_methods_rmse[_180/_289/_316/_341]` — per-round benchmark bar
  plots (used by the eval notebooks).
- `plot_shap_waterfall`, `print_umap` — SHAP and UMAP visualizations.

Warnings are globally silenced (`warnings.filterwarnings("ignore")`).

## Datasets & rounds

The experimental campaign grew the dataset across 4 rounds. Sizes encode the
round and appear throughout filenames:

- **180** points (round 1) → `Original Data_ML_20220829.xlsx`, `data/df_180.csv`
- **289** (round 2) → `data/df_289.csv`
- **316** (round 3) → `data/df_316.csv`
- **341** (round 4) → `data/df_341.csv`, `data/df_341_ML_methods.csv`

`*_verified_*.xlsx`, `*_red_*.xlsx`, and `ML_ei&pred*.xlsx` are intermediate /
prediction spreadsheets. The `round{1..4}_*_eval.ipynb` notebooks hold the
model benchmarking and correlation analysis for each round (they are large —
round 1 is ~9 MB with embedded outputs).

## `data_mining/` — protein bioinformatics

Independent pipeline that produces the initial 180 hydrogel formulas from
natural adhesive-protein sequences. See `data_mining/README.md`. Main flow:

1. `gb_analysis_200.ipynb` — read raw `.gp` protein sequences (in
   `data_mining/sequence/`), select top-200 species, write FASTA to `fa_files/`.
2. `fa_files/run_clustalo.sh` — multiple-sequence alignment with `clustalo`;
   outputs to `aln_files/` and consensus to `con_files/`.
3. `histogram_visualize_200.ipynb` — pairwise amino-acid composition analysis →
   the 180 monomer proportions used to seed the wet-lab dataset.

This directory contains many parallel `*_files_<taxon>/` folders (per protein
class: Bivalvia, Polychaeta, Ascidiacea, fibrous, resilin, etc.) plus vendored
clustal binaries and figure outputs. Treat it as a data/asset tree, not source
to refactor.

## Conventions & gotchas

- **Maximize via minimize**: `y` is always negated before optimization; remember
  to flip the sign (`-res.func_vals`) when reading results.
- **Formulas are L1-normalized** 6-vectors (ratios summing to 1). Always
  normalize before predicting.
- **Reproducibility seeds** are set explicitly and inconsistently across scripts
  (`np.random.seed(0)`, `random.seed(0)`, then `np.random.seed(237)`; `rs=929`
  or `1126`). Don't "clean these up" — exact values affect published results.
- **`importlib.reload`** of `hu_utils` and `plt_configs` at the top of every
  script is intentional (notebook iteration), not redundant.
- Scripts assume the **repo root as CWD** (relative paths like
  `"Original Data_ML_20220829.xlsx"`, `"data/df_289.csv"`, `pca_figs/...`). Some
  plotting helpers write to directories (e.g. `pca_figs/`) that may need to be
  created first.
- Code style uses **black-formatted** Python with Japanese comments in
  `plt_configs.py`. Match the existing style; keep diffs minimal.
- This is research code with low test coverage and no CI — verify changes by
  running the relevant script/notebook end to end rather than relying on tests.

## Git workflow

- License: MIT.
- Develop on the branch designated for your current task (or, absent one, a
  dedicated feature branch); commit and push there with
  `git push -u origin <branch>`. Do not push to `main` without explicit
  permission.
- Do not open a pull request unless explicitly asked.
- Large binary artifacts (`.xlsx`, `.log`, `.zip`, big notebooks) are committed
  to the repo; be deliberate before adding more.
