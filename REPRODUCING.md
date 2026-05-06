# Reproducing the Manuscript Pipeline

The canonical entry point is:

```sh
Rscript R/run_manuscript_pipeline.R
```

The runner:

- validates literal `source()` paths across `R/`
- links local input data from `DARTER_DATA_ROOT`
- runs scripts in the order of the current methods/results outline
- logs task status under `Outputs/manuscript_runs/<run_id>/`
- compares generated output filenames against `Documents/darter_figures`

## Input Data

By default, local input data are read from:

```text
~/Documents/darter_morphometrics
```

Override this with:

```sh
DARTER_DATA_ROOT=/path/to/darter_morphometrics Rscript R/run_manuscript_pipeline.R
```

The runner expects these inputs in the data root:

- `darter_curves.txt`
- `landmarks_ref.txt`
- `side_shapes/`

Optional analyses also use:

- `photos/`
- `trait_measurements/`
- `papillae.csv`

These inputs are symlinked into the repo working directory and ignored by Git.

## Checks

Static validation only:

```sh
Rscript R/run_manuscript_pipeline.R --check-only
```

Stop immediately on required failures:

```sh
Rscript R/run_manuscript_pipeline.R --stop-on-error
```

Include archive/exploratory analyses:

```sh
Rscript R/run_manuscript_pipeline.R --include-archive
```

## Organization

See `MANUSCRIPT_PIPELINE.md` for the mapping between the outline and script locations.

## Current Validation

On 2026-05-06, after reorganizing around `outline methods and results 12.19.16 PM.txt`:

- `Rscript R/run_manuscript_pipeline.R --check-only` passed.
- All 80 R files parsed successfully.
- The default task plan contains 59 runnable tasks.
- `R/01_methods/01_size_allometry/05_combined_allometry_procD.R` completed.
- `R/02_results/05_ct_reference_distribution/03_ct_reference_mahalanobis_rarefaction.R` completed with `DARTER_RAREFACTION_ITER=10` as a smoke test.
- `R/02_results/09_modularity_integration/99_modularity_multiple_testing_correction.R` completed using existing timestamped modularity outputs.
