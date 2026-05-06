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
- copies newly created or updated outputs into `manuscript_outputs/<run_id>/` by outline section
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

Non-used scripts are retained under `R/non_used_scripts/` and are excluded from default validation.
To include them in source-path validation:

```sh
Rscript R/run_manuscript_pipeline.R --include-non-used
```

## Organization

See `MANUSCRIPT_PIPELINE.md` for the mapping between the outline and script locations.

## Current Validation

On 2026-05-06, after reorganizing around `outline methods and results 12.19.16 PM.txt`:

- The repository was reorganized into `R/methods/`, `R/results/`, `R/supplemental/`, and `R/non_used_scripts/`.
- The default task plan contains the runnable scripts required by the current outline.
- New runs create sectioned output copies under `manuscript_outputs/<run_id>/`.
