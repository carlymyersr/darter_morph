# Reproducing Figures

The master reproduction entry point is:

```sh
Rscript R/master_reproduce_figures.R
```

The runner:

- validates literal `source()` paths across `R/`
- links local input data from `DARTER_DATA_ROOT`
- runs the canonical pipeline and figure/stat scripts in sequence
- logs task status under `Outputs/master_runs/<run_id>/`
- compares generated output filenames against `Documents/darter_figures`

## Input Data

By default, local input data are read from:

```text
~/Documents/darter_morphometrics
```

Override this with:

```sh
DARTER_DATA_ROOT=/path/to/darter_morphometrics Rscript R/master_reproduce_figures.R
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
Rscript R/master_reproduce_figures.R --check-only
```

Stop immediately on required failures:

```sh
Rscript R/master_reproduce_figures.R --stop-on-error
```

Interactive landmark digitizing is excluded by default. Include it only when running from an interactive R-capable environment:

```sh
Rscript R/master_reproduce_figures.R --include-interactive
```

## Current Validation

On 2026-05-05, the full master runner completed with:

- 65 tasks
- 64 successful tasks
- 0 required failures
- 1 optional failure: `R/modularity/Bonferroni_modularity_integration.R`

The optional Bonferroni script expects a pre-existing `Outputs/Clean_modularity_analysis` directory. The runner leaves this as optional because the upstream modularity test scripts generate timestamped raw/residual outputs instead.
