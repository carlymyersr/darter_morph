# Darter Morph

R scripts for the darter morphometrics manuscript pipeline.

The repository is organized around the current manuscript methods/results outline rather than by the original figure folders in `Documents/darter_figures`.

## Main Entry Point

```sh
Rscript R/run_manuscript_pipeline.R
```

Compatibility entry point:

```sh
Rscript R/master_reproduce_figures.R
```

## R Directory Layout

- `R/00_setup_morpho.R`, `R/01_build_metadata.R`, `R/02_subset_1950.R`, `R/03_subset_CT_timeseries.R`, `R/04_subset_CT_timeseries_plus_1950habitats.R`, `R/05_subset_1950_quabbin_swift.R` - canonical setup and analysis subsets.
- `R/01_methods/` - method validation and sensitivity analyses.
- `R/02_results/` - canonical scripts for Results sections 1-9.
- `R/03_supplemental/` - supplemental analyses that support, but are not central to, the Results sequence.
- `R/04_archive/` - older, SICB, unsigned, duplicate, or exploratory scripts excluded from the main pipeline.

See `MANUSCRIPT_PIPELINE.md` for the section-by-section script map.

## Data

Input data, photos, and generated outputs are not committed. By default, scripts expect local data under:

```text
~/Documents/darter_morphometrics
```

Override with:

```sh
DARTER_DATA_ROOT=/path/to/darter_morphometrics Rscript R/run_manuscript_pipeline.R
```
