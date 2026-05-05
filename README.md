# Darter Morph

R scripts for the darter morphometrics pipeline, figure reproduction, trait analyses, modularity analyses, and supplemental analyses.

## Contents

- `R/00_setup_morpho.R` - core morphometric initialization
- `R/01_build_metadata.R` - specimen metadata aligned to GPA coordinates
- `R/01_helpers_angles.R` - helper functions for angle measurements
- `R/01_build_angle_measurements.R` - angle measurement construction
- `R/02_subset_1950.R` - 1950 specimen subset
- `R/03_subset_CT_timeseries.R` - CT time-series subset
- `R/04_subset_CT_timeseries_plus_1950habitats.R` - CT time-series plus 1950 habitats subset
- `R/05_subset_1950_quabbin_swift.R` - Quabbin and Swift 1950 subset
- `R/master_reproduce_figures.R` - master reproduction runner
- `R/figures/` - main figure scripts
- `R/traits/` - trait quantification and trait figure scripts
- `R/stats/` - model/statistical analysis scripts
- `R/modularity/` - modularity and integration scripts
- `R/supplemental/` - supplemental analyses

## Notes

The repository includes R scripts imported from `Documents/darter_pipeline/R` and `Documents/darter_figures`.
Input data, photos, generated figures, and analysis outputs are not committed.

See `REPRODUCING.md` for the master figure reproduction workflow.
