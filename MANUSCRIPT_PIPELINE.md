# Manuscript Pipeline

This file maps the methods/results outline in `~/Documents/darter_figures/outline methods and results.txt` to the canonical scripts in this repository.

## Running

Static validation only:

```sh
Rscript R/run_manuscript_pipeline.R --check-only
```

Full run:

```sh
Rscript R/run_manuscript_pipeline.R
```

Stop at the first required task failure:

```sh
Rscript R/run_manuscript_pipeline.R --stop-on-error
```

Archive/exploratory tasks are excluded unless explicitly requested:

```sh
Rscript R/run_manuscript_pipeline.R --include-archive
```

## Canonical Script Map

### Methods 1-2: Sampling, GPA, and Subsets

- `R/00_setup_morpho.R`
- `R/01_build_metadata.R`
- `R/02_subset_1950.R`
- `R/03_subset_CT_timeseries.R`
- `R/04_subset_CT_timeseries_plus_1950habitats.R`
- `R/05_subset_1950_quabbin_swift.R`

### Methods 3: Size Correction and Allometry

- `R/01_methods/01_size_allometry/02_size_boxplots_and_tests.R`
- `R/01_methods/01_size_allometry/03_1950_raw_pca_size_colored.R`
- `R/01_methods/01_size_allometry/04_ct_timeseries_raw_pca_size_colored.R`
- `R/01_methods/01_size_allometry/01_combined_raw_pca_size_colored.R`

### Methods 4: Genital Papillae Sensitivity

- `R/01_methods/02_papillae_sensitivity/01_papillae_parallel_sensitivity.R`

### Results 1: Structure of Variation in Morphospace

- `R/02_results/01_morphospace_structure/01_1950_waterbody_pca_residual_pc12.R`
- `R/02_results/01_morphospace_structure/02_ct_timeseries_pca_residual_pc12.R`
- `R/02_results/01_morphospace_structure/03_1950_tps_residual_pc12_extremes.R`
- `R/02_results/01_morphospace_structure/04_ct_timeseries_tps_residual_pc12_extremes.R`

Supplemental PC-axis views are kept in the same folder with `supplemental` in the filename.

### Results 2: Mean Shape Differentiation

- `R/02_results/02_mean_shape_differentiation/01_mean_positions_and_tps_by_group.R`
- `R/02_results/02_mean_shape_differentiation/02_1950_waterbody_procD_pairwise.R`
- `R/02_results/02_mean_shape_differentiation/03_ct_timeseries_procD_pairwise.R`

### Results 3: Trait-Specific Patterns

- `R/02_results/03_trait_patterns/01_linear_trait_measurements.R`
- `R/02_results/03_trait_patterns/02_curve_trait_measurements.R`
- `R/02_results/03_trait_patterns/03_mouth_angle_measurements_and_tests.R`
- `R/02_results/03_trait_patterns/04_trait_boxplots_combined_groups.R`
- `R/02_results/03_trait_patterns/05_mouth_angle_network_figure.R`
- `R/02_results/03_trait_patterns/06_trait_faceted_summary_figure.R`

The outline says to keep hyoid deviation and remove tortuosity/snout-deviation interpretation from the main text. Older broader trait-PC scripts are under `R/03_supplemental/trait_pc_exploratory/` and are excluded by default.

### Results 4: Distribution of Within-Group Variation

- `R/02_results/04_within_group_variation/01_disparity_and_dispersion_tests.R`

This is the canonical script for Procrustes variance, pairwise disparity, and `betadisper`. The older no-letter version is archived.

### Results 5: CT Reference Distribution

- `R/02_results/05_ct_reference_distribution/01_full_landscape_residual_pc12.R`
- `R/02_results/05_ct_reference_distribution/02_full_landscape_residual_pc23.R`

Mahalanobis thresholding is currently implemented in the hydrology script because the existing code defines the pooled mainstem reference there. Rarefaction/balanced subsampling is not yet present as a separate canonical script.

### Results 6: Persistent Local Divergence

- `R/02_results/06_persistent_local_divergence/01_quabbin_swift_context_dependence_procD.R`

### Results 7: Internal Reservoir Structure

- `R/02_results/07_reservoir_internal_structure/01_quabbin_swift_sampling_location_pca.R`

### Results 8: Hydrology-Based Structure

- `R/02_results/08_hydrology_structure/01_hydrology_groups_pca_procD_mahalanobis.R`

This script uses `variation_source_alt` for Mainstem, Reservoir System, and Tributaries.

### Results 9: Modularity and Integration

- `R/02_results/09_modularity_integration/TestA_5modules.R`
- `R/02_results/09_modularity_integration/TestB_4modules.R`
- `R/02_results/09_modularity_integration/TestC_3modules.R`
- `R/02_results/09_modularity_integration/TestD_2modules_AP.R`
- `R/02_results/09_modularity_integration/TestE_2modules_DV.R`
- CT versions of the same tests are suffixed `_CT3.R`.

The old multiple-testing correction script is retained as `R/02_results/09_modularity_integration/99_multiple_testing_correction_optional.R`, but it is excluded from the default master run because it expects a manually curated `Outputs/Clean_modularity_analysis` directory rather than the generated modularity output folders.

## Excluded From Main Pipeline

- `R/04_archive/sicb/` contains presentation-specific figure scripts.
- `R/04_archive/deprecated_traits/` contains unsigned or standalone trait variants replaced by the canonical trait section.
- `R/04_archive/deprecated_stats/` contains the older disparity script without compact letters.
- `R/03_supplemental/trait_pc_exploratory/` is excluded by default because the outline says those analyses are not included.

## Known Gap

The outline asks for rarefaction or balanced subsampling for the CT reference distribution. I did not find a dedicated existing script for that analysis. The current master runner documents this gap instead of silently substituting a different analysis.

The outline also asks for Bonferroni/BH correction across modularity and integration tests. The existing correction script is preserved, but it needs either a curated clean-output directory or a refactor to consume the timestamped output folders generated by the modularity test scripts.
