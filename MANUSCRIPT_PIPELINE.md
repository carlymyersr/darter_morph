# Manuscript Pipeline

This file maps the current methods/results outline to the executable scripts in this repository.

Outline used:

```text
~/Documents/darter_figures/outline methods and results 12.19.16 PM.txt
```

## Running

```sh
Rscript R/run_manuscript_pipeline.R
```

For validation without running analyses:

```sh
Rscript R/run_manuscript_pipeline.R --check-only
```

The master script writes two kinds of run records:

- `Outputs/manuscript_runs/<run_id>/` contains logs, task plans, and manifests.
- `manuscript_outputs/<run_id>/` contains a clean section-specific copy of outputs generated or updated during that run.

## Script Map

### Methods 1: Specimen Sampling and Study Design

- `R/methods/01_specimen_sampling_study_design/01_build_metadata.R`
- `R/methods/01_specimen_sampling_study_design/02_subset_1950.R`
- `R/methods/01_specimen_sampling_study_design/03_subset_CT_timeseries.R`
- `R/methods/01_specimen_sampling_study_design/04_subset_CT_timeseries_plus_1950habitats.R`
- `R/methods/01_specimen_sampling_study_design/05_subset_1950_quabbin_swift.R`

### Methods 2: Landmark Acquisition and Geometric Morphometrics

- `R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R`

### Methods 3: Size Correction and Allometry

- `R/methods/03_size_correction_allometry/01_combined_raw_pca_size_colored.R`
- `R/methods/03_size_correction_allometry/02_size_boxplots_and_tests.R`
- `R/methods/03_size_correction_allometry/03_1950_raw_pca_size_colored.R`
- `R/methods/03_size_correction_allometry/04_ct_timeseries_raw_pca_size_colored.R`
- `R/methods/03_size_correction_allometry/05_combined_allometry_procD.R`

### Methods 4: Genital Papillae Sensitivity

- `R/methods/04_genital_papillae_sensitivity/01_papillae_parallel_sensitivity.R`

### Results 1: Structure of Variation in Morphospace

- `R/results/01_structure_variation_morphospace/01_1950_waterbody_pca_residual_pc12.R`
- `R/results/01_structure_variation_morphospace/02_ct_timeseries_pca_residual_pc12.R`
- `R/results/01_structure_variation_morphospace/03_1950_tps_residual_pc12_extremes.R`
- `R/results/01_structure_variation_morphospace/04_ct_timeseries_tps_residual_pc12_extremes.R`

### Results 2: Mean Shape Differentiation

- `R/results/02_mean_shape_differentiation/01_mean_positions_and_tps_by_group.R`
- `R/results/02_mean_shape_differentiation/02_1950_waterbody_procD_pairwise.R`
- `R/results/02_mean_shape_differentiation/03_ct_timeseries_procD_pairwise.R`

### Results 3: Trait-Specific Patterns

- `R/results/03_trait_specific_patterns/01_linear_trait_measurements.R`
- `R/results/03_trait_specific_patterns/02_curve_trait_measurements.R`
- `R/results/03_trait_specific_patterns/03_mouth_angle_measurements_and_tests.R`
- `R/results/03_trait_specific_patterns/04_trait_boxplots_combined_groups.R`
- `R/results/03_trait_specific_patterns/05_mouth_angle_network_figure.R`
- `R/results/03_trait_specific_patterns/06_trait_faceted_summary_figure.R`

### Results 4: Distribution of Within-Group Variation

- `R/results/04_within_group_variation/01_disparity_and_dispersion_tests.R`

### Results 5: Shared Empirical Landscape / CT Reference Distribution

- `R/results/05_ct_reference_distribution/01_full_landscape_residual_pc12.R`
- `R/results/05_ct_reference_distribution/02_full_landscape_residual_pc23.R`
- `R/results/05_ct_reference_distribution/03_ct_reference_mahalanobis_rarefaction.R`

### Results 6: Persistent Local Divergence

- `R/results/06_persistent_local_divergence/01_quabbin_swift_context_dependence_procD.R`

### Results 7: Hydrology-Based Structure

- `R/results/07_hydrology_based_structure/01_hydrology_groups_pca_procD_mahalanobis.R`

This script uses `variation_source_alt` for Mainstem, Reservoir System, and Tributaries.

### Results 8: Internal Reservoir Structure

- `R/results/08_reservoir_internal_structure/01_quabbin_swift_sampling_location_pca.R`

### Results 9: Modularity and Integration

- `R/results/09_modularity_integration/TestA_5modules.R`
- `R/results/09_modularity_integration/TestB_4modules.R`
- `R/results/09_modularity_integration/TestC_3modules.R`
- `R/results/09_modularity_integration/TestD_2modules_AP.R`
- `R/results/09_modularity_integration/TestE_2modules_DV.R`
- CT versions of the same tests are suffixed `_CT3.R`.
- `R/results/09_modularity_integration/99_modularity_multiple_testing_correction.R`

## Supplemental Scripts

- `R/supplemental/01_structure_variation_extra_axes/` contains PC2-PC3 and PC3-PC4 views and extra TPS outputs.
- `R/supplemental/03_size_correction_allometry/` contains supporting size/allometry checks.
- `R/supplemental/06_persistent_local_divergence/` contains the clean Swift/Quabbin PCA figure.
- `R/supplemental/07_hydrology_alternatives/` contains 1950-only and alternative hydrology grouping checks.

## Non-Used Scripts

Scripts excluded by the current outline are in:

```text
R/non_used_scripts/
```

This includes old compatibility wrappers, deprecated trait/stat scripts, SICB presentation scripts, and exploratory trait-PC analyses.
