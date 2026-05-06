# Darter Morph

R scripts for the darter morphometrics manuscript pipeline.

The repository is organized around the current manuscript outline:

```text
~/Documents/darter_figures/outline methods and results 12.19.16 PM.txt
```

## Run

```sh
Rscript R/run_manuscript_pipeline.R
```

Static path validation only:

```sh
Rscript R/run_manuscript_pipeline.R --check-only
```

## Where Scripts Live

- `R/methods/01_specimen_sampling_study_design/`
- `R/methods/02_landmark_acquisition_geometric_morphometrics/`
- `R/methods/03_size_correction_allometry/`
- `R/methods/04_genital_papillae_sensitivity/`
- `R/results/01_structure_variation_morphospace/`
- `R/results/02_mean_shape_differentiation/`
- `R/results/03_trait_specific_patterns/`
- `R/results/04_within_group_variation/`
- `R/results/05_ct_reference_distribution/`
- `R/results/06_persistent_local_divergence/`
- `R/results/07_hydrology_based_structure/`
- `R/results/08_reservoir_internal_structure/`
- `R/results/09_modularity_integration/`
- `R/supplemental/` for supporting analyses retained outside the main results sequence.
- `R/non_used_scripts/` for old wrappers, deprecated scripts, SICB presentation scripts, and analyses excluded by the current outline.

## Where Outputs Go

Legacy scripts still write to their original working folders such as `figures/`, `Outputs/`, `papillae/`, and `outputs_mouth_to_body_angle/`.

The master runner now also copies newly created or updated outputs into a clean sectioned folder for each run:

```text
manuscript_outputs/<run_id>/<outline_section>/<figures|tables|models|other>/
```

Run logs are written separately to:

```text
Outputs/manuscript_runs/<run_id>/
```

## Data

By default, scripts expect local input data under:

```text
~/Documents/darter_morphometrics
```

Override with:

```sh
DARTER_DATA_ROOT=/path/to/darter_morphometrics Rscript R/run_manuscript_pipeline.R
```
