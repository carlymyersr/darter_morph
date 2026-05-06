# ============================================================
# R/build_trait_measurements.R
#
# Compute mouth angle from separate trait landmarks and merge
# with existing morphometric metadata.
#
# Trait landmarks expected:
#   mouth_tip
#   mouth_base
#   body_axis_anterior
#   body_axis_posterior
#
# Outputs:
#   trait_measurements/mouth_angle_raw.csv
#   trait_measurements/mouth_angle_with_metadata.csv
# ============================================================

# ---------------------------
# Root detection
# ---------------------------
get_script_dir <- function() {
  this <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (is.null(this) || is.na(this)) return(NULL)
  dirname(this)
}

script_dir <- get_script_dir()
if (!is.null(script_dir)) {
  project_root <- normalizePath(file.path(script_dir, ".."))
  if (file.exists(file.path(project_root, "darter_curves.txt"))) {
    if (getwd() != project_root) setwd(project_root)
    cat("Project root set to:", project_root, "\n")
  }
}

# ---------------------------
# Libraries
# ---------------------------
suppressPackageStartupMessages({
  library(StereoMorph)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

# ---------------------------
# Inputs / outputs
# ---------------------------
trait_shapes_dir <- "trait_measurements/mouth_angle_shapes"
out_dir <- "trait_measurements"

raw_out <- file.path(out_dir, "mouth_angle_raw.csv")
merged_out <- file.path(out_dir, "mouth_angle_with_metadata.csv")

if (!dir.exists(trait_shapes_dir)) {
  stop("Trait shapes folder not found: ", trait_shapes_dir,
       "\nRun R/10_digitize_trait_landmarks.R first.")
}

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ---------------------------
# Load main morphometric metadata
# ---------------------------
source("R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
source("R/methods/01_specimen_sampling_study_design/01_build_metadata.R")

# ---------------------------
# Read trait landmark shapes
# ---------------------------
trait_shapes <- StereoMorph::readShapes(trait_shapes_dir)

if (is.null(trait_shapes$landmarks.scaled)) {
  stop("No scaled landmarks found in: ", trait_shapes_dir)
}

coords <- trait_shapes$landmarks.scaled

required_lms <- c(
  "mouth_tip",
  "mouth_base",
  "body_axis_anterior",
  "body_axis_posterior"
)

missing_lms <- setdiff(required_lms, dimnames(coords)[[1]])

if (length(missing_lms) > 0) {
  stop(
    "Missing required trait landmarks: ",
    paste(missing_lms, collapse = ", ")
  )
}

ids <- dimnames(coords)[[3]]

# ---------------------------
# Helper functions
# ---------------------------
get_point <- function(coords_arr, landmark, specimen) {
  coords_arr[landmark, , specimen]
}

angle_between_degrees <- function(v1, v2) {
  denom <- sqrt(sum(v1^2)) * sqrt(sum(v2^2))
  
  if (denom == 0) return(NA_real_)
  
  cos_theta <- sum(v1 * v2) / denom
  
  # numerical protection
  cos_theta <- max(min(cos_theta, 1), -1)
  
  theta_rad <- acos(cos_theta)
  theta_rad * 180 / pi
}

# ---------------------------
# Compute mouth angle
# ---------------------------
angle_df <- purrr::map_dfr(ids, function(id) {
  
  mouth_tip  <- get_point(coords, "mouth_tip", id)
  mouth_base <- get_point(coords, "mouth_base", id)
  
  axis_ant  <- get_point(coords, "body_axis_anterior", id)
  axis_post <- get_point(coords, "body_axis_posterior", id)
  
  # Vectors
  # Mouth vector = posterior/base of mouth toward mouth tip
  mouth_vec <- mouth_tip - mouth_base
  
  # Body vector = posterior body-axis point toward anterior body-axis point
  body_vec <- axis_ant - axis_post
  
  mouth_angle_deg <- angle_between_degrees(mouth_vec, body_vec)
  
  tibble(
    specimen = id,
    mouth_angle_deg = mouth_angle_deg
  )
})

# ---------------------------
# Save raw mouth angle file
# ---------------------------
write.csv(angle_df, raw_out, row.names = FALSE)

cat("Wrote raw mouth angle data to:\n", raw_out, "\n")

# ---------------------------
# Merge with existing metadata
# ---------------------------
gdf_with_mouth_angle <- gdf %>%
  left_join(angle_df, by = "specimen")

write.csv(gdf_with_mouth_angle, merged_out, row.names = FALSE)

cat("Wrote merged metadata + mouth angle data to:\n", merged_out, "\n")

# ---------------------------
# Quick diagnostics
# ---------------------------
cat("\nMouth angle summary:\n")
print(summary(gdf_with_mouth_angle$mouth_angle_deg))

cat("\nMissing mouth angle values:\n")
print(table(is.na(gdf_with_mouth_angle$mouth_angle_deg)))

# ---------------------------
# Optional diagnostic histogram
# ---------------------------
p_hist <- ggplot(
  gdf_with_mouth_angle,
  aes(x = mouth_angle_deg)
) +
  geom_histogram(bins = 30, color = "black", fill = "gray80") +
  theme_classic() +
  labs(
    x = "Mouth angle relative to body axis (degrees)",
    y = "Count",
    title = "Distribution of mouth angle measurements"
  )

ggsave(
  file.path(out_dir, "mouth_angle_histogram.png"),
  p_hist,
  width = 5,
  height = 4,
  dpi = 400
)

cat("\nSaved diagnostic histogram to:\n",
    file.path(out_dir, "mouth_angle_histogram.png"), "\n")

cat("\nDone.\n")