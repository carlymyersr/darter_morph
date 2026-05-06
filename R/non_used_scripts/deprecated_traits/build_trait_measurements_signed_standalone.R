# ============================================================
# R/build_trait_measurements_signed.R
#
# Compute SIGNED mouth angle from separate trait landmarks and
# merge with existing morphometric metadata.
#
# Trait landmarks expected:
#   mouth_tip
#   mouth_base
#   body_axis_anterior     # cranial/head anterior reference point, e.g. nostril
#   body_axis_posterior    # cranial/head posterior reference point, e.g. operculum tip
#
# Outputs:
#   trait_measurements/mouth_angle_signed_raw.csv
#   trait_measurements/mouth_angle_signed_with_metadata.csv
#   trait_measurements/mouth_angle_signed_extremes.csv
#   trait_measurements/mouth_angle_signed_histogram.png
#
# Interpretation goal:
#   positive signed_mouth_angle_deg = more upturned
#   negative signed_mouth_angle_deg = more downturned
#
# IMPORTANT:
#   Image coordinate conventions can flip signs. After running,
#   inspect mouth_angle_signed_extremes.csv. If positive/negative
#   are reversed, change ANGLE_SIGN_MULTIPLIER from 1 to -1.
# ============================================================

# ---------------------------
# User setting: sign convention
# ---------------------------
ANGLE_SIGN_MULTIPLIER <- -1
# If positive angles look DOWNturned after inspecting extremes,
# change this to:
# ANGLE_SIGN_MULTIPLIER <- -1

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
  if (dir.exists(file.path(project_root, "trait_measurements"))) {
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
  library(purrr)
})

# ---------------------------
# Inputs / outputs
# ---------------------------
trait_shapes_dir <- "trait_measurements/mouth_angle_shapes"
out_dir <- "trait_measurements"

raw_out       <- file.path(out_dir, "mouth_angle_signed_raw.csv")
merged_out    <- file.path(out_dir, "mouth_angle_signed_with_metadata.csv")
extremes_out  <- file.path(out_dir, "mouth_angle_signed_extremes.csv")
missing_out   <- file.path(out_dir, "mouth_angle_signed_missing_landmark_check.csv")
hist_out      <- file.path(out_dir, "mouth_angle_signed_histogram.png")

if (!dir.exists(trait_shapes_dir)) {
  stop(
    "Trait shapes folder not found: ", trait_shapes_dir,
    "\nRun R/10_digitize_trait_landmarks.R first."
  )
}

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ---------------------------
# Load main morphometric metadata
# ---------------------------
# This is only needed to build the merged metadata file.
# The later analysis script reads the merged CSV and does not rerun GPA/curves.
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
# Diagnostics: missing landmarks by specimen
# ---------------------------
missing_trait_check <- purrr::map_dfr(ids, function(id) {
  vals <- coords[required_lms, , id, drop = FALSE]
  missing_by_lm <- apply(vals, 1, function(z) any(is.na(z)))

  tibble(
    specimen = id,
    has_missing_trait_landmarks = any(missing_by_lm),
    missing_landmarks = paste(required_lms[missing_by_lm], collapse = "; ")
  )
})

write.csv(missing_trait_check, missing_out, row.names = FALSE)

if (any(missing_trait_check$has_missing_trait_landmarks)) {
  cat("\nSpecimens with missing trait landmarks:\n")
  print(missing_trait_check %>% filter(has_missing_trait_landmarks))
}

# ---------------------------
# Helper functions
# ---------------------------
get_point <- function(coords_arr, landmark, specimen) {
  coords_arr[landmark, , specimen]
}

signed_angle_degrees <- function(reference_vec, target_vec) {
  # Signed angle from reference_vec to target_vec in degrees.
  # Uses atan2(cross, dot), preserving direction.

  if (any(is.na(reference_vec)) || any(is.na(target_vec))) return(NA_real_)

  ref_norm <- sqrt(sum(reference_vec^2))
  tar_norm <- sqrt(sum(target_vec^2))

  if (is.na(ref_norm) || is.na(tar_norm) || ref_norm == 0 || tar_norm == 0) {
    return(NA_real_)
  }

  dot <- sum(reference_vec * target_vec)

  # 2D scalar cross product
  cross <- reference_vec[1] * target_vec[2] - reference_vec[2] * target_vec[1]

  angle_rad <- atan2(cross, dot)
  angle_rad * 180 / pi
}

angle_magnitude_degrees <- function(reference_vec, target_vec) {
  # Unsigned angle magnitude, useful as a diagnostic.

  if (any(is.na(reference_vec)) || any(is.na(target_vec))) return(NA_real_)

  denom <- sqrt(sum(reference_vec^2)) * sqrt(sum(target_vec^2))
  if (is.na(denom) || denom == 0) return(NA_real_)

  cos_theta <- sum(reference_vec * target_vec) / denom
  cos_theta <- max(min(cos_theta, 1), -1)

  acos(cos_theta) * 180 / pi
}

# ---------------------------
# Compute signed mouth angle
# ---------------------------
angle_df <- purrr::map_dfr(ids, function(id) {

  mouth_tip  <- get_point(coords, "mouth_tip", id)
  mouth_base <- get_point(coords, "mouth_base", id)

  axis_ant  <- get_point(coords, "body_axis_anterior", id)
  axis_post <- get_point(coords, "body_axis_posterior", id)

  # Mouth vector:
  # posterior/base of upper jaw -> anterior tip of upper jaw
  mouth_vec <- mouth_tip - mouth_base

  # Cranial/head reference vector:
  # posterior operculum/reference point -> anterior nostril/reference point
  body_vec <- axis_ant - axis_post

  raw_signed_angle <- signed_angle_degrees(
    reference_vec = body_vec,
    target_vec    = mouth_vec
  )

  signed_angle <- ANGLE_SIGN_MULTIPLIER * raw_signed_angle

  unsigned_angle <- angle_magnitude_degrees(
    reference_vec = body_vec,
    target_vec    = mouth_vec
  )

  tibble(
    specimen = id,
    raw_signed_mouth_angle_deg = raw_signed_angle,
    signed_mouth_angle_deg = signed_angle,
    mouth_angle_abs_deg = unsigned_angle,
    angle_sign_multiplier = ANGLE_SIGN_MULTIPLIER
  )
})

# Backward-compatible alias for older scripts if needed.
# NOTE: This is now SIGNED.
angle_df$mouth_angle_deg <- angle_df$signed_mouth_angle_deg

# ---------------------------
# Save raw signed mouth angle file
# ---------------------------
write.csv(angle_df, raw_out, row.names = FALSE)

cat("Wrote raw signed mouth angle data to:\n", raw_out, "\n")

# ---------------------------
# Merge with existing metadata
# ---------------------------
gdf_with_mouth_angle <- gdf %>%
  left_join(angle_df, by = "specimen")

write.csv(gdf_with_mouth_angle, merged_out, row.names = FALSE)

cat("Wrote merged metadata + signed mouth angle data to:\n", merged_out, "\n")

# ---------------------------
# Extreme specimens for visual sign check
# ---------------------------
extremes_df <- bind_rows(
  gdf_with_mouth_angle %>%
    filter(!is.na(signed_mouth_angle_deg)) %>%
    arrange(signed_mouth_angle_deg) %>%
    slice_head(n = 10) %>%
    mutate(extreme_type = "most_negative"),

  gdf_with_mouth_angle %>%
    filter(!is.na(signed_mouth_angle_deg)) %>%
    arrange(desc(signed_mouth_angle_deg)) %>%
    slice_head(n = 10) %>%
    mutate(extreme_type = "most_positive")
) %>%
  select(
    extreme_type,
    specimen,
    habitat,
    year,
    signed_mouth_angle_deg,
    raw_signed_mouth_angle_deg,
    mouth_angle_abs_deg,
    angle_sign_multiplier
  )

write.csv(extremes_df, extremes_out, row.names = FALSE)

cat("\nWrote extreme specimens for visual sign check to:\n", extremes_out, "\n")
cat("\nInspect these specimens. Desired convention:\n")
cat("  positive signed_mouth_angle_deg = more upturned\n")
cat("  negative signed_mouth_angle_deg = more downturned\n")
cat("If reversed, set ANGLE_SIGN_MULTIPLIER <- -1 and rerun.\n")

# ---------------------------
# Quick diagnostics
# ---------------------------
cat("\nSigned mouth angle summary:\n")
print(summary(gdf_with_mouth_angle$signed_mouth_angle_deg))

cat("\nUnsigned mouth angle magnitude summary:\n")
print(summary(gdf_with_mouth_angle$mouth_angle_abs_deg))

cat("\nMissing signed mouth angle values:\n")
print(table(is.na(gdf_with_mouth_angle$signed_mouth_angle_deg)))

# ---------------------------
# Diagnostic histogram with zero line
# ---------------------------
p_hist <- ggplot(
  gdf_with_mouth_angle,
  aes(x = signed_mouth_angle_deg)
) +
  geom_histogram(bins = 30, color = "black", fill = "gray80") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(
    x = "Signed mouth angle relative to cranial axis (degrees)",
    y = "Count",
    title = "Distribution of signed mouth angle measurements",
    subtitle = "Positive should be upturned; negative should be downturned after sign check"
  )

ggsave(
  hist_out,
  p_hist,
  width = 5.5,
  height = 4,
  dpi = 400
)

cat("\nSaved diagnostic histogram to:\n", hist_out, "\n")
cat("\nDone.\n")
