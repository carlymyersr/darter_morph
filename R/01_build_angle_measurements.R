# ============================================================
# R/build_angle_measurements.R
#
# Compute angle-based trait measurements from separate trait landmarks
# and merge them with existing morphometric metadata.
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
#   trait_measurements/mouth_angle_histogram.png
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
  library(purrr)
})

# ---------------------------
# Source project setup + helpers
# ---------------------------
source("R/00_setup_morpho.R")
source("R/01_build_metadata.R")
source("R/helpers_angles.R")

# ---------------------------
# Inputs / outputs
# ---------------------------
trait_shapes_dir <- "trait_measurements/mouth_angle_shapes"
out_dir <- "trait_measurements"

raw_out    <- file.path(out_dir, "mouth_angle_raw.csv")
merged_out <- file.path(out_dir, "mouth_angle_with_metadata.csv")
hist_out   <- file.path(out_dir, "mouth_angle_histogram.png")

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
# Read trait landmark shapes
# ---------------------------
trait_shapes <- StereoMorph::readShapes(trait_shapes_dir)

if (is.null(trait_shapes$landmarks.scaled)) {
  stop("No scaled landmarks found in: ", trait_shapes_dir)
}

coords_angle <- trait_shapes$landmarks.scaled

required_lms <- c(
  "mouth_tip",
  "mouth_base",
  "body_axis_anterior",
  "body_axis_posterior"
)

missing_lms <- setdiff(required_lms, dimnames(coords_angle)[[1]])

if (length(missing_lms) > 0) {
  stop(
    "Missing required trait landmarks: ",
    paste(missing_lms, collapse = ", ")
  )
}

ids_angle <- dimnames(coords_angle)[[3]]

# ---------------------------
# Compute mouth angle
# ---------------------------
angle_df <- purrr::map_dfr(ids_angle, function(id) {
  tibble(
    specimen = id,
    mouth_angle_deg = compute_mouth_angle(coords_angle, id)
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
# Diagnostics
# ---------------------------
cat("\nMouth angle summary:\n")
print(summary(gdf_with_mouth_angle$mouth_angle_deg))

cat("\nMissing mouth angle values:\n")
print(table(is.na(gdf_with_mouth_angle$mouth_angle_deg)))

# ---------------------------
# Diagnostic histogram
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
  hist_out,
  p_hist,
  width = 5,
  height = 4,
  dpi = 400
)

cat("\nSaved diagnostic histogram to:\n", hist_out, "\n")
cat("\nDone.\n")
