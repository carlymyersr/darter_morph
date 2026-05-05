# ============================================================
# R/digitize_trait_landmarks.R
# Digitize landmarks for mouth angle measurement
# ============================================================

suppressPackageStartupMessages({
  library(StereoMorph)
})

VERBOSE <- TRUE

img_dir     <- "photos"
shapes_dir  <- "trait_measurements/mouth_angle_shapes"
lm_ref_file <- "trait_measurements/mouth_angle_landmarks_ref.txt"

if (!dir.exists(img_dir)) stop("photos folder not found")

if (!dir.exists(shapes_dir)) {
  dir.create(shapes_dir, recursive = TRUE)
  if (VERBOSE) cat("Created:", shapes_dir, "\n")
}

if (!file.exists(lm_ref_file)) stop("Landmark file not found")

if (VERBOSE) {
  cat("Launching trait digitizer...\n")
}

digitizeImages(
  image.file    = img_dir,
  shapes.file   = shapes_dir,
  landmarks.ref = lm_ref_file
)

if (VERBOSE) cat("Done.\n")