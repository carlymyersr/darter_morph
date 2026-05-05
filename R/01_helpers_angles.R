# ============================================================
# R/helpers_angles.R
#
# Helper functions for angle-based trait measurements.
# This file should define functions only. It should not read files,
# source setup scripts, write outputs, or make plots.
# ============================================================

# ---------------------------
# Extract one landmark point
# ---------------------------
get_point <- function(coords_arr, landmark, specimen) {
  coords_arr[landmark, , specimen]
}

# ---------------------------
# Angle between two vectors, in degrees
# ---------------------------
angle_between_degrees <- function(v1, v2) {
  denom <- sqrt(sum(v1^2)) * sqrt(sum(v2^2))

  if (denom == 0) return(NA_real_)

  cos_theta <- sum(v1 * v2) / denom

  # Numerical protection
  cos_theta <- max(min(cos_theta, 1), -1)

  theta_rad <- acos(cos_theta)
  theta_rad * 180 / pi
}

# ---------------------------
# Compute mouth angle for one specimen
#
# Current definition:
#   mouth vector = mouth_base -> mouth_tip
#   reference vector = body_axis_posterior -> body_axis_anterior
#
# Interpretation depends on landmark geometry:
#   larger angle = mouth direction more horizontal relative to snout/body axis
#   smaller angle = mouth direction more aligned/upturned toward snout axis
# ---------------------------
compute_mouth_angle <- function(coords_arr, specimen) {

  mouth_tip  <- get_point(coords_arr, "mouth_tip", specimen)
  mouth_base <- get_point(coords_arr, "mouth_base", specimen)

  axis_ant  <- get_point(coords_arr, "body_axis_anterior", specimen)
  axis_post <- get_point(coords_arr, "body_axis_posterior", specimen)

  mouth_vec <- mouth_tip - mouth_base
  body_vec  <- axis_ant - axis_post

  angle_between_degrees(mouth_vec, body_vec)
}
