# ============================================================
# Scripts/facet_mouth_angle_SICB_figure.R
#
# SICB-style faceted mouth angle figure.
#
# Goal:
#   Match the formatting/sizing of the SICB curve plots where
#   snout + hyoid are combined via ggplot faceting.
#
# Here:
#   - One facet panel contains signed mouth angle.
#   - One facet panel is intentionally empty, preserving the
#     same two-column facet layout used for Snout/Hyoid plots.
#
# Group order:
#   CT 1970, CT 1956, CT 1950, Quabbin, Swift, Fort, Sawmill
#
# Input:
#   trait_measurements/mouth_angle_signed_with_metadata.csv
#
# Before running:
#   source("R/build_trait_measurements_signed.R")
#
# Output:
#   Figures/SICB_signed_mouth_angle_facet/
#     facet_mouth_angle_SICB_figure.pdf
#     facet_mouth_angle_SICB_figure.png
#
# Interpretation:
#   signed_mouth_angle_deg should be oriented so that:
#     positive = more upturned
#     negative = more downturned
#   If reversed, change ANGLE_SIGN_MULTIPLIER in
#   R/build_trait_measurements_signed.R and rerun first.
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
  if (dir.exists(file.path(project_root, "trait_measurements"))) {
    if (getwd() != project_root) setwd(project_root)
    cat("Project root set to:", project_root, "\n")
  }
}

# ---------------------------
# Libraries
# ---------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ---------------------------
# Inputs / outputs
# ---------------------------
INFILE <- file.path("trait_measurements", "mouth_angle_signed_with_metadata.csv")

if (!file.exists(INFILE)) {
  stop(
    "Could not find: ", INFILE,
    "\nRun source('R/build_trait_measurements_signed.R') first."
  )
}

OUTDIR <- file.path("Figures", "SICB_signed_mouth_angle_facet")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Read signed mouth angle data
# ---------------------------
mouth_df <- read.csv(INFILE, stringsAsFactors = FALSE)

required_cols <- c("specimen", "habitat", "year", "signed_mouth_angle_deg")
missing_cols <- setdiff(required_cols, names(mouth_df))

if (length(missing_cols) > 0) {
  stop("Input file missing required columns: ", paste(missing_cols, collapse = ", "))
}

mouth_df <- mouth_df %>%
  mutate(
    year = as.integer(year),
    signed_mouth_angle_deg = as.numeric(signed_mouth_angle_deg)
  )

cat("\nSigned mouth angle missing check:\n")
print(table(is.na(mouth_df$signed_mouth_angle_deg), useNA = "ifany"))

# ============================================================
# Build combined dataset
# ============================================================

combined_levels <- c(
  "CT 1970",
  "CT 1956",
  "CT 1950",
  "Quabbin",
  "Swift",
  "Fort",
  "Sawmill"
)

mouth_plot_dat <- mouth_df %>%
  dplyr::filter(
    (
      habitat == "Connecticut River" & year %in% c(1950, 1956, 1970)
    ) |
      (
        year == 1950 & habitat %in% c(
          "Quabbin",
          "Swift River",
          "Fort River",
          "Sawmill River"
        )
      ),
    is.finite(signed_mouth_angle_deg)
  ) %>%
  dplyr::mutate(
    group = dplyr::case_when(
      habitat == "Connecticut River" & year == 1970 ~ "CT 1970",
      habitat == "Connecticut River" & year == 1956 ~ "CT 1956",
      habitat == "Connecticut River" & year == 1950 ~ "CT 1950",
      habitat == "Quabbin" ~ "Quabbin",
      habitat == "Swift River" ~ "Swift",
      habitat == "Fort River" ~ "Fort",
      habitat == "Sawmill River" ~ "Sawmill",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = combined_levels),
    facet_label = factor(
      "Mouth angle",
      levels = c("Mouth angle", " ")
    )
  ) %>%
  dplyr::filter(!is.na(group))

cat("\nGroup counts:\n")
print(table(mouth_plot_dat$group, useNA = "ifany"))

# Create an empty facet panel so the figure has the same 2-column
# visual structure as the Snout/Hyoid curve plots.
empty_facet_dat <- mouth_plot_dat %>%
  dplyr::slice(0) %>%
  dplyr::mutate(
    facet_label = factor(" ", levels = c("Mouth angle", " "))
  )

mouth_plot_dat_faceted <- dplyr::bind_rows(
  mouth_plot_dat,
  empty_facet_dat
)

# ============================================================
# SICB-style faceted plot
# ============================================================

p_mouth_angle_facet <- ggplot(
  mouth_plot_dat_faceted,
  aes(x = group, y = signed_mouth_angle_deg)
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.25,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.25
  ) +
  geom_jitter(
    width = 0.10,
    height = 0,
    alpha = 0.8,
    size = 0.8
  ) +
  facet_wrap(
    ~ facet_label,
    scales = "free_y",
    ncol = 2,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Signed mouth angle (degrees)"
  ) +
  theme_bw(base_size = 6, base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.15),
    panel.border = element_rect(linewidth = 0.25),
    strip.background = element_rect(fill = "grey95", linewidth = 0.25),
    strip.text = element_text(size = 6),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 5),
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.ticks = element_line(linewidth = 0.25),
    plot.margin = margin(2, 2, 2, 2, unit = "pt")
  )

# ---------------------------
# Save with exact combined curve plot sizing
# ---------------------------
# Matches:
#   curve_shape_metrics_combined_CT_timepoints_1950_habitats_no_tortuosity_SICB.pdf
# which is saved at width = 4, height = 2.6 inches.

ggsave(
  filename = file.path(OUTDIR, "facet_mouth_angle_SICB_figure.pdf"),
  plot = p_mouth_angle_facet,
  width = 4,
  height = 2.6,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = file.path(OUTDIR, "facet_mouth_angle_SICB_figure.png"),
  plot = p_mouth_angle_facet,
  width = 4,
  height = 2.6,
  units = "in",
  dpi = 400,
  bg = "white"
)

cat("\nSaved faceted mouth angle SICB figure to:\n")
cat(normalizePath(OUTDIR), "\n")
print(list.files(OUTDIR))

cat("\nDone.\n")
