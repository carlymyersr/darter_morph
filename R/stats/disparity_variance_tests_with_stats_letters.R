# ============================================================
# Scripts/disparity_variance_tests_with_stats_letters.R
#
# Run:
#   1) Morphological disparity analysis
#   2) Variance / dispersion test across groups
#
# Groups included:
#   - All 1950 habitats
#   - Connecticut River 1956
#   - Connecticut River 1970
#
# Uses:
#   - residual shape coordinates from
#     R/04_subset_CT_timeseries_plus_1950habitats.R
#
# Outputs:
#   Figures/disparity_boxplot_Fig6groups.png/.pdf
#   Figures/variance_boxplot_Fig6groups.png/.pdf
#   Outputs/disparity_variance_tests_<runid>/...
#
# Notes:
#   - "Disparity" here = Procrustes variance within group
#   - "Variance test" here = multivariate dispersion test based on
#     specimen distance to group centroid in tangent space
# disparity analysis on residualized shape data
#For the variance test, I interpreted that as a multivariate dispersion test: 
#are some groups more spread out around their own centroid than others? 
#That is different from mean-shape testing.

#The disparity boxplot is based on bootstrap distributions of within-group Procrustes variance, with the observed disparity overlaid as white points.
#The variance boxplot is based on individual distances to their group centroid.

#morphol.disparity() (Procrustes variance)
#betadisper() (distance to centroid)
#How wide or diffuse each group is around its own centroid.

#betadisper() is from vegan, not geomorph.
#It’s not part of the standard morphometrics workflow
#It comes from ecology (beta diversity logic)


#For each group:
#Find its centroid
#Measure how far each specimen is from that centroid
#Summarize that spread

#Then ask for example:
  
#  Is CT more spread out than tributaries?
#  Is Quabbin tighter than everything else?


#Mean shape analyses test for differences in group centroids in morphospace (procD.lm), 
#whereas disparity and dispersion analyses quantify the distribution of variation 
#within groups, allowing us to distinguish shifts in mean phenotype from changes in 
#the structure of phenotypic variation.
# ============================================================

# ---------------------------
# Root detection (source-safe)
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
  } else {
    warning("Detected project_root but darter_curves.txt not found there; leaving getwd() unchanged.")
  }
} else {
  message("Could not detect script directory (likely run interactively). Leaving getwd() unchanged.")
}

# ---------------------------
# Libraries
# ---------------------------
suppressPackageStartupMessages({
  library(geomorph)
  library(RRPP)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(stringr)
  library(vegan)
  library(multcompView)
})

# ---------------------------
# Parameters
# ---------------------------
set.seed(123)

run_id        <- format(Sys.time(), "%Y%m%d_%H%M%S")
fig_dir       <- "Figures"
out_root      <- file.path("Outputs", paste0("disparity_variance_tests_", run_id))
tables_dir    <- file.path(out_root, "tables")
text_dir      <- file.path(out_root, "text")
rds_dir       <- file.path(out_root, "rds")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(text_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Helper functions for disparity/dispersion stats + letters
# ============================================================

calc_group_disparity <- function(X, groups) {
  groups <- droplevels(factor(groups))
  levs <- levels(groups)
  purrr::map_dfr(levs, function(g) {
    Xg <- X[groups == g, , drop = FALSE]
    ctr <- colMeans(Xg)
    sqd <- rowSums((Xg - matrix(ctr, nrow = nrow(Xg), ncol = ncol(Xg), byrow = TRUE))^2)
    tibble(
      group = g,
      n = nrow(Xg),
      disparity = mean(sqd)
    )
  })
}

bootstrap_group_disparity <- function(Xg, nboot = 999) {
  n <- nrow(Xg)
  out <- numeric(nboot)
  for (b in seq_len(nboot)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    Xb  <- Xg[idx, , drop = FALSE]
    ctr <- colMeans(Xb)
    sqd <- rowSums((Xb - matrix(ctr, nrow = nrow(Xb), ncol = ncol(Xb), byrow = TRUE))^2)
    out[b] <- mean(sqd)
  }
  out
}

make_pairwise_perm_disparity <- function(coords2d, groups, nperm = 999, p_adjust = "BH") {
  groups <- droplevels(factor(groups))
  levs <- levels(groups)
  obs <- calc_group_disparity(coords2d, groups)
  
  if (length(levs) < 2) {
    return(tibble(group1 = character(), group2 = character(), obs_diff = numeric(), p_raw = numeric(), p_adj = numeric()))
  }
  
  pairwise <- combn(levs, 2, simplify = FALSE) %>%
    purrr::map_dfr(function(pair) {
      g1 <- pair[1]
      g2 <- pair[2]
      obs_diff <- abs(obs$disparity[obs$group == g1] - obs$disparity[obs$group == g2])
      
      keep <- groups %in% c(g1, g2)
      X_sub <- coords2d[keep, , drop = FALSE]
      g_sub <- droplevels(groups[keep])
      
      perm_diffs <- replicate(nperm, {
        shuffled <- sample(g_sub)
        vals <- sapply(c(g1, g2), function(g) {
          Xg <- X_sub[shuffled == g, , drop = FALSE]
          ctr <- colMeans(Xg)
          sqd <- rowSums((Xg - matrix(ctr, nrow = nrow(Xg), ncol = ncol(Xg), byrow = TRUE))^2)
          mean(sqd)
        })
        abs(diff(vals))
      })
      
      tibble(
        group1 = g1,
        group2 = g2,
        obs_diff = obs_diff,
        p_raw = (sum(perm_diffs >= obs_diff) + 1) / (nperm + 1)
      )
    }) %>%
    mutate(p_adj = p.adjust(p_raw, method = p_adjust))
  
  pairwise
}

make_pairwise_wilcox_table <- function(values, groups, p_adjust = "BH") {
  groups <- droplevels(factor(groups))
  pw <- pairwise.wilcox.test(
    x = values,
    g = groups,
    p.adjust.method = p_adjust,
    exact = FALSE
  )
  
  out <- as.data.frame(as.table(pw$p.value), stringsAsFactors = FALSE) %>%
    filter(!is.na(Freq)) %>%
    transmute(
      group1 = as.character(Var1),
      group2 = as.character(Var2),
      p_adj = as.numeric(Freq)
    )
  
  out
}

make_letters_from_pairwise <- function(pairwise_df, group_levels, p_col = "p_adj", alpha = 0.05) {
  group_levels <- as.character(group_levels)
  
  if (is.null(pairwise_df) || nrow(pairwise_df) == 0) {
    return(tibble(group = factor(group_levels, levels = group_levels), letters = "a"))
  }
  
  pvals <- pairwise_df[[p_col]]
  names(pvals) <- paste(pairwise_df$group1, pairwise_df$group2, sep = "-")
  
  letters <- multcompView::multcompLetters(
    pvals,
    threshold = alpha
  )$Letters
  
  missing_groups <- setdiff(group_levels, names(letters))
  if (length(missing_groups) > 0) {
    letters <- c(letters, stats::setNames(rep("a", length(missing_groups)), missing_groups))
  }
  
  tibble(
    group = factor(names(letters), levels = group_levels),
    letters = as.character(letters)
  ) %>%
    arrange(group)
}

make_letter_positions <- function(plot_df, group_col = "group", y_col, pad_frac = 0.08) {
  y_range <- range(plot_df[[y_col]], na.rm = TRUE)
  pad <- diff(y_range) * pad_frac
  if (!is.finite(pad) || pad == 0) pad <- max(abs(y_range), na.rm = TRUE) * 0.08
  if (!is.finite(pad) || pad == 0) pad <- 0.001
  
  plot_df %>%
    group_by(.data[[group_col]]) %>%
    summarize(y_pos = max(.data[[y_col]], na.rm = TRUE) + pad, .groups = "drop") %>%
    rename(group = all_of(group_col))
}

add_letters_to_plot <- function(p, letters_df, text_size = 4) {
  p + geom_text(
    data = letters_df,
    aes(x = group, y = y_pos, label = letters),
    inherit.aes = FALSE,
    size = text_size
  )
}

# ---------------------------
# Source pipeline
# ---------------------------
source("R/00_setup_morpho.R")
source("R/01_build_metadata.R")

# helper in case it is not already defined upstream
if (!exists("subset_coords_to_gdf", inherits = TRUE)) {
  subset_coords_to_gdf <- function(coords_arr, gdf_sub, specimen_col = "specimen") {
    ids <- dimnames(coords_arr)[[3]]
    keep <- gdf_sub[[specimen_col]]
    stopifnot(all(keep %in% ids))
    coords_sub <- coords_arr[, , match(keep, ids), drop = FALSE]
    stopifnot(identical(dimnames(coords_sub)[[3]], keep))
    coords_sub
  }
}

source("R/04_subset_CT_timeseries_plus_1950habitats.R")

# ---------------------------
# Analysis dataset
# ---------------------------
# Use allometry-corrected residual shapes from Fig6/Fig7 subset
coords_use <- coords_resid_Fig6
meta_use   <- gdf_Fig6

stopifnot(identical(dimnames(coords_use)[[3]], meta_use$specimen))
stopifnot(!any(is.na(meta_use$group)))

meta_use <- meta_use %>%
  mutate(
    group = droplevels(group)
  )

group_levels <- levels(meta_use$group)

# palette: keep your habitat palette logic, but make CT years explicit
group_palette <- c(
  "CT_1950"            = "steelblue",
  "CT_1956"            = "deepskyblue3",
  "CT_1970"            = "navy",
  "Fort River_1950"    = "black",
  "Quabbin_1950"       = "darkgoldenrod2",
  "Sawmill River_1950" = "orchid3",
  "Swift River_1950"   = "tomato"
)

group_palette <- group_palette[group_levels]

# ============================================================
# 1) MORPHOLOGICAL DISPARITY
# ============================================================

cat("\n============================================================\n")
cat("Running disparity analysis...\n")
cat("============================================================\n")

# geomorph disparity model
fit_disp <- geomorph::procD.lm(coords_use ~ group, data = meta_use, iter = 999, RRPP = TRUE)

disp_obj <- geomorph::morphol.disparity(
  fit_disp,
  groups = meta_use$group,
  iter = 999
)

# ------------------------------------------------------------
# Make group-level observed disparity table manually as well
# (Procrustes variance within each group = mean squared distance
# to the group mean shape)
# ------------------------------------------------------------
coords2d <- two.d.array(coords_use)  # n x variables
rownames(coords2d) <- meta_use$specimen

group_disparity_obs <- calc_group_disparity(coords2d, meta_use$group) %>%
  mutate(group = factor(group, levels = group_levels))

# ------------------------------------------------------------
# Bootstrap disparity distributions for boxplot
# ------------------------------------------------------------

disp_boot_df <- map_dfr(group_levels, function(g) {
  idx <- which(meta_use$group == g)
  Xg  <- coords2d[idx, , drop = FALSE]
  tibble(
    group = g,
    disparity = bootstrap_group_disparity(Xg, nboot = 999)
  )
}) %>%
  mutate(group = factor(group, levels = group_levels))

disp_summary <- disp_boot_df %>%
  group_by(group) %>%
  summarize(
    mean_boot = mean(disparity),
    sd_boot   = sd(disparity),
    q025      = quantile(disparity, 0.025),
    q25       = quantile(disparity, 0.25),
    median    = quantile(disparity, 0.50),
    q75       = quantile(disparity, 0.75),
    q975      = quantile(disparity, 0.975),
    .groups = "drop"
  ) %>%
  left_join(group_disparity_obs, by = "group") %>%
  select(group, n, disparity_observed = disparity, everything(), -disparity)

# ------------------------------------------------------------
# Disparity pairwise stats + letters
# Letters are based on BH-adjusted permutation tests of observed
# pairwise differences in Procrustes variance.
# ------------------------------------------------------------
disp_pairwise <- make_pairwise_perm_disparity(
  coords2d = coords2d,
  groups = meta_use$group,
  nperm = 999
)

disp_letters <- make_letters_from_pairwise(
  pairwise_df = disp_pairwise,
  group_levels = group_levels,
  p_col = "p_adj",
  alpha = 0.05
)

disp_letter_pos <- make_letter_positions(
  disp_boot_df,
  group_col = "group",
  y_col = "disparity"
) %>%
  left_join(disp_letters, by = "group")

# ------------------------------------------------------------
# Disparity figure
# ------------------------------------------------------------
p_disp <- ggplot(disp_boot_df, aes(x = group, y = disparity, fill = group)) +
  geom_boxplot(width = 0.72, outlier.shape = NA, alpha = 0.9) +
  geom_point(
    data = group_disparity_obs,
    aes(x = group, y = disparity),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.8,
    stroke = 0.8,
    fill = "white",
    color = "black"
  ) +
  geom_text(
    data = disp_letter_pos,
    aes(x = group, y = y_pos, label = letters),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_fill_manual(values = group_palette, drop = FALSE) +
  labs(
    title = "Morphological disparity by group",
    subtitle = "Boxplots = bootstrap distributions; white points = observed Procrustes variance; letters = BH-adjusted pairwise tests",
    x = NULL,
    y = "Procrustes variance (within-group disparity)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave(
  filename = file.path(fig_dir, "disparity_boxplot_Fig6groups.png"),
  plot = p_disp, width = 8.5, height = 5.5, dpi = 400
)
ggsave(
  filename = file.path(fig_dir, "disparity_boxplot_Fig6groups.pdf"),
  plot = p_disp, width = 8.5, height = 5.5
)

# ============================================================
# 2) VARIANCE / DISPERSION TEST
# ============================================================

cat("\n============================================================\n")
cat("Running variance / dispersion test...\n")
cat("============================================================\n")

# Tangent-space coordinates already represented in coords_use after GPA/residualization
# Distances among specimens:
D <- dist(coords2d)

# Multivariate dispersion by group
bd <- vegan::betadisper(D, group = meta_use$group, type = "centroid", bias.adjust = FALSE)

# Permutation test for group differences in dispersion
bd_perm <- vegan::permutest(bd, permutations = 999)

# Pairwise permutation tests
bd_pair <- vegan::permutest(bd, pairwise = TRUE, permutations = 999)

# Specimen-level distances to group centroid for boxplot
variance_df <- tibble(
  specimen = names(bd$distances),
  group = factor(meta_use$group[match(names(bd$distances), meta_use$specimen)], levels = group_levels),
  dist_to_centroid = as.numeric(bd$distances)
)

variance_summary <- variance_df %>%
  group_by(group) %>%
  summarize(
    n = n(),
    mean_dist = mean(dist_to_centroid),
    sd_dist   = sd(dist_to_centroid),
    median_dist = median(dist_to_centroid),
    min_dist = min(dist_to_centroid),
    max_dist = max(dist_to_centroid),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Dispersion pairwise stats + letters
# Letters are based on BH-adjusted pairwise Wilcoxon tests on
# individual distances to group centroid.
# ------------------------------------------------------------
variance_pairwise_wilcox <- make_pairwise_wilcox_table(
  values = variance_df$dist_to_centroid,
  groups = variance_df$group,
  p_adjust = "BH"
)

variance_letters <- make_letters_from_pairwise(
  pairwise_df = variance_pairwise_wilcox,
  group_levels = group_levels,
  p_col = "p_adj",
  alpha = 0.05
)

variance_letter_pos <- make_letter_positions(
  variance_df,
  group_col = "group",
  y_col = "dist_to_centroid"
) %>%
  left_join(variance_letters, by = "group")

# ------------------------------------------------------------
# Variance figure
# ------------------------------------------------------------
p_var <- ggplot(variance_df, aes(x = group, y = dist_to_centroid, fill = group)) +
  geom_boxplot(width = 0.72, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.75) +
  geom_text(
    data = variance_letter_pos,
    aes(x = group, y = y_pos, label = letters),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_fill_manual(values = group_palette, drop = FALSE) +
  labs(
    title = "Multivariate dispersion by group",
    subtitle = "Distances of individuals to their group centroid",
    x = NULL,
    y = "Distance to group centroid"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave(
  filename = file.path(fig_dir, "variance_boxplot_Fig6groups.png"),
  plot = p_var, width = 8.5, height = 5.5, dpi = 400
)
ggsave(
  filename = file.path(fig_dir, "variance_boxplot_Fig6groups.pdf"),
  plot = p_var, width = 8.5, height = 5.5
)

# ============================================================
# 3) WRITE OUTPUTS
# ============================================================

# ---------------------------
# Save tables
# ---------------------------
write.csv(
  group_disparity_obs,
  file.path(tables_dir, "disparity_observed_by_group.csv"),
  row.names = FALSE
)

write.csv(
  disp_boot_df,
  file.path(tables_dir, "disparity_bootstrap_distributions.csv"),
  row.names = FALSE
)

write.csv(
  disp_summary,
  file.path(tables_dir, "disparity_summary.csv"),
  row.names = FALSE
)

write.csv(
  disp_pairwise,
  file.path(tables_dir, "disparity_pairwise_permutation_tests_BH.csv"),
  row.names = FALSE
)

write.csv(
  disp_letters,
  file.path(tables_dir, "disparity_letters_BH.csv"),
  row.names = FALSE
)

write.csv(
  variance_df,
  file.path(tables_dir, "variance_distances_to_centroid.csv"),
  row.names = FALSE
)

write.csv(
  variance_summary,
  file.path(tables_dir, "variance_summary.csv"),
  row.names = FALSE
)

write.csv(
  variance_pairwise_wilcox,
  file.path(tables_dir, "variance_pairwise_wilcox_BH.csv"),
  row.names = FALSE
)

write.csv(
  variance_letters,
  file.path(tables_dir, "variance_letters_BH.csv"),
  row.names = FALSE
)

# pairwise table extraction if available
pairwise_df <- NULL
if (!is.null(bd_pair$pairwise)) {
  pw <- bd_pair$pairwise
  
  # Try to coerce the pairwise object into a readable long table
  if (is.matrix(pw)) {
    pw_df <- as.data.frame(as.table(pw), stringsAsFactors = FALSE)
    names(pw_df) <- c("group1", "group2", "value")
    pairwise_df <- pw_df %>%
      filter(!is.na(value))
  } else if (is.data.frame(pw)) {
    pairwise_df <- pw
  }
}

manifest <- bind_rows(
  manifest,
  tibble::tribble(
    ~category, ~file,
    "table", file.path(tables_dir, "disparity_pairwise_permutation_tests_BH.csv"),
    "table", file.path(tables_dir, "disparity_letters_BH.csv"),
    "table", file.path(tables_dir, "variance_pairwise_wilcox_BH.csv"),
    "table", file.path(tables_dir, "variance_letters_BH.csv"),
    "table", file.path(tables_dir, "CT_1950_vs_CT_ALL_disparity_pairwise_BH.csv"),
    "table", file.path(tables_dir, "CT_1950_vs_CT_ALL_disparity_letters_BH.csv"),
    "table", file.path(tables_dir, "CT_1950_vs_CT_ALL_dispersion_pairwise_wilcox_BH.csv"),
    "table", file.path(tables_dir, "CT_1950_vs_CT_ALL_dispersion_letters_BH.csv"),
    "table", file.path(tables_dir, "CT_timeseries_vs_1950habitats_disparity_pairwise_BH.csv"),
    "table", file.path(tables_dir, "CT_timeseries_vs_1950habitats_disparity_letters_BH.csv"),
    "table", file.path(tables_dir, "CT_timeseries_vs_1950habitats_dispersion_pairwise_wilcox_BH.csv"),
    "table", file.path(tables_dir, "CT_timeseries_vs_1950habitats_dispersion_letters_BH.csv"),
    "text", file.path(text_dir, "CT_1950_vs_CT_ALL_disparity_dispersion_stats.txt"),
    "text", file.path(text_dir, "CT_timeseries_vs_1950habitats_disparity_dispersion.txt")
  )
)

if (!is.null(pairwise_df)) {
  write.csv(
    pairwise_df,
    file.path(tables_dir, "variance_pairwise_permutation_tests.csv"),
    row.names = FALSE
  )
}

# ---------------------------
# Save R objects
# ---------------------------
saveRDS(fit_disp,  file.path(rds_dir, "fit_disp_procDlm.rds"))
saveRDS(disp_obj,  file.path(rds_dir, "disp_obj_morphol_disparity.rds"))
saveRDS(bd,        file.path(rds_dir, "betadisper_obj.rds"))
saveRDS(bd_perm,   file.path(rds_dir, "betadisper_permutest.rds"))
saveRDS(bd_pair,   file.path(rds_dir, "betadisper_pairwise.rds"))

# ---------------------------
# Write text summaries
# ---------------------------
sink(file.path(text_dir, "disparity_analysis.txt"))
cat("============================================================\n")
cat("DISPARITY ANALYSIS\n")
cat("============================================================\n\n")
cat("Dataset: all 1950 habitats + Connecticut River 1956 and 1970\n")
cat("Shape data: allometry-corrected residual GPA coordinates (coords_resid_Fig6)\n\n")

cat("Group counts:\n")
print(table(meta_use$group))
cat("\nObserved disparity by group:\n")
print(group_disparity_obs)
cat("\nBootstrap summary by group:\n")
print(disp_summary)

cat("\nprocD.lm summary:\n")
print(summary(fit_disp))

cat("\nmorphol.disparity object:\n")
print(disp_obj)
sink()

sink(file.path(text_dir, "variance_analysis.txt"))
cat("============================================================\n")
cat("VARIANCE / DISPERSION ANALYSIS\n")
cat("============================================================\n\n")
cat("Dataset: all 1950 habitats + Connecticut River 1956 and 1970\n")
cat("Shape data: allometry-corrected residual GPA coordinates (coords_resid_Fig6)\n\n")

cat("Group counts:\n")
print(table(meta_use$group))

cat("\nDistance-to-centroid summary by group:\n")
print(variance_summary)

cat("\nANOVA on multivariate dispersion:\n")
print(anova(bd))

cat("\nPermutation test on multivariate dispersion:\n")
print(bd_perm)

cat("\nPairwise permutation tests:\n")
print(bd_pair)
sink()

# ============================================================
# 4) CT ONLY COMPARISON: CT1950 vs CT_ALL
# ============================================================

cat("\n============================================================\n")
cat("Running CT-only disparity and dispersion comparisons...\n")
cat("============================================================\n")

ct_meta <- meta_use %>%
  filter(grepl("^CT_", group)) %>%
  mutate(
    ct_group = ifelse(group == "CT_1950", "CT_1950", "CT_ALL")
  )

ct_coords <- subset_coords_to_gdf(coords_use, ct_meta)

coords2d_ct <- two.d.array(ct_coords)
rownames(coords2d_ct) <- ct_meta$specimen

#observed disparity
ct_levels <- c("CT_1950", "CT_ALL")

ct_disparity_obs <- calc_group_disparity(coords2d_ct, ct_meta$ct_group) %>%
  mutate(group = factor(group, levels = ct_levels))

#bootstrap disparity
ct_disp_boot <- map_dfr(ct_levels, function(g) {
  idx <- which(ct_meta$ct_group == g)
  Xg  <- coords2d_ct[idx, , drop = FALSE]
  
  tibble(
    group = g,
    disparity = bootstrap_group_disparity(Xg, nboot = 999)
  )
}) %>%
  mutate(group = factor(group, levels = ct_levels))

# pairwise disparity stats + letters
ct_disp_pairwise <- make_pairwise_perm_disparity(
  coords2d = coords2d_ct,
  groups = ct_meta$ct_group,
  nperm = 999
)

ct_disp_letters <- make_letters_from_pairwise(
  pairwise_df = ct_disp_pairwise,
  group_levels = ct_levels,
  p_col = "p_adj",
  alpha = 0.05
)

ct_disp_letter_pos <- make_letter_positions(
  ct_disp_boot,
  group_col = "group",
  y_col = "disparity"
) %>%
  left_join(ct_disp_letters, by = "group")

#disparity boxplot
p_ct_disp <- ggplot(ct_disp_boot, aes(x = group, y = disparity, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_point(
    data = ct_disparity_obs,
    aes(x = group, y = disparity),
    inherit.aes = FALSE,
    shape = 21,
    size = 3,
    fill = "white"
  ) +
  geom_text(
    data = ct_disp_letter_pos,
    aes(x = group, y = y_pos, label = letters),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_fill_manual(values = c("CT_1950" = "steelblue", "CT_ALL" = "navy")) +
  theme_classic() +
  labs(
    title = "Disparity: CT 1950 vs CT (1950 + 1956 + 1970)",
    y = "Procrustes variance"
  )

ggsave(
  file.path(fig_dir, "CT_disparity_1950_vs_all.png"),
  p_ct_disp, width = 6, height = 5, dpi = 400
)

#permutation test for disparity differences
ct_disp_diff <- abs(diff(ct_disparity_obs$disparity))

perm_diff <- replicate(999, {
  shuffled <- sample(ct_meta$ct_group)
  
  vals <- sapply(ct_levels, function(g) {
    idx <- which(shuffled == g)
    Xg <- coords2d_ct[idx, , drop = FALSE]
    ctr <- colMeans(Xg)
    sqd <- rowSums((Xg - matrix(ctr, nrow = nrow(Xg), ncol = ncol(Xg), byrow = TRUE))^2)
    mean(sqd)
  })
  
  abs(diff(vals))
})

ct_disp_p_value <- (sum(perm_diff >= ct_disp_diff) + 1) / (length(perm_diff) + 1)

cat("\nCT disparity difference test p-value:", ct_disp_p_value, "\n")


#dispersion comparison 1950 vs all
D_ct <- dist(coords2d_ct)

bd_ct <- betadisper(D_ct, group = ct_meta$ct_group)

bd_ct_perm <- permutest(bd_ct, permutations = 999)

#build dataframe
ct_variance_df <- tibble(
  specimen = names(bd_ct$distances),
  group = factor(ct_meta$ct_group[match(names(bd_ct$distances), ct_meta$specimen)], levels = ct_levels),
  dist_to_centroid = as.numeric(bd_ct$distances)
)

# pairwise dispersion stats + letters
ct_var_pairwise_wilcox <- make_pairwise_wilcox_table(
  values = ct_variance_df$dist_to_centroid,
  groups = ct_variance_df$group,
  p_adjust = "BH"
)

ct_var_letters <- make_letters_from_pairwise(
  pairwise_df = ct_var_pairwise_wilcox,
  group_levels = ct_levels,
  p_col = "p_adj",
  alpha = 0.05
)

ct_var_letter_pos <- make_letter_positions(
  ct_variance_df,
  group_col = "group",
  y_col = "dist_to_centroid"
) %>%
  left_join(ct_var_letters, by = "group")

#dispersion boxplot
p_ct_var <- ggplot(ct_variance_df, aes(x = group, y = dist_to_centroid, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  geom_text(
    data = ct_var_letter_pos,
    aes(x = group, y = y_pos, label = letters),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_fill_manual(values = c("CT_1950" = "steelblue", "CT_ALL" = "navy")) +
  theme_classic() +
  labs(
    title = "Dispersion: CT 1950 vs CT (1950 + 1956 + 1970)",
    y = "Distance to centroid"
  )

ggsave(
  file.path(fig_dir, "CT_dispersion_1950_vs_all.png"),
  p_ct_var, width = 6, height = 5, dpi = 400
)

#print stats
cat("\nCT dispersion test:\n")
print(anova(bd_ct))
print(bd_ct_perm)

write.csv(ct_disparity_obs, file.path(tables_dir, "CT_1950_vs_CT_ALL_disparity_observed.csv"), row.names = FALSE)
write.csv(ct_disp_boot, file.path(tables_dir, "CT_1950_vs_CT_ALL_disparity_bootstrap.csv"), row.names = FALSE)
write.csv(ct_disp_pairwise, file.path(tables_dir, "CT_1950_vs_CT_ALL_disparity_pairwise_BH.csv"), row.names = FALSE)
write.csv(ct_disp_letters, file.path(tables_dir, "CT_1950_vs_CT_ALL_disparity_letters_BH.csv"), row.names = FALSE)
write.csv(ct_variance_df, file.path(tables_dir, "CT_1950_vs_CT_ALL_distances_to_centroid.csv"), row.names = FALSE)
write.csv(ct_var_pairwise_wilcox, file.path(tables_dir, "CT_1950_vs_CT_ALL_dispersion_pairwise_wilcox_BH.csv"), row.names = FALSE)
write.csv(ct_var_letters, file.path(tables_dir, "CT_1950_vs_CT_ALL_dispersion_letters_BH.csv"), row.names = FALSE)

sink(file.path(text_dir, "CT_1950_vs_CT_ALL_disparity_dispersion_stats.txt"))
cat("============================================================
")
cat("CT 1950 VS. LATER CT TIMEPOINTS
")
cat("============================================================

")
cat("Dataset: CT_1950 vs CT_ALL, where CT_ALL = CT_1956 + CT_1970
")
cat("Shape data: allometry-corrected residual GPA coordinates

")
cat("Group counts:
")
print(table(ct_meta$ct_group))
cat("
Observed disparity by group:
")
print(ct_disparity_obs)
cat("
Pairwise disparity permutation test with BH adjustment:
")
print(ct_disp_pairwise)
cat("
Compact letter display for CT disparity:
")
print(ct_disp_letters)
cat("
Dispersion test:
")
print(anova(bd_ct))
print(bd_ct_perm)
cat("
Pairwise Wilcoxon tests on distance to centroid with BH adjustment:
")
print(ct_var_pairwise_wilcox)
cat("
Compact letter display for CT dispersion:
")
print(ct_var_letters)
sink()

# ============================================================
# 4B) CT TIME SERIES AS REFERENCE LANDSCAPE VS 1950 HABITATS
# ============================================================

cat("\n============================================================\n")
cat("Running CT time series vs. 1950 habitat disparity/dispersion...\n")
cat("============================================================\n")

# ------------------------------------------------------------
# Build comparison groups:
#   - CT_timeseries = CT_1950 + CT_1956 + CT_1970
#   - 1950 habitats remain separate
# ------------------------------------------------------------

ct_landscape_meta <- meta_use %>%
  mutate(
    landscape_group = case_when(
      group %in% c("CT_1950", "CT_1956", "CT_1970") ~ "CT_timeseries",
      group == "Quabbin_1950" ~ "Quabbin_1950",
      group == "Swift River_1950" ~ "Swift River_1950",
      group == "Fort River_1950" ~ "Fort River_1950",
      group == "Sawmill River_1950" ~ "Sawmill River_1950",
      TRUE ~ as.character(group)
    ),
    landscape_group = factor(
      landscape_group,
      levels = c(
        "CT_timeseries",
        "Quabbin_1950",
        "Swift River_1950",
        "Fort River_1950",
        "Sawmill River_1950"
      )
    )
  ) %>%
  filter(!is.na(landscape_group))

ct_landscape_coords <- subset_coords_to_gdf(coords_use, ct_landscape_meta)

stopifnot(identical(dimnames(ct_landscape_coords)[[3]], ct_landscape_meta$specimen))

coords2d_landscape <- two.d.array(ct_landscape_coords)
rownames(coords2d_landscape) <- ct_landscape_meta$specimen

landscape_levels <- levels(ct_landscape_meta$landscape_group)

landscape_palette <- c(
  "CT_timeseries"      = "steelblue4",
  "Quabbin_1950"       = "darkgoldenrod2",
  "Swift River_1950"   = "tomato",
  "Fort River_1950"    = "black",
  "Sawmill River_1950" = "orchid3"
)

landscape_palette <- landscape_palette[landscape_levels]

# ------------------------------------------------------------
# A) Disparity: Procrustes variance within each landscape group
# ------------------------------------------------------------

fit_disp_landscape <- geomorph::procD.lm(
  ct_landscape_coords ~ landscape_group,
  data = ct_landscape_meta,
  iter = 999,
  RRPP = TRUE
)

disp_obj_landscape <- geomorph::morphol.disparity(
  fit_disp_landscape,
  groups = ct_landscape_meta$landscape_group,
  iter = 999
)

landscape_disparity_obs <- calc_group_disparity(
  coords2d_landscape,
  ct_landscape_meta$landscape_group
) %>%
  mutate(group = factor(group, levels = landscape_levels))

landscape_disp_boot <- map_dfr(landscape_levels, function(g) {
  idx <- which(ct_landscape_meta$landscape_group == g)
  Xg  <- coords2d_landscape[idx, , drop = FALSE]
  
  tibble(
    group = g,
    disparity = bootstrap_group_disparity(Xg, nboot = 999)
  )
}) %>%
  mutate(group = factor(group, levels = landscape_levels))

landscape_disp_summary <- landscape_disp_boot %>%
  group_by(group) %>%
  summarize(
    mean_boot = mean(disparity),
    sd_boot   = sd(disparity),
    q025      = quantile(disparity, 0.025),
    q25       = quantile(disparity, 0.25),
    median    = quantile(disparity, 0.50),
    q75       = quantile(disparity, 0.75),
    q975      = quantile(disparity, 0.975),
    .groups = "drop"
  ) %>%
  left_join(landscape_disparity_obs, by = "group") %>%
  select(group, n, disparity_observed = disparity, everything(), -disparity)

landscape_disp_pairwise <- make_pairwise_perm_disparity(
  coords2d = coords2d_landscape,
  groups = ct_landscape_meta$landscape_group,
  nperm = 999
)

landscape_disp_letters <- make_letters_from_pairwise(
  pairwise_df = landscape_disp_pairwise,
  group_levels = landscape_levels,
  p_col = "p_adj",
  alpha = 0.05
)

landscape_disp_letter_pos <- make_letter_positions(
  landscape_disp_boot,
  group_col = "group",
  y_col = "disparity"
) %>%
  left_join(landscape_disp_letters, by = "group")

p_landscape_disp <- ggplot(
  landscape_disp_boot,
  aes(x = group, y = disparity, fill = group)
) +
  geom_boxplot(width = 0.72, outlier.shape = NA, alpha = 0.9) +
  geom_point(
    data = landscape_disparity_obs,
    aes(x = group, y = disparity),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.8,
    stroke = 0.8,
    fill = "white",
    color = "black"
  ) +
  geom_text(
    data = landscape_disp_letter_pos,
    aes(x = group, y = y_pos, label = letters),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_fill_manual(values = landscape_palette, drop = FALSE) +
  labs(
    title = "Disparity: CT time series vs. 1950 habitats",
    subtitle = "CT_timeseries = CT 1950 + 1956 + 1970; white points = observed Procrustes variance; letters = BH-adjusted pairwise tests",
    x = NULL,
    y = "Procrustes variance"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave(
  file.path(fig_dir, "CT_timeseries_vs_1950habitats_disparity.png"),
  p_landscape_disp,
  width = 8.5,
  height = 5.5,
  dpi = 400
)

ggsave(
  file.path(fig_dir, "CT_timeseries_vs_1950habitats_disparity.pdf"),
  p_landscape_disp,
  width = 8.5,
  height = 5.5
)

# ------------------------------------------------------------
# B) Dispersion: distance to group centroid
# ------------------------------------------------------------

D_landscape <- dist(coords2d_landscape)

bd_landscape <- vegan::betadisper(
  D_landscape,
  group = ct_landscape_meta$landscape_group,
  type = "centroid",
  bias.adjust = FALSE
)

bd_landscape_perm <- vegan::permutest(
  bd_landscape,
  permutations = 999
)

bd_landscape_pair <- vegan::permutest(
  bd_landscape,
  pairwise = TRUE,
  permutations = 999
)

landscape_variance_df <- tibble(
  specimen = names(bd_landscape$distances),
  group = factor(ct_landscape_meta$landscape_group[
    match(names(bd_landscape$distances), ct_landscape_meta$specimen)
  ], levels = landscape_levels),
  dist_to_centroid = as.numeric(bd_landscape$distances)
)

landscape_variance_summary <- landscape_variance_df %>%
  group_by(group) %>%
  summarize(
    n = n(),
    mean_dist = mean(dist_to_centroid),
    sd_dist = sd(dist_to_centroid),
    median_dist = median(dist_to_centroid),
    min_dist = min(dist_to_centroid),
    max_dist = max(dist_to_centroid),
    .groups = "drop"
  )

landscape_var_pairwise_wilcox <- make_pairwise_wilcox_table(
  values = landscape_variance_df$dist_to_centroid,
  groups = landscape_variance_df$group,
  p_adjust = "BH"
)

landscape_var_letters <- make_letters_from_pairwise(
  pairwise_df = landscape_var_pairwise_wilcox,
  group_levels = landscape_levels,
  p_col = "p_adj",
  alpha = 0.05
)

landscape_var_letter_pos <- make_letter_positions(
  landscape_variance_df,
  group_col = "group",
  y_col = "dist_to_centroid"
) %>%
  left_join(landscape_var_letters, by = "group")

p_landscape_var <- ggplot(
  landscape_variance_df,
  aes(x = group, y = dist_to_centroid, fill = group)
) +
  geom_boxplot(width = 0.72, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.75) +
  geom_text(
    data = landscape_var_letter_pos,
    aes(x = group, y = y_pos, label = letters),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_fill_manual(values = landscape_palette, drop = FALSE) +
  labs(
    title = "Dispersion: CT time series vs. 1950 habitats",
    subtitle = "Distances of individuals to their group centroid",
    x = NULL,
    y = "Distance to group centroid"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave(
  file.path(fig_dir, "CT_timeseries_vs_1950habitats_dispersion.png"),
  p_landscape_var,
  width = 8.5,
  height = 5.5,
  dpi = 400
)

ggsave(
  file.path(fig_dir, "CT_timeseries_vs_1950habitats_dispersion.pdf"),
  p_landscape_var,
  width = 8.5,
  height = 5.5
)

# ------------------------------------------------------------
# C) Save outputs
# ------------------------------------------------------------

write.csv(
  landscape_disparity_obs,
  file.path(tables_dir, "CT_timeseries_vs_1950habitats_disparity_observed.csv"),
  row.names = FALSE
)

write.csv(
  landscape_disp_boot,
  file.path(tables_dir, "CT_timeseries_vs_1950habitats_disparity_bootstrap.csv"),
  row.names = FALSE
)

write.csv(
  landscape_disp_summary,
  file.path(tables_dir, "CT_timeseries_vs_1950habitats_disparity_summary.csv"),
  row.names = FALSE
)

write.csv(
  landscape_disp_pairwise,
  file.path(tables_dir, "CT_timeseries_vs_1950habitats_disparity_pairwise_BH.csv"),
  row.names = FALSE
)

write.csv(
  landscape_disp_letters,
  file.path(tables_dir, "CT_timeseries_vs_1950habitats_disparity_letters_BH.csv"),
  row.names = FALSE
)

write.csv(
  landscape_variance_df,
  file.path(tables_dir, "CT_timeseries_vs_1950habitats_distances_to_centroid.csv"),
  row.names = FALSE
)

write.csv(
  landscape_variance_summary,
  file.path(tables_dir, "CT_timeseries_vs_1950habitats_dispersion_summary.csv"),
  row.names = FALSE
)

write.csv(
  landscape_var_pairwise_wilcox,
  file.path(tables_dir, "CT_timeseries_vs_1950habitats_dispersion_pairwise_wilcox_BH.csv"),
  row.names = FALSE
)

write.csv(
  landscape_var_letters,
  file.path(tables_dir, "CT_timeseries_vs_1950habitats_dispersion_letters_BH.csv"),
  row.names = FALSE
)

saveRDS(
  fit_disp_landscape,
  file.path(rds_dir, "fit_disp_CT_timeseries_vs_1950habitats.rds")
)

saveRDS(
  disp_obj_landscape,
  file.path(rds_dir, "disp_obj_CT_timeseries_vs_1950habitats.rds")
)

saveRDS(
  bd_landscape,
  file.path(rds_dir, "betadisper_CT_timeseries_vs_1950habitats.rds")
)

saveRDS(
  bd_landscape_perm,
  file.path(rds_dir, "betadisper_permutest_CT_timeseries_vs_1950habitats.rds")
)

saveRDS(
  bd_landscape_pair,
  file.path(rds_dir, "betadisper_pairwise_CT_timeseries_vs_1950habitats.rds")
)

sink(file.path(text_dir, "CT_timeseries_vs_1950habitats_disparity_dispersion.txt"))
cat("============================================================\n")
cat("CT TIME SERIES VS. 1950 HABITATS\n")
cat("============================================================\n\n")
cat("Dataset: CT 1950 + 1956 + 1970 pooled as CT_timeseries; 1950 habitats separate\n")
cat("Shape data: allometry-corrected residual GPA coordinates\n\n")

cat("Group counts:\n")
print(table(ct_landscape_meta$landscape_group))

cat("\nObserved disparity by group:\n")
print(landscape_disparity_obs)

cat("\nBootstrap disparity summary:\n")
print(landscape_disp_summary)

cat("\nprocD.lm summary for disparity model:\n")
print(summary(fit_disp_landscape))

cat("\nmorphol.disparity object:\n")
print(disp_obj_landscape)

cat("\nDispersion summary: distance to centroid\n")
print(landscape_variance_summary)

cat("\nANOVA on multivariate dispersion:\n")
print(anova(bd_landscape))

cat("\nPermutation test on multivariate dispersion:\n")
print(bd_landscape_perm)

cat("\nPairwise permutation tests:\n")
print(bd_landscape_pair)
sink()

# ============================================================
# 5) CT ONLY: MEAN SHAPE THROUGH TIME 
# ============================================================

cat("\n============================================================\n")
cat("Running CT-only mean shape through time and trajectory analysis...\n")
cat("============================================================\n")

# make sure output directories exist before writing files
dir.create(text_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir,  recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Build CT-only dataset with actual time labels
# ------------------------------------------------------------
ct_time_meta <- meta_use %>%
  filter(group %in% c("CT_1950", "CT_1956", "CT_1970")) %>%
  mutate(
    year_f   = factor(group, levels = c("CT_1950", "CT_1956", "CT_1970")),
    year_num = c(1950, 1956, 1970)[match(group, c("CT_1950", "CT_1956", "CT_1970"))]
  )

ct_time_coords <- subset_coords_to_gdf(coords_use, ct_time_meta)
stopifnot(identical(dimnames(ct_time_coords)[[3]], ct_time_meta$specimen))

# safer 2D version for RRPP functions
Y_ct <- two.d.array(ct_time_coords)
rownames(Y_ct) <- ct_time_meta$specimen

# ------------------------------------------------------------
# A) Mean shape through time
# ------------------------------------------------------------
fit_mean_ct_time <- geomorph::procD.lm(
  ct_time_coords ~ year_f,
  data = ct_time_meta,
  iter = 999,
  RRPP = TRUE
)

pw_mean_ct_time <- RRPP::pairwise(
  fit_mean_ct_time,
  groups = ct_time_meta$year_f
)

dir.create(text_dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(text_dir, "mean_shape_1950vsalltime.txt"))
cat("============================================================\n")
cat("MEAN SHAPE THROUGH TIME: CT 1950 vs all CT timepoints\n")
cat("============================================================\n\n")
cat("Dataset: Connecticut River only\n")
cat("Timepoints: 1950, 1956, 1970\n")
cat("Shape data: allometry-corrected residual GPA coordinates\n\n")

cat("Group counts:\n")
print(table(ct_time_meta$year_f))

cat("\nOverall Procrustes ANOVA:\n")
print(summary(fit_mean_ct_time))

cat("\nPairwise comparisons among CT timepoints:\n")
print(summary(pw_mean_ct_time))
sink()

saveRDS(fit_mean_ct_time, file.path(rds_dir, "fit_mean_ct_time.rds"))
saveRDS(pw_mean_ct_time,  file.path(rds_dir, "pairwise_mean_ct_time.rds"))

# ============================================================
# 5B) CT ONLY: MANUAL TRAJECTORY DISTANCES THROUGH TIME
# ============================================================

cat("\n============================================================\n")
cat("Running CT-only manual trajectory distance analysis...\n")
cat("============================================================\n")

dir.create(text_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# CT-only data already built above:
#   ct_time_meta
#   ct_time_coords
# ------------------------------------------------------------

# helper: mean shape by group
mean_shape_by_group <- function(coords_arr, groups, level) {
  idx <- which(groups == level)
  stopifnot(length(idx) > 0)
  geomorph::mshape(coords_arr[, , idx, drop = FALSE])
}

# helper: Procrustes distance between two shapes
shape_dist <- function(shape1, shape2) {
  sqrt(sum((shape1 - shape2)^2))
}

# ------------------------------------------------------------
# Compute CT mean shapes for each timepoint
# ------------------------------------------------------------
mean_1950 <- mean_shape_by_group(ct_time_coords, ct_time_meta$year_f, "CT_1950")
mean_1956 <- mean_shape_by_group(ct_time_coords, ct_time_meta$year_f, "CT_1956")
mean_1970 <- mean_shape_by_group(ct_time_coords, ct_time_meta$year_f, "CT_1970")

# distances
d_1950_1956 <- shape_dist(mean_1950, mean_1956)
d_1956_1970 <- shape_dist(mean_1956, mean_1970)
d_1950_1970 <- shape_dist(mean_1950, mean_1970)

path_length <- d_1950_1956 + d_1956_1970
path_to_net_ratio <- path_length / d_1950_1970

trajectory_distances <- tibble::tibble(
  comparison = c(
    "CT_1950_to_CT_1956",
    "CT_1956_to_CT_1970",
    "CT_1950_to_CT_1970_net",
    "CT_path_length",
    "CT_path_to_net_ratio"
  ),
  value = c(
    d_1950_1956,
    d_1956_1970,
    d_1950_1970,
    path_length,
    path_to_net_ratio
  )
)

write.csv(
  trajectory_distances,
  file.path(tables_dir, "trajectory_distances_1950vsalltime.csv"),
  row.names = FALSE
)

sink(file.path(text_dir, "trajectory_distances_1950vsalltime.txt"))
cat("============================================================\n")
cat("MANUAL TRAJECTORY DISTANCES: CT through time\n")
cat("============================================================\n\n")
cat("Dataset: Connecticut River only\n")
cat("Ordered timepoints: 1950 -> 1956 -> 1970\n")
cat("Shape data: allometry-corrected residual GPA coordinates\n\n")

cat("Group counts:\n")
print(table(ct_time_meta$year_f))

cat("\nTrajectory distances (Procrustes distances between mean shapes):\n")
print(trajectory_distances)

cat("\nInterpretation notes:\n")
cat("- CT_1950_to_CT_1956 and CT_1956_to_CT_1970 are stepwise changes.\n")
cat("- CT_1950_to_CT_1970_net is net displacement from start to end.\n")
cat("- CT_path_length is the total path traveled through morphospace.\n")
cat("- CT_path_to_net_ratio quantifies how direct the path is.\n")
cat("  * Ratio near 1 = straighter path\n")
cat("  * Larger ratio = more curved / less direct path\n")
sink()

saveRDS(
  list(
    mean_1950 = mean_1950,
    mean_1956 = mean_1956,
    mean_1970 = mean_1970,
    trajectory_distances = trajectory_distances
  ),
  file.path(rds_dir, "trajectory_distances_ct_time.rds")
)

# ------------------------------------------------------------
# PCA plot with CT trajectory + convex hulls
# ------------------------------------------------------------

# PCA on CT specimens
ct_pca <- geomorph::gm.prcomp(ct_time_coords)

# specimen scores
ct_scores <- as.data.frame(ct_pca$x[, 1:2, drop = FALSE])
colnames(ct_scores) <- c("PC1", "PC2")
ct_scores$specimen <- rownames(ct_pca$x)

ct_scores <- dplyr::left_join(
  ct_scores,
  ct_time_meta %>% dplyr::select(specimen, year_f),
  by = "specimen"
)

# compute mean PC scores directly
mean_scores <- ct_scores %>%
  dplyr::group_by(year_f) %>%
  dplyr::summarize(
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    year_f = factor(year_f, levels = c("CT_1950", "CT_1956", "CT_1970"))
  )

# convex hull helper
hull_df <- ct_scores %>%
  dplyr::group_by(year_f) %>%
  dplyr::slice(chull(PC1, PC2)) %>%
  dplyr::ungroup()

# plot
p_ct_traj <- ggplot() +
  geom_polygon(
    data = hull_df,
    aes(x = PC1, y = PC2, fill = year_f, group = year_f),
    alpha = 0.18,
    color = NA
  ) +
  geom_point(
    data = ct_scores,
    aes(x = PC1, y = PC2, color = year_f),
    alpha = 0.45,
    size = 2
  ) +
  geom_path(
    data = mean_scores,
    aes(x = PC1, y = PC2, group = 1),
    color = "black",
    linewidth = 0.9
  ) +
  geom_point(
    data = mean_scores,
    aes(x = PC1, y = PC2, fill = year_f),
    shape = 21,
    color = "black",
    size = 4,
    stroke = 0.8
  ) +
  geom_text(
    data = mean_scores,
    aes(x = PC1, y = PC2, label = year_f),
    nudge_y = 0.01,
    size = 3.5
  ) +
  scale_color_manual(values = c(
    "CT_1950" = "steelblue",
    "CT_1956" = "deepskyblue3",
    "CT_1970" = "navy"
  )) +
  scale_fill_manual(values = c(
    "CT_1950" = "steelblue",
    "CT_1956" = "deepskyblue3",
    "CT_1970" = "navy"
  )) +
  theme_classic() +
  labs(
    title = "CT trajectory through time (PCA)",
    subtitle = "Points = individuals; hulls = timepoint distributions; black path = mean trajectory",
    x = "PC1",
    y = "PC2"
  )

ggsave(
  file.path(fig_dir, "CT_trajectory_means_PC12.png"),
  p_ct_traj,
  width = 6.5, height = 5.5, dpi = 400
)

ggsave(
  file.path(fig_dir, "CT_trajectory_means_PC12.pdf"),
  p_ct_traj,
  width = 6.5, height = 5.5
)



# ------------------------------------------------------------
# PCA plot with CT trajectory + convex hulls (PC1 inverted)
# ------------------------------------------------------------

# PCA on CT specimens
ct_pca <- geomorph::gm.prcomp(ct_time_coords)

# specimen scores
ct_scores <- as.data.frame(ct_pca$x[, 1:2, drop = FALSE])
colnames(ct_scores) <- c("PC1", "PC2")
ct_scores$specimen <- rownames(ct_pca$x)

# >>> INVERT PC1 HERE <<<
ct_scores$PC1 <- -ct_scores$PC1

# attach metadata
ct_scores <- dplyr::left_join(
  ct_scores,
  ct_time_meta %>% dplyr::select(specimen, year_f),
  by = "specimen"
)

# compute mean PC scores
mean_scores <- ct_scores %>%
  dplyr::group_by(year_f) %>%
  dplyr::summarize(
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    year_f = factor(year_f, levels = c("CT_1950", "CT_1956", "CT_1970"))
  )

# convex hulls
hull_df <- ct_scores %>%
  dplyr::group_by(year_f) %>%
  dplyr::slice(chull(PC1, PC2)) %>%
  dplyr::ungroup()

# plot
p_ct_traj <- ggplot() +
  geom_polygon(
    data = hull_df,
    aes(x = PC1, y = PC2, fill = year_f, group = year_f),
    alpha = 0.18,
    color = NA
  ) +
  geom_point(
    data = ct_scores,
    aes(x = PC1, y = PC2, color = year_f),
    alpha = 0.45,
    size = 2
  ) +
  geom_path(
    data = mean_scores,
    aes(x = PC1, y = PC2, group = 1),
    color = "black",
    linewidth = 0.9
  ) +
  geom_point(
    data = mean_scores,
    aes(x = PC1, y = PC2, fill = year_f),
    shape = 21,
    color = "black",
    size = 4,
    stroke = 0.8
  ) +
  geom_text(
    data = mean_scores,
    aes(x = PC1, y = PC2, label = year_f),
    nudge_y = 0.01,
    size = 3.5
  ) +
  scale_color_manual(values = c(
    "CT_1950" = "steelblue",
    "CT_1956" = "deepskyblue3",
    "CT_1970" = "navy"
  )) +
  scale_fill_manual(values = c(
    "CT_1950" = "steelblue",
    "CT_1956" = "deepskyblue3",
    "CT_1970" = "navy"
  )) +
  theme_classic() +
  labs(
    title = "CT trajectory through time (PCA)",
    subtitle = "Points = individuals; hulls = timepoint distributions; black path = mean trajectory",
    x = "PC1 (inverted)",
    y = "PC2"
  )

ggsave(
  file.path(fig_dir, "CT_trajectory_means_PC12_invert.png"),
  p_ct_traj,
  width = 6.5, height = 5.5, dpi = 400
)

ggsave(
  file.path(fig_dir, "CT_trajectory_means_PC12_invert.pdf"),
  p_ct_traj,
  width = 6.5, height = 5.5
)





# ---------------------------
# Manifest / run summary
# ---------------------------
manifest <- tibble::tribble(
  ~category, ~file,
  "figure", file.path(fig_dir, "disparity_boxplot_Fig6groups.png"),
  "figure", file.path(fig_dir, "disparity_boxplot_Fig6groups.pdf"),
  "figure", file.path(fig_dir, "variance_boxplot_Fig6groups.png"),
  "figure", file.path(fig_dir, "variance_boxplot_Fig6groups.pdf"),
  "table",  file.path(tables_dir, "disparity_observed_by_group.csv"),
  "table",  file.path(tables_dir, "disparity_bootstrap_distributions.csv"),
  "table",  file.path(tables_dir, "disparity_summary.csv"),
  "table",  file.path(tables_dir, "variance_distances_to_centroid.csv"),
  "table",  file.path(tables_dir, "variance_summary.csv"),
  "text",   file.path(text_dir, "disparity_analysis.txt"),
  "text",   file.path(text_dir, "variance_analysis.txt"),
  "rds",    file.path(rds_dir, "fit_disp_procDlm.rds"),
  "rds",    file.path(rds_dir, "disp_obj_morphol_disparity.rds"),
  "rds",    file.path(rds_dir, "betadisper_obj.rds"),
  "rds",    file.path(rds_dir, "betadisper_permutest.rds"),
  "rds",    file.path(rds_dir, "betadisper_pairwise.rds"),
  "text",   file.path(text_dir, "mean_shape_1950vsalltime.txt"),
  "text",   file.path(text_dir, "trajectory_1950vsalltime.txt"),
  "rds",    file.path(rds_dir, "fit_mean_ct_time.rds"),
  "rds",    file.path(rds_dir, "pairwise_mean_ct_time.rds"),
  "table",  file.path(tables_dir, "trajectory_distances_1950vsalltime.csv"),
  "text",   file.path(text_dir, "trajectory_distances_1950vsalltime.txt"),
  "figure", file.path(fig_dir, "CT_trajectory_means_PC12.png"),
  "figure", file.path(fig_dir, "CT_trajectory_means_PC12.pdf"),
  "rds",    file.path(rds_dir, "trajectory_distances_ct_time.rds")
)

if (!is.null(pairwise_df)) {
  manifest <- bind_rows(
    manifest,
    tibble(category = "table", file = file.path(tables_dir, "variance_pairwise_permutation_tests.csv"))
  )
}

write.csv(manifest, file.path(out_root, "MANIFEST.csv"), row.names = FALSE)

sink(file.path(out_root, "RUN_INFO.txt"))
cat("Run ID:", run_id, "\n")
cat("Working directory:", getwd(), "\n")
cat("Date:", as.character(Sys.time()), "\n\n")
cat("Analysis dataset:\n")
print(table(meta_use$group))
cat("\nSession info:\n")
print(sessionInfo())
sink()

cat("\n============================================================\n")
cat("Done.\n")
cat("Outputs written to:\n")
cat("  ", out_root, "\n")
cat("Figures written to:\n")
cat("  ", fig_dir, "\n")
cat("============================================================\n")