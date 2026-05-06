# ============================================================
# Scripts/1950_stats_models.R
# 1950 STATISTICAL MODELS (RAW + RESIDUAL SHAPES)
#
# Canonical objects created upstream:
#   - R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R
#   - R/methods/01_specimen_sampling_study_design/01_build_metadata.R
#   - R/methods/01_specimen_sampling_study_design/02_subset_1950.R
# ============================================================

# ---------------------------
# Load canonical project objects
# ---------------------------
source("R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
source("R/methods/01_specimen_sampling_study_design/01_build_metadata.R")
source("R/methods/01_specimen_sampling_study_design/02_subset_1950.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(geomorph)
  library(RRPP)
  library(scales)
})

# ============================================================
# 0) Output directories + run ID + manifest
# ============================================================

run_id   <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_base <- file.path("Outputs", paste0("1950_models_", run_id))
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

manifest <- data.frame(
  run_id = character(),
  model_id = character(),
  artifact_type = character(),
  file = character(),
  stringsAsFactors = FALSE
)

add_manifest <- function(model_id, artifact_type, file) {
  manifest <<- rbind(
    manifest,
    data.frame(
      run_id = run_id,
      model_id = model_id,
      artifact_type = artifact_type,
      file = file,
      stringsAsFactors = FALSE
    )
  )
}

safe_id <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  sub("^_|_$", "", x)
}

model_dir <- function(model_id) {
  d <- file.path(out_base, model_id)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

write_lines_safe <- function(lines, path) {
  writeLines(lines, con = path)
  message("Wrote: ", normalizePath(path))
}

write_csv_safe <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE)
  message("Wrote: ", normalizePath(path))
}

# ============================================================
# 0a) Preconditions / sanity
# ============================================================

stopifnot(exists("gdf_1950"), exists("coords_1950"), exists("coords_resid_1950"))
stopifnot(identical(dimnames(coords_1950)[[3]], gdf_1950$specimen))
stopifnot(identical(dimnames(coords_resid_1950), dimnames(coords_1950)))

if (!("size_for_allometry" %in% names(gdf_1950))) {
  stop("gdf_1950 must contain size_for_allometry (set in R/methods/01_specimen_sampling_study_design/01_build_metadata.R).")
}

# stable habitat levels if palette exists
if (exists("hab_palette")) {
  gdf_1950$habitat <- factor(gdf_1950$habitat, levels = names(hab_palette))
} else {
  gdf_1950$habitat <- factor(gdf_1950$habitat)
}

gdf_1950$size_for_allometry <- as.numeric(gdf_1950$size_for_allometry)

size_label_txt <- if ("size_label" %in% names(gdf_1950)) unique(gdf_1950$size_label)[1] else "size_for_allometry"

# ============================================================
# Helper: run procD.lm (geomorph version expects f1=)
# ============================================================

run_procD <- function(coords_arr, rhs, gdf_sub, iter = 999, RRPP = TRUE) {
  stopifnot(is.character(rhs), length(rhs) == 1)
  f <- stats::as.formula(paste("coords_arr ~", rhs))
  environment(f) <- environment()
  geomorph::procD.lm(f1 = f, data = gdf_sub, iter = iter, RRPP = RRPP)
}

# ============================================================
# Helper: write model card
# ============================================================

write_model_card <- function(model_id, coords_type, formula_text, gdf_sub,
                             size_covariate = NULL, size_label = NULL, notes = NULL) {
  d <- model_dir(model_id)
  path <- file.path(d, paste0(model_id, "_MODEL_CARD.txt"))
  
  lines <- c(
    paste0("run_id: ", run_id),
    paste0("model_id: ", model_id),
    paste0("coords_type: ", coords_type),
    paste0("formula: ", formula_text),
    "",
    paste0("n: ", nrow(gdf_sub)),
    "habitat_counts:",
    capture.output(print(table(gdf_sub$habitat, useNA = "ifany")))
  )
  
  if (!is.null(size_covariate)) lines <- c(lines, "", paste0("size_covariate: ", size_covariate))
  if (!is.null(size_label))     lines <- c(lines, paste0("size_label: ", size_label))
  if (!is.null(notes))          lines <- c(lines, "", "notes:", notes)
  
  lines <- c(lines, "", "sessionInfo():", capture.output(sessionInfo()))
  
  write_lines_safe(lines, path)
  add_manifest(model_id, "model_card_txt", path)
  invisible(path)
}

# ============================================================
# Helper: save fit + anova summary
# ============================================================

save_fit_bundle <- function(model_id, fit) {
  d <- model_dir(model_id)
  
  txt_path <- file.path(d, paste0(model_id, "_summary.txt"))
  lines <- c(
    paste0("run_id: ", run_id),
    paste0("model_id: ", model_id),
    "",
    "===== fit print =====",
    capture.output(print(fit)),
    "",
    "===== anova(fit) =====",
    capture.output(anova(fit))
  )
  write_lines_safe(lines, txt_path)
  add_manifest(model_id, "summary_txt", txt_path)
  
  rds_path <- file.path(d, paste0(model_id, "_fit.rds"))
  saveRDS(fit, rds_path)
  add_manifest(model_id, "fit_rds", rds_path)
  
  invisible(list(txt = txt_path, rds = rds_path))
}

# ============================================================
# Helper: mean shapes + simple Procrustes mean distance
# ============================================================

hab_mean_shapes <- function(coords_arr, hab_vec) {
  hab_vec <- factor(hab_vec)
  levs <- levels(hab_vec)
  out <- setNames(vector("list", length(levs)), levs)
  
  for (h in levs) {
    idx <- which(hab_vec == h)
    M <- apply(coords_arr[, , idx, drop = FALSE], c(1, 2), mean)
    rownames(M) <- dimnames(coords_arr)[[1]]
    colnames(M) <- c("x", "y")
    out[[h]] <- M
  }
  out
}

procdist_mean <- function(M1, M2) sqrt(sum((M1 - M2)^2))

write_mean_shape_distance_tables <- function(model_id, coords_arr, gdf_sub) {
  d <- model_dir(model_id)
  
  means <- hab_mean_shapes(coords_arr, gdf_sub$habitat)
  grand <- geomorph::mshape(coords_arr)
  
  dist_to_grand <- data.frame(
    habitat = names(means),
    dist_to_grand = vapply(means, function(M) procdist_mean(M, grand), numeric(1)),
    stringsAsFactors = FALSE
  )
  csv1 <- file.path(d, paste0(model_id, "_meanShapeDist_toGrand.csv"))
  write_csv_safe(dist_to_grand, csv1)
  add_manifest(model_id, "csv_meanShapeDist_toGrand", csv1)
  
  levs <- names(means)
  pairs <- utils::combn(levs, 2)
  df_pairs <- data.frame(
    hab1 = pairs[1, ],
    hab2 = pairs[2, ],
    dist = NA_real_,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(df_pairs))) {
    df_pairs$dist[i] <- procdist_mean(means[[df_pairs$hab1[i]]], means[[df_pairs$hab2[i]]])
  }
  csv2 <- file.path(d, paste0(model_id, "_meanShapeDist_pairwise.csv"))
  write_csv_safe(df_pairs, csv2)
  add_manifest(model_id, "csv_meanShapeDist_pairwise", csv2)
  
  invisible(list(to_grand = csv1, pairwise = csv2))
}

# ============================================================
# Helper: morphol.disparity -> df
# ============================================================

disp_to_df <- function(disp_obj) {
  if (is.list(disp_obj) && !is.null(disp_obj$Procrustes.var)) {
    pv <- disp_obj$Procrustes.var
    return(data.frame(
      group = names(pv),
      Procrustes.var = as.numeric(pv),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(group = NA_character_, Procrustes.var = NA_real_, stringsAsFactors = FALSE)
}

# ============================================================
# Helper: morphol.disparity (ROBUST EVAL)
#   - assigns coords array into .GlobalEnv temporarily
# ============================================================

write_disparity_outputs <- function(model_id, coords_arr, gdf_sub, iter = 999) {
  d <- model_dir(model_id)
  
  tmp_name <- paste0("._tmp_coords_", safe_id(model_id), "_", as.integer(Sys.time()))
  assign(tmp_name, coords_arr, envir = .GlobalEnv)
  on.exit({
    if (exists(tmp_name, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = tmp_name, envir = .GlobalEnv)
    }
  }, add = TRUE)
  
  f <- stats::as.formula(paste(tmp_name, "~ habitat"))
  
  disp <- geomorph::morphol.disparity(
    f1     = f,
    groups = gdf_sub$habitat,
    data   = gdf_sub,
    iter   = iter
  )
  
  txt_path <- file.path(d, paste0(model_id, "_morpholDisparity_print.txt"))
  write_lines_safe(capture.output(print(disp)), txt_path)
  add_manifest(model_id, "disparity_txt", txt_path)
  
  df <- disp_to_df(disp)
  csv_path <- file.path(d, paste0(model_id, "_morpholDisparity_byHabitat.csv"))
  write_csv_safe(df, csv_path)
  add_manifest(model_id, "disparity_csv", csv_path)
  
  rds_path <- file.path(d, paste0(model_id, "_morpholDisparity.rds"))
  saveRDS(disp, rds_path)
  add_manifest(model_id, "disparity_rds", rds_path)
  
  invisible(list(txt = txt_path, csv = csv_path, rds = rds_path))
}

# ============================================================
# Helper: RRPP pairwise habitat tests (MEAN SHAPE differences)
#   - uses sink() instead of capture.output() (ROBUST)
# ============================================================

write_pairwise_rrpp <- function(model_id, fit, groups) {
  d <- model_dir(model_id)
  
  pw <- RRPP::pairwise(fit, groups = groups)
  
  # ---- Always save the object ----
  rds_path <- file.path(d, paste0(model_id, "_pairwise_habitat_RRPP.rds"))
  saveRDS(pw, rds_path)
  add_manifest(model_id, "pairwise_rds", rds_path)
  
  # ---- 1) Write a readable TXT (whatever print/summary gives us) ----
  txt_path <- file.path(d, paste0(model_id, "_pairwise_habitat_RRPP.txt"))
  sink(txt_path)
  on.exit(sink(), add = TRUE)
  cat("RRPP::pairwise output (print):\n\n")
  print(pw)
  cat("\n\nRRPP::pairwise output (summary):\n\n")
  s <- summary(pw, stat.table = TRUE, test.type = "dist")
  print(s)
  sink()
  add_manifest(model_id, "pairwise_txt", txt_path)
  
  # ---- 2) TRY HARD to extract the actual distance/stat table ----
  # RRPP versions differ; check common locations
  s_list <- s
  tab <- NULL
  
  # Common candidates across versions:
  candidates <- list(
    s_list$summary.table,
    s_list$stat.table,
    s_list$tables,
    s_list$pairwise.tables,
    s_list$pairwise.table,
    s_list$pairwise,
    pw$summary.table,
    pw$stat.table,
    pw$tables
  )
  
  for (obj in candidates) {
    if (is.null(obj)) next
    
    # Sometimes it's a data.frame directly
    if (is.data.frame(obj)) { tab <- obj; break }
    
    # Sometimes it's a list of tables; take the first data.frame found
    if (is.list(obj)) {
      for (x in obj) {
        if (is.data.frame(x)) { tab <- x; break }
      }
      if (!is.null(tab)) break
    }
  }
  
  # If we found a table, write it
  if (!is.null(tab)) {
    csv_path <- file.path(d, paste0(model_id, "_pairwise_habitat_RRPP_table.csv"))
    utils::write.csv(tab, csv_path, row.names = FALSE)
    message("Wrote: ", normalizePath(csv_path))
    add_manifest(model_id, "pairwise_csv", csv_path)
  } else {
    # Fallback: write structure so we can see where the table lives
    debug_path <- file.path(d, paste0(model_id, "_pairwise_habitat_RRPP_STR.txt"))
    write_lines_safe(c(
      "Could not auto-extract a data.frame stat table from summary(pw).",
      "Here is str(pw):",
      capture.output(str(pw, max.level = 3)),
      "",
      "Here is str(summary(pw)):",
      capture.output(str(s, max.level = 3)),
      "",
      "TIP: search these str() outputs for 'table', 'dist', 'p.value', 'Z', 'Pr', or similar."
    ), debug_path)
    add_manifest(model_id, "pairwise_str_txt", debug_path)
  }
  
  invisible(list(pw = pw, summary_obj = s))
}

# ============================================================
# SIZE BUNDLE (RESTORED)
# ============================================================

model_id_size <- safe_id("SIZE_summaries_1950")
write_model_card(
  model_id = model_id_size,
  coords_type = "RAW",
  formula_text = "N/A (descriptive size summaries / plots)",
  gdf_sub = gdf_1950,
  size_covariate = "logCsize",
  size_label = "logCsize / centroid size from GPA",
  notes = c("Contains descriptive size summaries and a size distribution plot; not a procD model.")
)

d_size <- model_dir(model_id_size)

size_summary <- gdf_1950 %>%
  group_by(habitat) %>%
  summarize(
    n = n(),
    mean_logCsize = mean(logCsize, na.rm = TRUE),
    sd_logCsize   = sd(logCsize, na.rm = TRUE),
    mean_Csize    = mean(Csize, na.rm = TRUE),
    sd_Csize      = sd(Csize, na.rm = TRUE),
    mean_SL_mm    = mean(SL_mm, na.rm = TRUE),
    sd_SL_mm      = sd(SL_mm, na.rm = TRUE),
    .groups = "drop"
  )

csv_size <- file.path(d_size, paste0(model_id_size, "_size_summary_by_habitat.csv"))
write_csv_safe(size_summary, csv_size)
add_manifest(model_id_size, "csv_size_summary", csv_size)

if (exists("hab_palette")) {
  p_size <- ggplot(gdf_1950, aes(habitat, logCsize, fill = habitat)) +
    geom_boxplot(alpha = 0.25, outlier.alpha = 0.6) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.75, size = 1.6) +
    scale_fill_manual(values = hab_palette, drop = FALSE) +
    labs(
      title = "1950 centroid size distribution by habitat",
      x = "Habitat",
      y = "log(Centroid size)"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
} else {
  p_size <- ggplot(gdf_1950, aes(habitat, logCsize)) +
    geom_boxplot(alpha = 0.25, outlier.alpha = 0.6) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.75, size = 1.6) +
    labs(
      title = "1950 centroid size distribution by habitat",
      x = "Habitat",
      y = "log(Centroid size)"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
}

png_size <- file.path(d_size, paste0(model_id_size, "_logCsize_by_habitat.png"))
ggsave(png_size, p_size, width = 8.5, height = 4.75, dpi = 300)
add_manifest(model_id_size, "figure_png_size", png_size)

# ============================================================
# Identify the largest Swift River specimen
# ============================================================

swift_outlier <- gdf_1950 %>%
  filter(habitat == "Swift River") %>%
  arrange(desc(logCsize)) %>%
  select(specimen, habitat, logCsize, Csize, SL_mm) %>%
  slice(1)

print(swift_outlier)

# ============================================================
# MODELS
# ============================================================

# MODEL 1: RAW shape ~ size
model_id_m1_raw <- safe_id("RAW_Model1_shape_by_size_1950")
write_model_card(
  model_id = model_id_m1_raw,
  coords_type = "RAW",
  formula_text = "coords_1950 ~ size_for_allometry",
  gdf_sub = gdf_1950,
  size_covariate = "size_for_allometry",
  size_label = size_label_txt
)
fit_raw_m1 <- run_procD(coords_1950, "size_for_allometry", gdf_1950, iter = 999)
save_fit_bundle(model_id_m1_raw, fit_raw_m1)

# MODEL 2: RAW shape ~ habitat
model_id_m2_raw <- safe_id("RAW_Model2_shape_by_habitat_1950")
write_model_card(
  model_id = model_id_m2_raw,
  coords_type = "RAW",
  formula_text = "coords_1950 ~ habitat",
  gdf_sub = gdf_1950,
  size_covariate = "size_for_allometry",
  size_label = size_label_txt,
  notes = c("Includes RRPP pairwise mean-shape comparisons, mean-shape distance tables, and morphol.disparity.")
)
fit_raw_m2 <- run_procD(coords_1950, "habitat", gdf_1950, iter = 999)
save_fit_bundle(model_id_m2_raw, fit_raw_m2)
write_pairwise_rrpp(model_id_m2_raw, fit_raw_m2, gdf_1950$habitat)
write_mean_shape_distance_tables(model_id_m2_raw, coords_1950, gdf_1950)
write_disparity_outputs(model_id_m2_raw, coords_1950, gdf_1950)

# MODEL 3: RAW shape ~ size * habitat
model_id_m3_raw <- safe_id("RAW_Model3_shape_by_sizeXhabitat_1950")
write_model_card(
  model_id = model_id_m3_raw,
  coords_type = "RAW",
  formula_text = "coords_1950 ~ size_for_allometry * habitat",
  gdf_sub = gdf_1950,
  size_covariate = "size_for_allometry",
  size_label = size_label_txt
)
fit_raw_m3 <- run_procD(coords_1950, "size_for_allometry * habitat", gdf_1950, iter = 999)
save_fit_bundle(model_id_m3_raw, fit_raw_m3)

# MODEL 2 (RESID): residual shape ~ habitat
model_id_m2_res <- safe_id("RESID_Model2_shape_by_habitat_1950")
write_model_card(
  model_id = model_id_m2_res,
  coords_type = "RESID",
  formula_text = "coords_resid_1950 ~ habitat",
  gdf_sub = gdf_1950,
  size_covariate = "size_for_allometry",
  size_label = size_label_txt,
  notes = c("Includes RRPP pairwise mean-shape comparisons, mean-shape distance tables, and morphol.disparity (on residual shapes).")
)
fit_res_m2 <- run_procD(coords_resid_1950, "habitat", gdf_1950, iter = 999)
save_fit_bundle(model_id_m2_res, fit_res_m2)
write_pairwise_rrpp(model_id_m2_res, fit_res_m2, gdf_1950$habitat)
write_mean_shape_distance_tables(model_id_m2_res, coords_resid_1950, gdf_1950)
write_disparity_outputs(model_id_m2_res, coords_resid_1950, gdf_1950)

# ============================================================
# MANIFEST
# ============================================================

manifest_path <- file.path(out_base, "MANIFEST.csv")
write_csv_safe(manifest, manifest_path)
message("DONE. Outputs in: ", normalizePath(out_base))
