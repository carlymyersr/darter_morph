# ============================================================
# Scripts/_run_modularity_integration_CT3.R
# CT3 Modularity + Integration test runner (RAW + RESID)
#
# Upstream canonical scripts:
#   - R/00_setup_morpho.R
#   - R/01_build_metadata.R
#   - R/03_subset_CT_timeseries.R
#
# Expected objects after sourcing:
#   coords_CT3, coords_resid_CT3, gdf_CT3
#
# Output structure:
#   Outputs/output_modularity_raw_CT3/<run_id>/<test_id>/...
#   Outputs/output_modularity_resid_CT3/<run_id>/<test_id>/...
# ============================================================

suppressPackageStartupMessages({
  library(geomorph)
})

# ---------------------------
# 0) Load canonical CT3 objects if missing
# ---------------------------
ensure_canonical_CT3 <- function() {
  need <- c("coords_CT3", "coords_resid_CT3", "gdf_CT3")
  missing <- need[!vapply(need, exists, logical(1), envir = .GlobalEnv)]
  if (length(missing) > 0) {
    message("Sourcing canonical scripts (missing: ", paste(missing, collapse = ", "), ")")
    source("R/00_setup_morpho.R")
    source("R/01_build_metadata.R")
    source("R/03_subset_CT_timeseries.R")
  }
  # sanity: alignment
  stopifnot(identical(dimnames(coords_CT3)[[3]], gdf_CT3$specimen))
  stopifnot(identical(dimnames(coords_resid_CT3), dimnames(coords_CT3)))
  invisible(TRUE)
}

# ---------------------------
# 1) IO helpers
# ---------------------------
write_lines_safe <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, con = path)
  message("Wrote: ", normalizePath(path, winslash = "/"))
}

write_csv_safe <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, path, row.names = FALSE)
  message("Wrote: ", normalizePath(path, winslash = "/"))
}

safe_id <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  sub("^_|_$", "", x)
}

# ---------------------------
# 2) Partition helpers
# ---------------------------
assert_points_exist <- function(module_list, pt_names, test_id, coords_label) {
  all_pts <- unique(unlist(module_list, use.names = FALSE))
  missing <- setdiff(all_pts, pt_names)
  if (length(missing) > 0) {
    stop(sprintf("[%s | %s] Missing %d landmark(s): %s",
                 test_id, coords_label, length(missing), paste(missing, collapse = ", ")))
  }
  invisible(TRUE)
}

assert_semilandmark_curves_intact <- function(module_list, test_id) {
  required_curves <- list(
    cranium_orbital = c(
      "cranium_orbital_start",
      paste0("cranium_orbital_sl", 1:8),
      "cranium_orbital_end"
    ),
    hyoid_pelvic = c(
      "hyoid_pelvic_start",
      paste0("hyoid_pelvic_sl", 1:8),
      "hyoid_pelvic_end"
    )
  )
  point_module <- unlist(lapply(names(module_list), function(module_name) {
    setNames(rep(module_name, length(module_list[[module_name]])), module_list[[module_name]])
  }))
  for (curve_name in names(required_curves)) {
    curve_points <- required_curves[[curve_name]]
    modules <- unique(unname(point_module[curve_points]))
    modules <- modules[!is.na(modules)]
    if (length(modules) != 1) {
      stop(sprintf(
        "[%s] Semilandmark curve '%s' is split across modules: %s",
        test_id, curve_name, paste(modules, collapse = ", ")
      ))
    }
  }
  invisible(TRUE)
}

partition_from_modules <- function(module_list, pt_names) {
  part <- rep(NA_character_, length(pt_names))
  names(part) <- pt_names
  for (mn in names(module_list)) part[module_list[[mn]]] <- mn
  if (anyNA(part)) {
    stop("partition_from_modules(): some points were not assigned to any module (this should not happen after subsetting).")
  }
  factor(part, levels = names(module_list))
}

# ---------------------------
# 3) Version-aware integration + PLS
# ---------------------------
run_integration_partitioned <- function(A, part, iter = 999, seed = 1, print.progress = FALSE) {
  if (!exists("integration.test", where = asNamespace("geomorph"), inherits = FALSE)) {
    warning("geomorph::integration.test() not found; using integration.Vrel(A) (unpartitioned).")
    return(geomorph::integration.Vrel(A))
  }
  geomorph::integration.test(
    A = A,
    partition.gp = part,
    iter = iter,
    seed = seed,
    print.progress = print.progress
  )
}

run_two_b_pls_pair <- function(A1, A2, iter = 999, seed = 1, print.progress = FALSE) {
  f <- geomorph::two.b.pls
  fml <- names(formals(f))
  
  # Your 1950 framework uses A1/A2 interface — keep that
  if (all(c("A1", "A2") %in% fml)) {
    return(f(A1 = A1, A2 = A2, iter = iter, seed = seed, print.progress = print.progress))
  }
  
  stop("geomorph::two.b.pls() signature is not A1/A2 in this environment; tell me the error and I’ll patch.")
}

# ---------------------------
# 4) Init run (CT3 outputs)
# ---------------------------
init_modularity_run_CT3 <- function() {
  run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out_root_raw   <- file.path("Outputs", "output_modularity_raw_CT3",   run_id)
  out_root_resid <- file.path("Outputs", "output_modularity_resid_CT3", run_id)
  dir.create(out_root_raw,   recursive = TRUE, showWarnings = FALSE)
  dir.create(out_root_resid, recursive = TRUE, showWarnings = FALSE)
  list(run_id = run_id, out_root_raw = out_root_raw, out_root_resid = out_root_resid)
}

# ---------------------------
# 5) Core: run one test on one coords array
# ---------------------------
run_one_test <- function(coords_arr, test_id, module_list, out_root, coords_label,
                         iter = 999, CI = TRUE, seed = 1, print_progress = FALSE) {
  
  set.seed(seed)
  
  test_id <- safe_id(test_id)
  out_dir <- file.path(out_root, test_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # validate points exist
  pt_names <- dimnames(coords_arr)[[1]]
  assert_points_exist(module_list, pt_names, test_id, coords_label)
  assert_semilandmark_curves_intact(module_list, test_id)
  
  # subset coords to points in this test
  pts_in_test <- unique(unlist(module_list, use.names = FALSE))
  coords_arr  <- coords_arr[pts_in_test, , , drop = FALSE]
  
  # build partition (must cover all points now)
  part <- partition_from_modules(module_list, pt_names = pts_in_test)
  
  # MODEL CARD
  card_path <- file.path(out_dir, paste0(test_id, "_MODEL_CARD.txt"))
  lines <- c(
    paste0("run_id: ", basename(out_root)),
    paste0("test_id: ", test_id),
    paste0("coords_type: ", coords_label),
    "",
    paste0("n_specimens: ", dim(coords_arr)[3]),
    paste0("p_landmarks_used: ", dim(coords_arr)[1]),
    "",
    "module_definitions:"
  )
  for (mn in names(module_list)) {
    lines <- c(lines, paste0("  - ", mn, " (n=", length(module_list[[mn]]), ")"))
  }
  lines <- c(lines, "", "sessionInfo():", capture.output(sessionInfo()))
  write_lines_safe(lines, card_path)
  
  # ---- Modularity (CR) ----
  mod_obj <- geomorph::modularity.test(
    A = coords_arr,
    partition.gp = part,
    iter = iter,
    CI = CI
  )
  mod_txt <- file.path(out_dir, paste0(test_id, "_modularity_test_CR.txt"))
  write_lines_safe(capture.output(print(mod_obj)), mod_txt)
  saveRDS(mod_obj, file.path(out_dir, paste0(test_id, "_modularity_test_CR.rds")))
  
  # ---- Integration ----
  int_obj <- run_integration_partitioned(
    coords_arr, part,
    iter = iter, seed = seed,
    print.progress = print_progress
  )
  int_txt <- file.path(out_dir, paste0(test_id, "_integration_test.txt"))
  write_lines_safe(capture.output(print(int_obj)), int_txt)
  saveRDS(int_obj, file.path(out_dir, paste0(test_id, "_integration_test.rds")))
  
  # ---- Pairwise two.b.pls across module pairs ----
  mods <- names(module_list)
  if (length(mods) >= 2) {
    pairs <- utils::combn(mods, 2)
    sum_rows <- vector("list", ncol(pairs))
    pls_objs <- vector("list", ncol(pairs))
    
    for (i in seq_len(ncol(pairs))) {
      m1 <- pairs[1, i]
      m2 <- pairs[2, i]
      A1 <- coords_arr[module_list[[m1]], , , drop = FALSE]
      A2 <- coords_arr[module_list[[m2]], , , drop = FALSE]
      
      pls <- run_two_b_pls_pair(A1, A2, iter = iter, seed = seed, print.progress = print_progress)
      pls_objs[[i]] <- pls
      
      r_val <- tryCatch(as.numeric(pls$r.pls), error = function(e) NA_real_)
      z_val <- tryCatch(as.numeric(pls$Z), error = function(e) NA_real_)
      p_val <- tryCatch(as.numeric(pls$P.value), error = function(e) NA_real_)
      
      sum_rows[[i]] <- data.frame(
        module1 = m1, module2 = m2,
        r_pls = r_val, Z = z_val, P_value = p_val,
        stringsAsFactors = FALSE
      )
    }
    
    pls_txt <- file.path(out_dir, paste0(test_id, "_pairwise_two_b_pls.txt"))
    write_lines_safe(
      unlist(lapply(seq_along(pls_objs), function(i) {
        hdr <- paste0("\n===== two.b.pls: ",
                      sum_rows[[i]]$module1[1], " vs ", sum_rows[[i]]$module2[1],
                      " =====")
        c(hdr, capture.output(print(pls_objs[[i]])))
      })),
      pls_txt
    )
    saveRDS(pls_objs, file.path(out_dir, paste0(test_id, "_pairwise_two_b_pls.rds")))
    
    pls_sum <- do.call(rbind, sum_rows)
    write_csv_safe(pls_sum, file.path(out_dir, paste0(test_id, "_pairwise_two_b_pls_summary.csv")))
  }
  
  invisible(out_dir)
}

# ---------------------------
# 6) Public entrypoint: run RAW + RESID for CT3
# ---------------------------
run_test_raw_resid_CT3 <- function(test_id, module_list, iter = 999, CI = TRUE, seed = 1, print_progress = FALSE) {
  ensure_canonical_CT3()
  run <- init_modularity_run_CT3()
  
  run_one_test(coords_CT3, test_id, module_list, run$out_root_raw, "RAW",
               iter = iter, CI = CI, seed = seed, print_progress = print_progress)
  
  run_one_test(coords_resid_CT3, test_id, module_list, run$out_root_resid, "RESID",
               iter = iter, CI = CI, seed = seed, print_progress = print_progress)
  
  message("DONE.")
  message("RAW outputs:   ", normalizePath(run$out_root_raw, winslash = "/"))
  message("RESID outputs: ", normalizePath(run$out_root_resid, winslash = "/"))
  invisible(run)
}
