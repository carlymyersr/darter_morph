# ============================================================
# Scripts/_run_modularity_integration_1950.R
# Shared runner: Modularity + Integration + pairwise PLS (RAW + RESID)
# for the 1950 subset.
#
# This file is SAFE to source: it defines functions only (no tests run).
#
# Expected canonical upstream scripts:
#   R/00_setup_morpho.R
#   R/01_build_metadata.R
#   R/02_subset_1950.R
#
# It will auto-source them if the required objects are not present.
# ============================================================

suppressPackageStartupMessages({
  library(geomorph)
})



# ---------------------------
# 0) Canonical object loader (idempotent)
# ---------------------------
ensure_canonical_1950 <- function() {
  need <- c("coords_1950", "coords_resid_1950", "gdf_1950")
  missing <- need[!vapply(need, exists, logical(1), envir = .GlobalEnv)]
  if (length(missing) > 0) {
    message("Loading canonical upstream scripts because these objects were missing: ",
            paste(missing, collapse = ", "))
    source("R/00_setup_morpho.R")
    source("R/01_build_metadata.R")
    source("R/02_subset_1950.R")
  }
  invisible(TRUE)
}

# ---------------------------
# 1) Small IO helpers (parent-dir safe)
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
# 2) Validation + partition helpers
# ---------------------------
assert_points_exist <- function(module_list, pt_names, test_id, coords_label) {
  all_pts <- unique(unlist(module_list, use.names = FALSE))
  missing <- setdiff(all_pts, pt_names)
  if (length(missing) > 0) {
    stop(sprintf(
      "[%s | %s] Missing %d landmarks in coords dimnames: %s",
      test_id, coords_label, length(missing), paste(missing, collapse = ", ")
    ))
  }
  invisible(TRUE)
}

assert_no_dupes_across_modules <- function(module_list, test_id) {
  all_pts <- unlist(module_list, use.names = FALSE)
  dupes <- unique(all_pts[duplicated(all_pts)])
  if (length(dupes) > 0) {
    warning(sprintf(
      "[%s] The following points appear in multiple modules (allowed but usually unintended): %s",
      test_id, paste(dupes, collapse = ", ")
    ))
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

partition_from_modules <- function(module_list) {
  # Named factor: names = landmarks; values = module labels
  pt_names <- unique(unlist(module_list, use.names = FALSE))
  part <- rep(NA_character_, length(pt_names))
  names(part) <- pt_names
  for (mn in names(module_list)) {
    part[module_list[[mn]]] <- mn
  }
  factor(part, levels = names(module_list))
}

# ---------------------------
# 3) Integration (version: uses integration.test if present)
# ---------------------------
run_integration_partitioned <- function(A, part, iter = 999, seed = 1, print.progress = FALSE) {
  if (!exists("integration.test", where = asNamespace("geomorph"), inherits = FALSE)) {
    warning("geomorph::integration.test() not found. Falling back to integration.Vrel(A) (unpartitioned).")
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

# ---------------------------
# 4) Pairwise PLS across modules (geomorph version-aware)
#    - Your current geomorph has two.b.pls(A1, A2, iter, seed, print.progress)
#    - This wrapper also supports newer signatures if you ever update geomorph.
# ---------------------------
run_two_b_pls_pair <- function(A1, A2, iter = 999, seed = 1, print.progress = FALSE) {
  f <- geomorph::two.b.pls
  fml <- names(formals(f))

  # Older geomorph: expects A1/A2 blocks
  if (all(c("A1", "A2") %in% fml)) {
    return(f(A1 = A1, A2 = A2, iter = iter, seed = seed, print.progress = print.progress))
  }

  # Newer geomorph: can accept a combined array + partition.gp
  if (all(c("A", "partition.gp") %in% fml)) {
    # Combine blocks and build a 2-module partition
    pts1 <- dimnames(A1)[[1]]
    pts2 <- dimnames(A2)[[1]]
    A <- abind::abind(A1, A2, along = 1)  # requires abind if you ever hit this branch
    part <- c(rep("M1", length(pts1)), rep("M2", length(pts2)))
    names(part) <- c(pts1, pts2)
    part <- factor(part, levels = c("M1","M2"))
    args <- list(A = A, partition.gp = part, iter = iter)
    if ("seed" %in% fml) args$seed <- seed
    if ("print.progress" %in% fml) args$print.progress <- print.progress
    return(do.call(f, args))
  }

  stop("Unsupported geomorph::two.b.pls signature in this environment.")
}

# ---------------------------
# 5) Run initialization
# ---------------------------
init_modularity_run <- function() {
  run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out_root_raw   <- file.path("Outputs", "output_modularity_raw",   run_id)
  out_root_resid <- file.path("Outputs", "output_modularity_resid", run_id)
  dir.create(out_root_raw,   recursive = TRUE, showWarnings = FALSE)
  dir.create(out_root_resid, recursive = TRUE, showWarnings = FALSE)
  list(run_id = run_id, out_root_raw = out_root_raw, out_root_resid = out_root_resid)
}

# ---------------------------
# 6) Core runner: one test on one coords array
# ---------------------------
run_one_test <- function(coords_arr, test_id, module_list, out_root, coords_label,
                         iter = 999, CI = TRUE, seed = 1,
                         print_progress = FALSE) {

  set.seed(seed)

  test_id <- safe_id(test_id)
  out_dir <- file.path(out_root, test_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  pt_names <- dimnames(coords_arr)[[1]]

  # Validate module names against coords
  assert_points_exist(module_list, pt_names, test_id, coords_label)
  assert_no_dupes_across_modules(module_list, test_id)
  assert_semilandmark_curves_intact(module_list, test_id)

  # Subset to the landmarks used in this test so EVERY landmark is assigned
  pts_in_test <- unique(unlist(module_list, use.names = FALSE))
  coords_arr  <- coords_arr[pts_in_test, , , drop = FALSE]

  # Build partition (named factor) for ALL landmarks in coords_arr
  part <- partition_from_modules(module_list)

  # MANIFEST
  manifest <- data.frame(
    run_id = character(),
    test_id = character(),
    coords_type = character(),
    artifact_type = character(),
    file = character(),
    stringsAsFactors = FALSE
  )
  add_manifest <- function(run_id, artifact_type, file) {
    manifest <<- rbind(manifest, data.frame(
      run_id = run_id,
      test_id = test_id,
      coords_type = coords_label,
      artifact_type = artifact_type,
      file = file,
      stringsAsFactors = FALSE
    ))
  }

  # MODEL CARD
  card_path <- file.path(out_dir, paste0(test_id, "_MODEL_CARD.txt"))
  lines <- c(
    paste0("run_id: ", basename(out_root)),
    paste0("test_id: ", test_id),
    paste0("coords_type: ", coords_label),
    "",
    "module_definitions:"
  )
  for (mn in names(module_list)) {
    lines <- c(lines, paste0("  - ", mn, " (n=", length(module_list[[mn]]), "): ",
                             paste(module_list[[mn]], collapse = ", ")))
  }
  lines <- c(lines,
             "",
             paste0("n_specimens: ", dim(coords_arr)[3]),
             paste0("p_landmarks_used: ", dim(coords_arr)[1]),
             "",
             "sessionInfo():",
             capture.output(sessionInfo())
  )
  write_lines_safe(lines, card_path)
  add_manifest(basename(out_root), "model_card_txt", card_path)

  # -------------- Modularity (CR) --------------
  mod_obj <- geomorph::modularity.test(
    A = coords_arr,
    partition.gp = part,
    iter = iter,
    CI = CI
  )

  mod_txt <- file.path(out_dir, paste0(test_id, "_modularity_test_CR.txt"))
  write_lines_safe(capture.output(print(mod_obj)), mod_txt)
  add_manifest(basename(out_root), "modularity_txt", mod_txt)

  mod_rds <- file.path(out_dir, paste0(test_id, "_modularity_test_CR.rds"))
  saveRDS(mod_obj, mod_rds)
  add_manifest(basename(out_root), "modularity_rds", mod_rds)

  # -------------- Integration (partitioned) --------------
  int_obj <- run_integration_partitioned(coords_arr, part, iter = iter, seed = seed,
                                        print.progress = print_progress)

  int_txt <- file.path(out_dir, paste0(test_id, "_integration_test.txt"))
  write_lines_safe(capture.output(print(int_obj)), int_txt)
  add_manifest(basename(out_root), "integration_txt", int_txt)

  int_rds <- file.path(out_dir, paste0(test_id, "_integration_test.rds"))
  saveRDS(int_obj, int_rds)
  add_manifest(basename(out_root), "integration_rds", int_rds)

  # -------------- Pairwise two.b.pls across module pairs --------------
  mods <- names(module_list)
  if (length(mods) >= 2) {
    pairs <- utils::combn(mods, 2)
    pls_list <- vector("list", ncol(pairs))
    sum_rows <- vector("list", ncol(pairs))

    for (i in seq_len(ncol(pairs))) {
      m1 <- pairs[1, i]
      m2 <- pairs[2, i]

      A1 <- coords_arr[module_list[[m1]], , , drop = FALSE]
      A2 <- coords_arr[module_list[[m2]], , , drop = FALSE]

      pls_obj <- run_two_b_pls_pair(A1 = A1, A2 = A2, iter = iter, seed = seed,
                                    print.progress = print_progress)

      pls_list[[i]] <- pls_obj

      r_val <- tryCatch(as.numeric(pls_obj$r.pls),     error = function(e) NA_real_)
      z_val <- tryCatch(as.numeric(pls_obj$Z),         error = function(e) NA_real_)
      p_val <- tryCatch(as.numeric(pls_obj$P.value),   error = function(e) NA_real_)

      sum_rows[[i]] <- data.frame(
        module1 = m1,
        module2 = m2,
        r_pls = r_val,
        Z = z_val,
        P_value = p_val,
        stringsAsFactors = FALSE
      )
    }

    pls_txt <- file.path(out_dir, paste0(test_id, "_pairwise_two_b_pls.txt"))
    write_lines_safe(
      unlist(lapply(seq_along(pls_list), function(i) {
        hdr <- paste0("\n===== two.b.pls: ", sum_rows[[i]]$module1[1], " vs ", sum_rows[[i]]$module2[1], " =====")
        c(hdr, capture.output(print(pls_list[[i]])))
      })),
      pls_txt
    )
    add_manifest(basename(out_root), "pairwise_pls_txt", pls_txt)

    pls_rds <- file.path(out_dir, paste0(test_id, "_pairwise_two_b_pls.rds"))
    saveRDS(pls_list, pls_rds)
    add_manifest(basename(out_root), "pairwise_pls_rds", pls_rds)

    pls_sum <- do.call(rbind, sum_rows)
    pls_csv <- file.path(out_dir, paste0(test_id, "_pairwise_two_b_pls_summary.csv"))
    write_csv_safe(pls_sum, pls_csv)
    add_manifest(basename(out_root), "pairwise_pls_csv", pls_csv)
  }

  # -------------- MANIFEST --------------
  manifest_path <- file.path(out_dir, "MANIFEST.csv")
  write_csv_safe(manifest, manifest_path)

  invisible(list(out_dir = out_dir, modularity = mod_obj, integration = int_obj))
}

# ---------------------------
# 7) Convenience: run RAW + RESID in one call
# ---------------------------
run_test_raw_resid_1950 <- function(test_id, module_list, iter = 999, CI = TRUE, seed = 1,
                                    print_progress = FALSE) {
  ensure_canonical_1950()
  run <- init_modularity_run()

  run_one_test(coords_1950, test_id, module_list,
               out_root = run$out_root_raw, coords_label = "RAW",
               iter = iter, CI = CI, seed = seed, print_progress = print_progress)

  run_one_test(coords_resid_1950, test_id, module_list,
               out_root = run$out_root_resid, coords_label = "RESID",
               iter = iter, CI = CI, seed = seed, print_progress = print_progress)

  message("DONE. Outputs:")
  message("  RAW:   ", normalizePath(run$out_root_raw, winslash = "/"))
  message("  RESID: ", normalizePath(run$out_root_resid, winslash = "/"))

  invisible(run)
}
