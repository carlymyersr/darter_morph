# ============================================================
# Multiple-testing correction for modularity/integration outputs
#
# Reads the timestamped output folders produced by the modularity test scripts:
#   Outputs/output_modularity_raw*
#   Outputs/output_modularity_resid*
#   Outputs/output_modularity_raw_CT3*
#   Outputs/output_modularity_resid_CT3*
#
# Produces Bonferroni and BH corrections for:
#   - modularity.test CR p-values
#   - integration.test p-values
#   - pairwise two.b.pls p-values
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

OUT_DIR <- file.path("Outputs", "modularity_multiple_testing_correction")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

roots <- c(
  "Outputs/output_modularity_raw",
  "Outputs/output_modularity_resid",
  "Outputs/output_modularity_raw_CT3",
  "Outputs/output_modularity_resid_CT3"
)
roots <- roots[dir.exists(roots)]
if (length(roots) == 0) stop("No generated modularity output roots found.")

root_info <- function(root) {
  base <- basename(root)
  data_context <- if (grepl("_CT3$", base)) "CT_time_series" else "habitats_1950"
  shape_type <- if (grepl("resid", base)) "size_corrected_residual" else "raw"
  data.frame(root = root, data_context = data_context, shape_type = shape_type)
}

extract_first_number <- function(pattern, lines) {
  hit <- grep(pattern, lines, value = TRUE)[1]
  if (is.na(hit)) return(NA_real_)
  as.numeric(sub(paste0(".*", pattern, "\\s*:?\\s*"), "", hit))
}

parse_overall_file <- function(path, analysis_type) {
  lines <- readLines(path, warn = FALSE)
  data.frame(
    analysis_type = analysis_type,
    statistic = if (analysis_type == "modularity_CR") {
      extract_first_number("CR", lines)
    } else {
      extract_first_number("r-PLS", lines)
    },
    z = extract_first_number("Effect Size( \\(Z\\))?", lines),
    p_value = extract_first_number("P-value", lines),
    stringsAsFactors = FALSE
  )
}

collect_overall <- function(root) {
  info <- root_info(root)
  files <- c(
    list.files(root, recursive = TRUE, full.names = TRUE, pattern = "_modularity_test_CR[.]txt$"),
    list.files(root, recursive = TRUE, full.names = TRUE, pattern = "_integration_test[.]txt$")
  )
  if (length(files) == 0) return(NULL)
  rows <- lapply(files, function(path) {
    test_id <- basename(dirname(path))
    run_id <- basename(dirname(dirname(path)))
    analysis_type <- if (grepl("_modularity_test_CR[.]txt$", path)) "modularity_CR" else "integration_PLS"
    parsed <- parse_overall_file(path, analysis_type)
    cbind(info, run_id = run_id, test_id = test_id, file = path, parsed)
  })
  do.call(rbind, rows)
}

latest_by_test <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  group_cols <- intersect(
    c("root", "data_context", "shape_type", "analysis_type", "test_id"),
    names(df)
  )
  df %>%
    mutate(file_mtime = as.numeric(file.info(file)$mtime)) %>%
    group_by(across(all_of(group_cols))) %>%
    slice_max(file_mtime, n = 1, with_ties = FALSE) %>%
    ungroup()
}

overall <- latest_by_test(do.call(rbind, lapply(roots, collect_overall)))

if (!is.null(overall) && nrow(overall) > 0) {
  overall_corrected <- overall %>%
    group_by(data_context, shape_type, analysis_type) %>%
    mutate(
      p_bh = p.adjust(p_value, method = "BH"),
      p_bonferroni = p.adjust(p_value, method = "bonferroni")
    ) %>%
    ungroup() %>%
    arrange(data_context, shape_type, analysis_type, test_id)
  write.csv(overall_corrected, file.path(OUT_DIR, "modularity_integration_overall_pvalue_corrections.csv"), row.names = FALSE)
} else {
  overall_corrected <- data.frame()
}

collect_pls <- function(root) {
  info <- root_info(root)
  files <- list.files(root, recursive = TRUE, full.names = TRUE, pattern = "_pairwise_two_b_pls_summary[.]csv$")
  if (length(files) == 0) return(NULL)
  rows <- lapply(files, function(path) {
    test_id <- basename(dirname(path))
    run_id <- basename(dirname(dirname(path)))
    dat <- read.csv(path, stringsAsFactors = FALSE)
    p_col <- intersect(names(dat), c("P_value", "p_value", "Pr...d", "P.value"))[1]
    if (is.na(p_col)) stop("No p-value column found in ", path)
    dat$p_value <- dat[[p_col]]
    cbind(info, run_id = run_id, test_id = test_id, file = path, dat)
  })
  do.call(rbind, rows)
}

pls <- latest_by_test(do.call(rbind, lapply(roots, collect_pls)))

if (!is.null(pls) && nrow(pls) > 0) {
  pls_corrected <- pls %>%
    group_by(data_context, shape_type, test_id) %>%
    mutate(
      p_bh = p.adjust(p_value, method = "BH"),
      p_bonferroni = p.adjust(p_value, method = "bonferroni")
    ) %>%
    ungroup() %>%
    arrange(data_context, shape_type, test_id, module1, module2)
  write.csv(pls_corrected, file.path(OUT_DIR, "modularity_pairwise_pls_pvalue_corrections.csv"), row.names = FALSE)
} else {
  pls_corrected <- data.frame()
}

capture.output(
  {
    cat("Modularity/integration multiple-testing correction\n\n")
    cat("Input roots:\n")
    print(roots)
    cat("\nOverall modularity/integration corrections:\n")
    print(overall_corrected)
    cat("\nPairwise PLS corrections:\n")
    print(pls_corrected)
  },
  file = file.path(OUT_DIR, "modularity_multiple_testing_correction_summary.txt")
)

cat("\nModularity multiple-testing corrections saved to:\n")
cat("  ", normalizePath(OUT_DIR), "\n")
