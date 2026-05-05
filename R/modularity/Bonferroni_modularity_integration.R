# ============================================================
# Bonferroni_modularity_integration_clean.R
#
# Reads clean folder structure:
#
# Outputs/Clean_modularity_analysis/
#   1950_RAW/
#     Test_A/
#       TestA_5modules_modularity_test_CR.txt
#       TestA_5modules_integration_test.txt
#     Test_B/
#     Test_C/
#     Test_D/
#     Test_E/
#
#   1950_RESID/
#   CT3_RAW/
#   CT3_RESID/
#
# Each test folder should contain:
#   - one modularity test txt
#   - one integration test txt
#
# Output:
#   Figures/modularity_integration_bonferroni_summary.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

# ------------------------------------------------------------
# 1) Base directory
# ------------------------------------------------------------
base_dir <- file.path("Outputs", "Clean_modularity_analysis")

if (!dir.exists(base_dir)) {
  stop("Base directory not found: ",
       normalizePath(base_dir, winslash = "/", mustWork = FALSE))
}

# ------------------------------------------------------------
# 2) Helpers to extract values from txt files
# ------------------------------------------------------------
extract_p_value <- function(file) {
  txt <- readLines(file, warn = FALSE)
  
  # Only grab the overall p-value line, not pairwise sections
  p_line <- txt[grepl("^P-value:", txt)]
  if (length(p_line) == 0) return(NA_real_)
  
  as.numeric(sub("^P-value:\\s*", "", p_line[1]))
}

extract_statistic_value <- function(file, test_type) {
  txt <- readLines(file, warn = FALSE)
  
  if (test_type == "modularity") {
    line <- txt[grepl("^CR:", txt)]
    if (length(line) == 0) return(NA_real_)
    return(as.numeric(sub("^CR:\\s*", "", line[1])))
  }
  
  if (test_type == "integration") {
    line <- txt[grepl("^r-PLS:", txt)]
    if (length(line) == 0) return(NA_real_)
    return(as.numeric(sub("^r-PLS:\\s*", "", line[1])))
  }
  
  NA_real_
}

extract_effect_size <- function(file, test_type) {
  txt <- readLines(file, warn = FALSE)
  
  if (test_type == "modularity") {
    line <- txt[grepl("^Effect Size:", txt)]
    if (length(line) == 0) return(NA_real_)
    return(as.numeric(sub("^Effect Size:\\s*", "", line[1])))
  }
  
  if (test_type == "integration") {
    line <- txt[grepl("^Effect Size \\(Z\\):", txt)]
    if (length(line) == 0) return(NA_real_)
    return(as.numeric(sub("^Effect Size \\(Z\\):\\s*", "", line[1])))
  }
  
  NA_real_
}

# ------------------------------------------------------------
# 3) Helpers to parse metadata from folder names
# ------------------------------------------------------------
extract_test_label <- function(test_folder_name) {
  case_when(
    str_detect(test_folder_name, "Test[_]?A") ~ "TestA",
    str_detect(test_folder_name, "Test[_]?B") ~ "TestB",
    str_detect(test_folder_name, "Test[_]?C") ~ "TestC",
    str_detect(test_folder_name, "Test[_]?D") ~ "TestD",
    str_detect(test_folder_name, "Test[_]?E") ~ "TestE",
    TRUE ~ NA_character_
  )
}

extract_n_modules <- function(test_folder_name, file_name = NULL) {
  # Prefer file name if available, because it often contains "5modules"
  source_text <- paste(c(test_folder_name, file_name), collapse = " ")
  out <- str_extract(source_text, "\\d+(?=modules)")
  as.numeric(out)
}

extract_hypothesis <- function(test_folder_name) {
  case_when(
    str_detect(test_folder_name, "Test[_]?A") ~ "5 modules (fine partition)",
    str_detect(test_folder_name, "Test[_]?B") ~ "4 modules (curves combined)",
    str_detect(test_folder_name, "Test[_]?C") ~ "3 modules",
    str_detect(test_folder_name, "Test[_]?D") ~ "2 modules (Anterior-Posterior)",
    str_detect(test_folder_name, "Test[_]?E") ~ "2 modules (Dorsal-Ventral)",
    TRUE ~ NA_character_
  )
}

# ------------------------------------------------------------
# 4) Find files flexibly inside each test folder
# ------------------------------------------------------------
find_one_file <- function(folder, pattern, label_for_warning) {
  files <- list.files(folder, pattern = pattern, full.names = TRUE)
  
  if (length(files) == 0) {
    warning("Missing ", label_for_warning, " file in ", folder)
    return(NA_character_)
  }
  
  if (length(files) > 1) {
    warning("Multiple ", label_for_warning, " files found in ",
            folder, ". Using the first one: ", basename(files[1]))
  }
  
  files[1]
}

# ------------------------------------------------------------
# 5) Read one condition folder
# ------------------------------------------------------------
read_condition_results <- function(condition_dir) {
  condition_name <- basename(condition_dir)
  
  parts <- strsplit(condition_name, "_")[[1]]
  if (length(parts) != 2) {
    warning("Skipping unexpected condition folder name: ", condition_name)
    return(NULL)
  }
  
  dataset <- parts[1]
  coords_type <- parts[2]
  
  test_dirs <- list.dirs(condition_dir, recursive = FALSE, full.names = TRUE)
  test_dirs <- test_dirs[grepl("^Test", basename(test_dirs))]
  
  if (length(test_dirs) == 0) {
    warning("No test folders found in: ", condition_dir)
    return(NULL)
  }
  
  out <- lapply(test_dirs, function(td) {
    test_folder <- basename(td)
    
    mod_file <- find_one_file(td, "modularity_test_CR\\.txt$", "modularity")
    int_file <- find_one_file(td, "integration_test\\.txt$", "integration")
    
    # Use whichever file exists to help parse n_modules
    example_file <- basename(c(mod_file, int_file)[!is.na(c(mod_file, int_file))][1])
    
    bind_rows(
      data.frame(
        source_dir = td,
        source_file = if (!is.na(mod_file)) basename(mod_file) else NA_character_,
        dataset = dataset,
        coords_type = coords_type,
        test_folder = test_folder,
        test_label = extract_test_label(test_folder),
        n_modules = extract_n_modules(test_folder, example_file),
        hypothesis = extract_hypothesis(test_folder),
        test_type = "modularity",
        statistic_name = "CR",
        statistic_value = if (!is.na(mod_file)) extract_statistic_value(mod_file, "modularity") else NA_real_,
        effect_name = "Effect Size",
        effect_value = if (!is.na(mod_file)) extract_effect_size(mod_file, "modularity") else NA_real_,
        p_value = if (!is.na(mod_file)) extract_p_value(mod_file) else NA_real_,
        stringsAsFactors = FALSE
      ),
      data.frame(
        source_dir = td,
        source_file = if (!is.na(int_file)) basename(int_file) else NA_character_,
        dataset = dataset,
        coords_type = coords_type,
        test_folder = test_folder,
        test_label = extract_test_label(test_folder),
        n_modules = extract_n_modules(test_folder, example_file),
        hypothesis = extract_hypothesis(test_folder),
        test_type = "integration",
        statistic_name = "r-PLS",
        statistic_value = if (!is.na(int_file)) extract_statistic_value(int_file, "integration") else NA_real_,
        effect_name = "Z",
        effect_value = if (!is.na(int_file)) extract_effect_size(int_file, "integration") else NA_real_,
        p_value = if (!is.na(int_file)) extract_p_value(int_file) else NA_real_,
        stringsAsFactors = FALSE
      )
    )
  })
  
  bind_rows(out)
}

# ------------------------------------------------------------
# 6) Read all condition folders
# ------------------------------------------------------------
condition_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)

results <- lapply(condition_dirs, read_condition_results) %>%
  bind_rows()

if (nrow(results) == 0) {
  stop("No results were found in: ",
       normalizePath(base_dir, winslash = "/", mustWork = FALSE))
}

# ------------------------------------------------------------
# 7) Validate expected structure
# ------------------------------------------------------------
cat("\nRows per family before correction:\n")
print(results %>% count(dataset, coords_type, test_type, name = "n_rows"))

# ------------------------------------------------------------
# 8) Bonferroni correction
#    Family = dataset x coords_type x test_type
# ------------------------------------------------------------
results <- results %>%
  group_by(dataset, coords_type, test_type) %>%
  mutate(
    n_tests = sum(!is.na(p_value)),
    alpha_nominal = 0.05,
    alpha_bonf = ifelse(n_tests > 0, alpha_nominal / n_tests, NA_real_),
    p_value_bonf = p.adjust(p_value, method = "bonferroni"),
    significant_nominal = !is.na(p_value) & p_value < alpha_nominal,
    significant_bonf = !is.na(p_value_bonf) & p_value_bonf < alpha_nominal
  ) %>%
  ungroup()

# ------------------------------------------------------------
# 9) Nice ordering
# ------------------------------------------------------------
results <- results %>%
  mutate(
    dataset = factor(dataset, levels = c("1950", "CT3")),
    coords_type = factor(coords_type, levels = c("RAW", "RESID")),
    test_type = factor(test_type, levels = c("modularity", "integration")),
    test_label = factor(test_label, levels = c("TestA", "TestB", "TestC", "TestD", "TestE"))
  ) %>%
  arrange(dataset, coords_type, test_type, test_label)

# ------------------------------------------------------------
# 10) Save summary
# ------------------------------------------------------------
out_csv <- file.path("Figures", "modularity_integration_bonferroni_summary.csv")
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

write.csv(results, out_csv, row.names = FALSE)

cat("\nSaved Bonferroni summary to:\n  ",
    normalizePath(out_csv, winslash = "/"), "\n", sep = "")

cat("\nPreview:\n")
print(
  results %>%
    select(
      dataset,
      coords_type,
      test_label,
      n_modules,
      hypothesis,
      test_type,
      statistic_name,
      statistic_value,
      effect_name,
      effect_value,
      p_value,
      p_value_bonf,
      alpha_bonf,
      significant_nominal,
      significant_bonf
    ),
  n = nrow(results)
)