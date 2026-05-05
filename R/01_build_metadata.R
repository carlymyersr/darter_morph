# ============================================================
# R/01_metadata.R
# Build ALL-specimen metadata aligned to coords_gpa order
#
# Requires (from R/00_setup_morpho.R):
#   coords_gpa, gpa
#
# Outputs:
#   gdf, hab_palette, collection_lookup, ids_gpa
#   + METADATA_OBJECTS (character vector of object names created)
# ============================================================

# ---------------------------
# Flags (inherits from 00_GPA.R if already defined)
# ---------------------------
if (!exists("DO_SANITY_PLOTS")) DO_SANITY_PLOTS <- FALSE
if (!exists("VERBOSE"))         VERBOSE         <- TRUE

if (VERBOSE) cat("Running R/01_GPA.R...\n")

# ---------------------------
# Preconditions
# ---------------------------
if (!exists("coords_gpa")) stop("coords_gpa not found. Run source('R/00_setup_morpho.R') first.")
if (!exists("gpa"))        stop("gpa not found. Run source('R/00_setup_morpho.R') first.")


suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

patch_rrpp_parallel_setup <- function(fallback = 2L) {
  if (!requireNamespace("RRPP", quietly = TRUE)) return(invisible(FALSE))

  ns <- asNamespace("RRPP")

  safe_cores <- function() {
    n <- tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
    if (length(n) == 0 || is.na(n) || n < 2) fallback else n
  }

  safe_parallel_setup <- function(Parallel) {
    ParLog <- is.logical(Parallel)
    usecluster <- inherits(Parallel, "cluster")
    if (usecluster) {
      cluster <- Parallel
      ParLog <- Parallel <- TRUE
    } else {
      cluster <- FALSE
    }

    ParCores <- NULL
    if (is.numeric(Parallel)) {
      ParCores <- Parallel
      ParLog <- TRUE
      Parallel <- TRUE
    }

    n_cores <- safe_cores()
    if (ParLog && is.null(ParCores)) {
      ParCores <- n_cores - 1
      if (usecluster) ParCores <- length(cluster)
      if (!Parallel) ParCores <- 1
    }
    if (is.numeric(ParCores) && ParCores > n_cores - 1) {
      ParCores <- n_cores - 1
    }

    Unix <- .Platform$OS.type == "unix"
    forking <- Unix && !usecluster
    if (is.null(ParCores)) ParCores <- 1
    if (ParCores == 1) {
      ParLog <- FALSE
      forking <- FALSE
      usecluster <- FALSE
      cluster <- NULL
    }
    if (ParCores > 1) {
      if (!Unix && !usecluster) cluster <- parallel::makeCluster(ParCores)
      if (Unix && usecluster) Unix <- FALSE
    }

    list(
      Parallel = ParLog,
      Unix = Unix,
      forking = forking,
      ParCores = ParCores,
      usecluster = usecluster,
      cluster = cluster
    )
  }

  was_locked <- bindingIsLocked("Parallel.setup", ns)
  if (was_locked) unlockBinding("Parallel.setup", ns)
  assign("Parallel.setup", safe_parallel_setup, envir = ns)
  if (was_locked) lockBinding("Parallel.setup", ns)

  invisible(TRUE)
}

patch_rrpp_parallel_setup()

# ============================================================
# 10) Build metadata aligned to coords_gpa order (ALL specimens)
# Outputs:
#   gdf
# Notes:
#   - Centroid size (Csize) is ALWAYS available from gpagen()
#   - Default allometry size variable is logCsize
# ============================================================

ids_gpa <- dimnames(coords_gpa)[[3]]

# Always-available size from GPA
gdf <- data.frame(
  specimen = ids_gpa,
  Csize    = as.numeric(gpa$Csize),
  logCsize = log(as.numeric(gpa$Csize)),
  stringsAsFactors = FALSE
)

stopifnot(identical(gdf$specimen, ids_gpa))

gdf$SL_mm <- NA_real_
gdf$logSL <- NA_real_

# ---- Canonical choice for allometry correction ----
# This is what you should use in procD.lm shape ~ size
gdf$size_for_allometry <- gdf$logCsize
gdf$size_label <- "logCsize (centroid size from GPA)"

# ============================================================
# Helper: Allometry correction via procD.lm (robust evaluation)
#   Regress shape ~ size
#   Returns residual shape array aligned to input coords
# ============================================================

allometry_residuals <- function(coords_arr, size_vec) {
  
  stopifnot(length(dim(coords_arr)) == 3)
  stopifnot(length(size_vec) == dim(coords_arr)[3])
  
  # enforce numeric
  size_vec <- as.numeric(size_vec)
  
  # Ensure specimen order consistency (if named)
  if (!is.null(names(size_vec))) {
    stopifnot(identical(names(size_vec), dimnames(coords_arr)[[3]]))
  }
  
  # IMPORTANT: procD.lm often wants predictors in `data=`
  dat <- data.frame(size = size_vec)
  
  fit <- geomorph::procD.lm(coords_arr ~ size, data = dat)
  
  coords_resid <- geomorph::arrayspecs(
    residuals(fit),
    p = dim(coords_arr)[1],
    k = dim(coords_arr)[2]
  )
  
  dimnames(coords_resid) <- dimnames(coords_arr)
  
  return(list(
    residuals = coords_resid,
    fit       = fit
  ))
}

# ============================================================
# 11) Add grouping variables (habitat, year)
# ============================================================

gdf$habitat <- dplyr::case_when(
  grepl("^F2380_|^F2383_|^F2381_|^F2423_|^F2466_", gdf$specimen) ~ "Connecticut River",
  grepl("^F2406_|^F2388_",                         gdf$specimen) ~ "Sawmill River",
  grepl("^F2384_",                                 gdf$specimen) ~ "Swift River",
  grepl("^F2386_|^F2377_|^F2382_",                 gdf$specimen) ~ "Fort River",
  grepl("^F2374_|^F2379_|^F2385_",                 gdf$specimen) ~ "Quabbin",
  TRUE                                                           ~ NA_character_
)

gdf$year <- dplyr::case_when(
  grepl("^F2466_", gdf$specimen) ~ 1970,
  grepl("^F2423_", gdf$specimen) ~ 1956,
  grepl("^F2406_", gdf$specimen) ~ 1954,
  grepl("^F2376_", gdf$specimen) ~ 1948,  # Mill River specimens will be NA habitat anyway
  TRUE                           ~ 1950
)


# Global habitat palette (used everywhere)
hab_palette <- c(
  "Connecticut River" = "steelblue",
  "Sawmill River"     = "orchid3",
  "Swift River"       = "tomato",
  "Fort River"        = "black",
  "Quabbin"           = "darkgoldenrod2"
)

# ============================================================
# 11b) Collection ID + month/date lookup
# ============================================================

gdf$collection_id <- sub("^([A-Za-z]\\d+)_.*$", "\\1", gdf$specimen)

if (VERBOSE) {
  cat("Unique collection_id count:", length(unique(gdf$collection_id)), "\n")
}

collection_lookup <- tibble::tribble(
  ~collection_id, ~date_str,    ~month,
  "F2384",        "8/16/1950",  8,
  "F2377",        "6/23/1950",  6,
  "F2382",        "8/5/1950",   8,
  "F2386",        "10/21/1950", 10,
  "F2380",        "6/6/1950",   6,
  "F2381",        "7/12/1950",  7,
  "F2383",        "8/7/1950",   8,
  "F2388",        "11/25/1950", 11,
  "F2374",        "5/23/1950",  5,
  "F2379",        "6/27/1950",  6,
  "F2385",        "8/12/1950",  8,
  "F2423",        "10/21/1956", 10,
  "F2466",        "9/22/1970",  9,
) %>%
  dplyr::mutate(date = as.Date(date_str, format = "%m/%d/%Y"))

gdf <- gdf %>%
  dplyr::left_join(
    collection_lookup %>% dplyr::select(collection_id, month, date),
    by = "collection_id"
  )

if (VERBOSE) {
  cat("\nMonth counts (including NA):\n")
  print(table(gdf$month, useNA = "ifany"))
  
  cat("\nCollection IDs missing month info:\n")
  print(sort(unique(gdf$collection_id[is.na(gdf$month)])))
}

#Collection IDs missing month info:
# "F2376" "F2406" "F2485"
# These are collections not used in this analysis 


# ============================================================
# Canonical object list (setup-style)
# ============================================================

METADATA_OBJECTS <- c(
  "ids_gpa",
  "gdf",
  "hab_palette",
  "collection_lookup"
)

if (VERBOSE) {
  cat("\nObjects created by R/01_build_metadata.R:\n")
  print(METADATA_OBJECTS)
  cat("\nMetadata build complete.\n")
}
