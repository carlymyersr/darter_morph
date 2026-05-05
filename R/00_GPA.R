# ============================================================
# 00_GPA.R
# Core morphometric initialization for darter dataset
# ============================================================


# ============================================================
# Determine project root automatically (source-safe)
#   - If sourced, sets wd to project root (parent of /R/)
#   - If run line-by-line, does NOT error (just leaves wd alone)
# ============================================================

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
# User flags
# ---------------------------
DO_SANITY_PLOTS <- FALSE   # TRUE = draw diagnostic plots
VERBOSE         <- TRUE    # FALSE = silence most cat() output

if (VERBOSE) cat("Running 00_setup_morpho.R...\n")

# ---------------------------
# Required libraries
# ---------------------------
suppressPackageStartupMessages({
  library(StereoMorph)
  library(geomorph)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(tidyr)
  library(stringr)
})

# ============================================================
# 1) Inputs
# ============================================================

img_dir       <- "photos"
shapes_dir    <- "side_shapes"
lm_ref_file   <- "landmarks_ref.txt"

# MANDATORY drops for all downstream analyses
#these are body landmarks that are not to be included 
DROP_POINTS <- c("1","2","3","4","5","6","7")

# ============================================================
# 2) Curve definitions
# ============================================================

curves_file <- "darter_curves.txt"

curves_ref <- as.matrix(
  read.table(
    curves_file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
)

stopifnot(ncol(curves_ref) == 3)
colnames(curves_ref) <- c("curve_name", "start", "end")

# ============================================================
# 6) trim whitespace in landmark names (prevents "premaxilla ")
#    - handles Unicode whitespace (incl NBSP) robustly
#    NOTE: Defined early so it can be used for curve refs + template
# ============================================================

trim_names <- function(x) {
  x <- as.character(x)
  
  # normalize encoding so regex sees characters consistently
  x <- enc2utf8(x)
  
  # convert any Unicode "separator" spaces (Zs) and any whitespace to a normal space
  # (NBSP is in Zs)
  x <- gsub("\\p{Zs}+", " ", x, perl = TRUE)
  x <- gsub("[[:space:]]+", " ", x)
  
  # hard trim leading/trailing whitespace
  x <- sub("^\\s+", "", x, perl = TRUE)
  x <- sub("\\s+$", "", x, perl = TRUE)
  
  x
}

trim_dimnames_landmarks <- function(coords_arr) {
  dn <- dimnames(coords_arr)
  dn[[1]] <- trim_names(dn[[1]])
  dimnames(coords_arr) <- dn
  coords_arr
}

# Clean curve ref inputs
curves_ref[, "curve_name"] <- trim_names(curves_ref[, "curve_name"])
curves_ref[, "start"]      <- trim_names(curves_ref[, "start"])
curves_ref[, "end"]        <- trim_names(curves_ref[, "end"])

stopifnot(!anyDuplicated(curves_ref[, "curve_name"]))

curve_endpoints <- setNames(
  lapply(seq_len(nrow(curves_ref)), function(i) curves_ref[i, c("start", "end")]),
  curves_ref[, "curve_name"]
)

# ---- Semilandmark density (TOTAL per curve incl endpoints) ----
# If you want 10 semilandmarks BETWEEN endpoints, set total to 12.
# Right now total=10 => internal=8.
n_curve_total <- c(
  hyoid_pelvic    = 10,
  cranium_orbital = 10
)

missing_counts <- setdiff(curves_ref[, "curve_name"], names(n_curve_total))
extra_counts   <- setdiff(names(n_curve_total), curves_ref[, "curve_name"])
if (length(missing_counts) > 0) stop("Missing n_curve_total for: ", paste(missing_counts, collapse=", "))
if (length(extra_counts)   > 0) stop("Extra curves in n_curve_total: ", paste(extra_counts, collapse=", "))

n_curve_total <- n_curve_total[curves_ref[, "curve_name"]]
n_internal    <- n_curve_total - 2
stopifnot(all(n_internal >= 0))

nCurvePts <- unname(n_curve_total)

# ============================================================
# 3) Read shapes + ID sanity checks
# ============================================================

shapes <- readShapes(shapes_dir)

stopifnot(
  !is.null(shapes$landmarks.scaled),
  !is.null(shapes$curves.scaled),
  !is.null(shapes$scaling)
)

lm_ids <- dimnames(shapes$landmarks.scaled)[[3]]
cv_ids <- names(shapes$curves.scaled)
sc_ids <- names(shapes$scaling)

ids <- Reduce(intersect, list(lm_ids, cv_ids, sc_ids))
stopifnot(length(ids) > 0)

# Fixed landmark names from template (trimmed!)
lm_names_fixed <- trim_names(dimnames(shapes$landmarks.scaled)[[1]])
stopifnot(dim(shapes$landmarks.scaled)[2] == 2)

endpoints_all <- trim_names(unlist(curve_endpoints, use.names = FALSE))
missing_endpoints <- setdiff(endpoints_all, lm_names_fixed)
if (length(missing_endpoints) > 0) {
  stop("Curve endpoints missing from template: ",
       paste(missing_endpoints, collapse=", "))
}

if (VERBOSE) cat("Specimens with all data:", length(ids), "\n")

# ============================================================
# 4) Helper functions
# ============================================================

resample_curve <- function(M, n) {
  stopifnot(is.matrix(M), ncol(M) == 2, n >= 2)
  if (nrow(M) == n) return(M)
  d <- sqrt(rowSums((M[-1, , drop=FALSE] - M[-nrow(M), , drop=FALSE])^2))
  s <- c(0, cumsum(d))
  if (max(s) == 0) stop("Curve has zero length (all points identical).")
  t <- seq(0, max(s), length.out = n)
  x <- approx(s, M[,1], xout = t)$y
  y <- approx(s, M[,2], xout = t)$y
  cbind(x, y)
}

make_curve_index <- function(curve_name, endpoints, pt_names) {
  endpoints  <- trim_names(endpoints)
  curve_name <- trim_names(curve_name)
  
  start <- match(endpoints[1], pt_names)
  end   <- match(endpoints[2], pt_names)
  stopifnot(!is.na(start), !is.na(end))
  
  sl_names <- grep(paste0("^", curve_name, "_sl[0-9]+$"), pt_names, value=TRUE)
  if (length(sl_names) > 0) {
    sl_names <- sl_names[order(as.integer(sub(paste0("^",curve_name,"_sl"),"",sl_names)))]
    sl_idx <- match(sl_names, pt_names)
    stopifnot(all(!is.na(sl_idx)))
    return(c(start, sl_idx, end))
  }
  
  c(start, end)
}

make_sliders_from_index <- function(idx_seq) {
  m <- length(idx_seq)
  if (m < 3) return(matrix(numeric(0), ncol=3))
  out <- cbind(idx_seq[1:(m-2)], idx_seq[2:(m-1)], idx_seq[3:m])
  colnames(out) <- c("before","slide","after")
  out
}

make_links_from_index <- function(idx_seq) {
  cbind(idx_seq[-length(idx_seq)], idx_seq[-1])
}

# ============================================================
# 5) Build combined coordinate array
# ============================================================

pt_names <- lm_names_fixed
for (cn in curves_ref[, "curve_name"]) {
  if (n_internal[cn] > 0) {
    pt_names <- c(pt_names, paste0(cn, "_sl", seq_len(n_internal[cn])))
  }
}

coords_all <- array(
  NA_real_,
  dim = c(length(pt_names), 2, length(ids)),
  dimnames = list(pt_names, c("x","y"), ids)
)

for (id in ids) {
  
  LM <- shapes$landmarks.scaled[, , id]
  rownames(LM) <- lm_names_fixed
  coords_all[lm_names_fixed, , id] <- LM
  
  for (cn in curves_ref[, "curve_name"]) {
    M <- shapes$curves.scaled[[id]][[cn]]
    if (is.null(M)) stop(sprintf("Specimen %s missing curve '%s'", id, cn))
    
    R <- resample_curve(M, n_curve_total[cn])
    
    # Anchor endpoints to fixed landmarks
    R[1,]       <- coords_all[curve_endpoints[[cn]][1], , id]
    R[nrow(R),] <- coords_all[curve_endpoints[[cn]][2], , id]
    
    if (n_internal[cn] > 0) {
      rows <- paste0(cn, "_sl", seq_len(n_internal[cn]))
      coords_all[rows, , id] <- R[2:(nrow(R)-1), , drop=FALSE]
    }
  }
}

stopifnot(!any(is.na(coords_all)))

if (DO_SANITY_PLOTS) {
  id0 <- ids[1]
  plot(coords_all[, , id0], asp=1, pch=16, main=paste("Pre-DROP coords_all:", id0))
}


# ============================================================
# 6b) Trim dimnames (coords_all) + FAIL FAST if any trailing WS
# ============================================================

coords_all <- trim_dimnames_landmarks(coords_all)
stopifnot(!any(grepl("\\s$", dimnames(coords_all)[[1]], perl = TRUE)))

# ============================================================
# 7) Drop low-confidence FIXED landmarks (MANDATORY)
# ============================================================

pt_names_now <- dimnames(coords_all)[[1]]

missing_drop <- setdiff(DROP_POINTS, pt_names_now)
if (length(missing_drop) > 0) {
  stop("Requested drop point(s) not found in coords_all: ",
       paste(missing_drop, collapse = ", "))
}

bad_endpoints <- intersect(DROP_POINTS, endpoints_all)
if (length(bad_endpoints) > 0) {
  stop("Do not drop curve endpoints (required for anchoring/sliding): ",
       paste(bad_endpoints, collapse = ", "))
}

if (any(grepl("_sl", DROP_POINTS))) {
  stop("Do not drop *_sl semilandmarks here. Only drop FIXED landmarks.")
}

keep_pts <- setdiff(pt_names_now, DROP_POINTS)

coords_all <- coords_all[keep_pts, , , drop = FALSE]
dimnames(coords_all)[[1]] <- keep_pts

# Canonical kept point names (AFTER drop)
KEEP_POINTS <- dimnames(coords_all)[[1]]

# Canonical fixed landmark names (AFTER trimming + drop)
lm_names <- intersect(lm_names_fixed, KEEP_POINTS)
lm_names <- trim_names(lm_names)
stopifnot(!any(grepl("\\s$", lm_names, perl = TRUE)))

if (VERBOSE) {
  cat("Dropped fixed landmarks:", paste(DROP_POINTS, collapse = ", "), "\n")
  cat("Number of kept landmarks:", length(KEEP_POINTS), "\n")
  cat("coords_all new dim:", paste(dim(coords_all), collapse = " x "), "\n")
}

stopifnot(!any(is.na(coords_all)))

if (DO_SANITY_PLOTS) {
  plot(coords_all[, , 1], asp = 1, pch = 16,
       main = "Pre-GPA coords_all after DROP_POINTS (specimen 1)")
}

# ============================================================
# 8) Sliders + links (built AFTER drops)
# ============================================================

pt_names_now <- dimnames(coords_all)[[1]]
p_now <- length(pt_names_now)

# Endpoints must still exist
stopifnot(all(endpoints_all %in% pt_names_now))

sliders_list <- lapply(curves_ref[, "curve_name"], function(cn) {
  idx_seq <- make_curve_index(cn, curve_endpoints[[cn]], pt_names_now)
  make_sliders_from_index(idx_seq)
})
curves_sliders <- do.call(rbind, sliders_list)

links_list <- lapply(curves_ref[, "curve_name"], function(cn) {
  idx_seq <- make_curve_index(cn, curve_endpoints[[cn]], pt_names_now)
  make_links_from_index(idx_seq)
})
links_all <- do.call(rbind, links_list)

# Basic sanity
stopifnot(is.matrix(curves_sliders), ncol(curves_sliders) == 3)
stopifnot(all(curves_sliders >= 1), all(curves_sliders <= p_now))
stopifnot(max(links_all) <= p_now)

# Ensure slide column are semis only
slide_names <- pt_names_now[curves_sliders[, 2]]
stopifnot(all(grepl("_sl", slide_names)))
stopifnot(!any(slide_names %in% endpoints_all))

# ============================================================
# 9) Sliding GPA
# ============================================================

gpa <- gpagen(coords_all, curves = curves_sliders, PrinAxes = FALSE)

# Restore dimnames explicitly (geomorph sometimes drops them)
dimnames(gpa$coords)[[1]] <- dimnames(coords_all)[[1]]
dimnames(gpa$coords)[[2]] <- c("x","y")
dimnames(gpa$coords)[[3]] <- dimnames(coords_all)[[3]]

coords_gpa <- gpa$coords

# Final trim + FAIL FAST
coords_gpa <- trim_dimnames_landmarks(coords_gpa)
stopifnot(!any(grepl("\\s$", dimnames(coords_gpa)[[1]], perl = TRUE)))

# Keep lm_names synced to the canonical coords_gpa (no stale objects!)
lm_names <- dimnames(coords_gpa)[[1]]
lm_names <- lm_names[!grepl("_sl[0-9]+$", lm_names)]
lm_names <- setdiff(lm_names, DROP_POINTS)
lm_names <- trim_names(lm_names)
stopifnot(!any(grepl("\\s$", lm_names, perl = TRUE)))

if (DO_SANITY_PLOTS) {
  plot(coords_gpa[, , 1], asp=1, pch=16, main="Post-GPA (specimen 1)")
}

# ============================================================
# 10) Canonical object list
# ============================================================

MORPHO_OBJECTS <- c(
  "curves_ref",
  "curve_endpoints",
  "endpoints_all",
  "n_curve_total",
  "n_internal",
  "nCurvePts",
  "shapes",
  "ids",
  "lm_names",
  "DROP_POINTS",
  "KEEP_POINTS",
  "coords_all",
  "curves_sliders",
  "links_all",
  "gpa",
  "coords_gpa"
)

if (VERBOSE) {
  cat("\nObjects created by 00_setup_morpho.R:\n")
  print(MORPHO_OBJECTS)
  cat("\nSetup complete.\n")
}
