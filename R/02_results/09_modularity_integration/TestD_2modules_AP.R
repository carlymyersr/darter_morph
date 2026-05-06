# ============================================================
# Scripts/tests/TestD_2modules_AP.R
# 1950 Modularity + Integration tests
# Test D: 2 modules (Anterior–Posterior)
#   - Module 1: full cranium_orbital curve
#   - Module 2: full hyoid_pelvic curve + orbit + jaw/operculum + pectoral/jaw
# ============================================================

source("R/02_results/09_modularity_integration/_run_modularity_integration_1950.R")

# ---- Test D modules ----
testD_modules <- list(
  M1_anterior = c(
    "cranium_orbital_start",
    "cranium_orbital_sl1", "cranium_orbital_sl2", "cranium_orbital_sl3", "cranium_orbital_sl4",
    "cranium_orbital_sl5", "cranium_orbital_sl6", "cranium_orbital_sl7", "cranium_orbital_sl8",
    "cranium_orbital_end"
  ),
  M2_posterior = c(
    "hyoid_pelvic_start","hyoid_pelvic_sl1","hyoid_pelvic_sl2","hyoid_pelvic_sl3", "hyoid_pelvic_sl4", "hyoid_pelvic_sl5", "hyoid_pelvic_sl6",
    "hyoid_pelvic_sl7", "hyoid_pelvic_sl8",
    "hyoid_pelvic_end",
    "orbit_1", "orbit_2",
    "ab_mandibulae", "preoperculum", "max_curve_preoperculum", "operculum",
    "dorsal_pect", "ventral_pect", "procoracoid", "premaxilla", "maxilla"
  )
)

# ---- Run settings ----
ITER <- 999
SEED <- 1
CI   <- TRUE

run_test_raw_resid_1950(
  test_id     = "TestD_2modules_AP",
  module_list = testD_modules,
  iter        = ITER,
  CI          = CI,
  seed        = SEED
)
