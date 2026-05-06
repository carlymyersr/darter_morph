# ============================================================
# Scripts/tests_CT3/TestD_2modules_AP_CT3.R
# CT3 Modularity + Integration tests
# Test D: 2 modules (Anterior–Posterior)
#   - Module 1: cranium_orbital + anterior hyoid_pelvic (start + sl1-2)
#   - Module 2: posterior hyoid_pelvic (sl3-8 + end) + orbit + jaw/operculum + pectoral/jaw
# ============================================================

source("R/results/09_modularity_integration/_run_modularity_integration_CT3.R")

# ---- Test D modules ----
testD_modules <- list(
  M1_anterior = c(
    "cranium_orbital_start",
    "cranium_orbital_sl1", "cranium_orbital_sl2", "cranium_orbital_sl3", "cranium_orbital_sl4",
    "cranium_orbital_sl5", "cranium_orbital_sl6", "cranium_orbital_sl7", "cranium_orbital_sl8",
    "cranium_orbital_end"
  ),
  M2_posterior = c(
    "hyoid_pelvic_start",
    "hyoid_pelvic_sl1", "hyoid_pelvic_sl2", "hyoid_pelvic_sl3", "hyoid_pelvic_sl4", "hyoid_pelvic_sl5", "hyoid_pelvic_sl6",
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

run_test_raw_resid_CT3(
  test_id     = "TestD_2modules_AP",
  module_list = testD_modules,
  iter        = ITER,
  CI          = CI,
  seed        = SEED
)