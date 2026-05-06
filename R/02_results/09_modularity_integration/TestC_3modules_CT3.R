# ============================================================
# Scripts/tests_CT3/TestC_3modules_CT3.R
# CT3 Modularity + Integration tests
# Test C: 3 modules
#   - All hyoid_pelvic curve 
#   - All upper jaw 
#   - everything else
# ============================================================

source("R/02_results/09_modularity_integration/_run_modularity_integration_CT3.R")

# ---- Test C modules ----
testC_modules <- list(
  M1_cranium_orbital_plus_hp_anterior = c(
    "cranium_orbital_start",
    "cranium_orbital_sl1", "cranium_orbital_sl2", "cranium_orbital_sl3", "cranium_orbital_sl4",
    "cranium_orbital_sl5", "cranium_orbital_sl6", "cranium_orbital_sl7", "cranium_orbital_sl8",
    "cranium_orbital_end"
  ),
  M2_hp_posterior_plus_jaw_pectoral = c(
    "hyoid_pelvic_start",
    "hyoid_pelvic_sl1", "hyoid_pelvic_sl2","hyoid_pelvic_sl3", "hyoid_pelvic_sl4", "hyoid_pelvic_sl5", "hyoid_pelvic_sl6",
    "hyoid_pelvic_sl7", "hyoid_pelvic_sl8",
    "hyoid_pelvic_end"
  ),
  M3_orbit = c(
    "orbit_1", "orbit_2", "ab_mandibulae", "preoperculum", "max_curve_preoperculum", "operculum",
    "dorsal_pect", "ventral_pect", "procoracoid", "premaxilla", "maxilla"
  )
)

# ---- Run settings ----
ITER <- 999
SEED <- 1
CI   <- TRUE

run_test_raw_resid_CT3(
  test_id     = "TestC_3modules",
  module_list = testC_modules,
  iter        = ITER,
  CI          = CI,
  seed        = SEED
)