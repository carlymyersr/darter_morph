# ============================================================
# Scripts/tests_CT3/TestA_5modules_CT3.R
# CT3 Modularity + Integration tests
# Test A: 5 modules (as specified)
# ============================================================

source("R/modularity/_run_modularity_integration_CT3.R")

# ---- Test A modules ----
testA_modules <- list(
  M1_cranium_orbital = c(
    "cranium_orbital_start",
    "cranium_orbital_sl1", "cranium_orbital_sl2", "cranium_orbital_sl3", "cranium_orbital_sl4",
    "cranium_orbital_sl5", "cranium_orbital_sl6", "cranium_orbital_sl7", "cranium_orbital_sl8",
    "cranium_orbital_end"
  ),
  M2_hyoid_pelvic = c(
    "hyoid_pelvic_start",
    "hyoid_pelvic_sl1", "hyoid_pelvic_sl2", "hyoid_pelvic_sl3", "hyoid_pelvic_sl4",
    "hyoid_pelvic_sl5", "hyoid_pelvic_sl6", "hyoid_pelvic_sl7", "hyoid_pelvic_sl8",
    "hyoid_pelvic_end"
  ),
  M3_orbit = c("orbit_1", "orbit_2"),
  M4_jaw_operculum = c("ab_mandibulae", "preoperculum", "max_curve_preoperculum", "operculum"),
  M5_pectoral_jaw = c("dorsal_pect", "ventral_pect", "procoracoid", "premaxilla", "maxilla")
)

# ---- Run settings ----
ITER <- 999
SEED <- 1
CI   <- TRUE

run_test_raw_resid_CT3(
  test_id     = "TestA_5modules",
  module_list = testA_modules,
  iter        = ITER,
  CI          = CI,
  seed        = SEED
)