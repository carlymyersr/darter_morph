# ============================================================
# Scripts/tests/TestE_2modules_DV.R
# 1950 Modularity + Integration tests
# Test E: 2 modules (Dorsal–Ventral)
#   - Module 1 (dorsal-ish): cranium_orbital curve + orbit + anterior jaw points (incl premaxilla)
#   - Module 2 (ventral-ish): hyoid_pelvic curve + operculum/procoracoid/pectoral + maxilla
# ============================================================

source("R/02_results/09_modularity_integration/_run_modularity_integration_1950.R")

# ---- Test E modules ----
testE_modules <- list(
  M1_dorsal = c(
    "cranium_orbital_start",
    "cranium_orbital_sl1", "cranium_orbital_sl2", "cranium_orbital_sl3", "cranium_orbital_sl4",
    "cranium_orbital_sl5", "cranium_orbital_sl6", "cranium_orbital_sl7", "cranium_orbital_sl8",
    "cranium_orbital_end",
    "orbit_1", "orbit_2",
    "ab_mandibulae", "preoperculum", "operculum", "premaxilla"
  ),
  M2_ventral = c(
    "hyoid_pelvic_start",
    "hyoid_pelvic_sl1", "hyoid_pelvic_sl2", "hyoid_pelvic_sl3", "hyoid_pelvic_sl4",
    "hyoid_pelvic_sl5", "hyoid_pelvic_sl6", "hyoid_pelvic_sl7", "hyoid_pelvic_sl8",
    "hyoid_pelvic_end",
    "max_curve_preoperculum", "procoracoid",
    "dorsal_pect", "ventral_pect",
    "maxilla"
  )
)

# ---- Run settings ----
ITER <- 999
SEED <- 1
CI   <- TRUE

run_test_raw_resid_1950(
  test_id     = "TestE_2modules_DV",
  module_list = testE_modules,
  iter        = ITER,
  CI          = CI,
  seed        = SEED
)
