# =============================================================================
# JULES Model Configuration File
# Description: Centralized configuration for JULES model calibration and analysis
# =============================================================================

# =============================================================================
# MAIN CONFIGURATION
# =============================================================================
CONFIG <- list(
  # === PROJECT SETTINGS ===
  project_name = "JULES_Emulation_Clean",
  version = "1.0",
  
  # === SCENARIO SETTINGS ===
  rcp = '85',                          # RCP scenario: '26', '45', '60', '85'
  ensemble = '15',                     # Ensemble ID: '01', '04', '06', '15'
  
  # === FILE PATHS (RELATIVE TO PROJECT ROOT) ===
  base_folder = "/Users/mingda/Downloads/Jules-dev/jules_clean",
  data_folder = "data",
  results_folder = "results", 
  lib_folder = "lib",
  
  # === DATA FILES ===
  design_data_pattern = "design_data_rcp_{rcp}_{ensemble}.rds",  # Will be filled with rcp/ensemble
  calibration_ensemble_file = "calibration_ensembleREDHARP.csv",
  base_ensemble_file = "calibration_ensembleAnna.csv",
  
  # === MODEL PARAMETERS ===
  variables = c("oneoveralpha", "hcon", "vcrit", "slope", "mean_pr1", 
               "spring_tas1", "winter_tasmin1", "mean_rsds_dec1", "vpd1"),
  transforms = c("identity", "identity", "identity", "log10", "log10", 
                "identity", "identity", "identity", "identity"),
  
  # === SAMPLING PARAMETERS ===
  n_samples = 400,                     # Number of design points
  n_test_samples = 100,                # Number of test samples
  random_seed = 50,                    # Random seed for reproducibility
  
  # === TIME PARAMETERS ===
  min_year = 2025,
  max_year = 2050,
  n_time_points = 10,                  # Number of time points for time-varying analysis
  
  # === EMULATION PARAMETERS ===
  pca_components = 1,                  # Number of PCA components
  emulator_names = c("matern2.5", "sexp"), # Kernel names for DGP
  use_vecchia = TRUE,                  # Use Vecchia approximation
  validation_seed = 99,                # Seed for emulator validation
  
  # === OUTPUT SETTINGS ===
  save_intermediate = TRUE,            # Save intermediate results
  create_plots = TRUE,                 # Generate diagnostic plots
  verbose = TRUE,                      # Enable verbose logging
  
  # === COORDINATE SETTINGS ===
  default_x = 352500,
  default_y = 497500,
  
  # === ADVANCED SETTINGS ===
  use_time_varying = TRUE,             # Enable time-varying emulation
  perform_sequential_design = FALSE,   # Enable sequential design (advanced)
  max_sequential_waves = 3,            # Maximum sequential design waves
  
  # === COMPUTATIONAL SETTINGS ===
  use_parallel = FALSE,                # Enable parallel processing
  n_cores = NULL,                      # Number of cores (NULL = auto-detect)
  chunk_size = 460                     # Chunk size for predictions
)

# =============================================================================
# PARAMETER RANGES (EXAMPLE - USER MUST DEFINE BASED ON THEIR STUDY)
# =============================================================================
PARAMETER_RANGES <- data.frame(
  # Define parameter ranges here - these are examples
  # mort_base_nt = c(0.001, 0.01),      # min, max values
  # alpha_red_nt = c(0.005, 0.02),      # min, max values  
  # crwn_area0_nt = c(0.1, 0.2),        # min, max values
  # lai_bal0_nt = c(0.8, 1.2),          # min, max values
  # mass0 = c(0.05, 0.15)               # min, max values (optional)
  stringsAsFactors = FALSE
)
# rownames(PARAMETER_RANGES) <- c('min', 'max')

# =============================================================================
# DERIVED SETTINGS (AUTO-CALCULATED FROM CONFIG)
# =============================================================================
CONFIG$id <- paste0(CONFIG$rcp, '_', CONFIG$ensemble)
CONFIG$design_data_file <- gsub("\\{rcp\\}", CONFIG$rcp, 
                               gsub("\\{ensemble\\}", CONFIG$ensemble, 
                                   CONFIG$design_data_pattern))

# Transform function mapping
TRANSFORM_FUNCTIONS <- list(
  "identity" = identity,
  "log10" = log10,
  "sqrt" = sqrt,
  "log" = log
)

CONFIG$transform_functions <- lapply(CONFIG$transforms, function(x) TRANSFORM_FUNCTIONS[[x]])

# =============================================================================
# VALIDATION FUNCTIONS
# =============================================================================

#' Validate configuration settings
validate_config <- function() {
  errors <- c()
  warnings <- c()
  
  # Check required directories
  required_dirs <- c(
    file.path(CONFIG$base_folder, CONFIG$data_folder),
    file.path(CONFIG$base_folder, CONFIG$results_folder),
    file.path(CONFIG$base_folder, CONFIG$lib_folder)
  )
  
  for (dir in required_dirs) {
    if (!dir.exists(dir)) {
      warnings <- c(warnings, paste("Directory does not exist:", dir))
    }
  }
  
  # Check parameter consistency
  if (length(CONFIG$variables) != length(CONFIG$transforms)) {
    errors <- c(errors, "Length of variables and transforms must match")
  }
  
  # Check RCP/ensemble values
  valid_rcps <- c("26", "45", "60", "85")
  valid_ensembles <- c("01", "04", "06", "15")
  
  if (!(CONFIG$rcp %in% valid_rcps)) {
    warnings <- c(warnings, paste("Unusual RCP value:", CONFIG$rcp))
  }
  
  if (!(CONFIG$ensemble %in% valid_ensembles)) {
    warnings <- c(warnings, paste("Unusual ensemble value:", CONFIG$ensemble))
  }
  
  # Report results
  if (length(errors) > 0) {
    cat("=== CONFIGURATION ERRORS ===\n")
    for (error in errors) {
      cat("ERROR:", error, "\n")
    }
    stop("Configuration validation failed")
  }
  
  if (length(warnings) > 0) {
    cat("=== CONFIGURATION WARNINGS ===\n")
    for (warning in warnings) {
      cat("WARNING:", warning, "\n")
    }
  }
  
  cat("Configuration validation completed\n")
  return(TRUE)
}

#' Print configuration summary
print_config_summary <- function() {
  cat("=== JULES CONFIGURATION SUMMARY ===\n")
  cat("Project:", CONFIG$project_name, "v", CONFIG$version, "\n")
  cat("Scenario: RCP", CONFIG$rcp, "Ensemble", CONFIG$ensemble, "\n")
  cat("ID:", CONFIG$id, "\n")
  cat("Variables:", length(CONFIG$variables), "variables\n")
  cat("Samples:", CONFIG$n_samples, "design points\n")
  cat("Time range:", CONFIG$min_year, "-", CONFIG$max_year, "\n")
  cat("Base folder:", CONFIG$base_folder, "\n")
  cat("=====================================\n")
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Get full path for data file
get_data_path <- function(filename) {
  return(file.path(CONFIG$base_folder, CONFIG$data_folder, filename))
}

#' Get full path for results file
get_results_path <- function(filename) {
  return(file.path(CONFIG$base_folder, CONFIG$results_folder, filename))
}

#' Get full path for lib file  
get_lib_path <- function(filename) {
  return(file.path(CONFIG$base_folder, CONFIG$lib_folder, filename))
}

#' Log message with timestamp
log_message <- function(message, level = "INFO") {
  if (CONFIG$verbose) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s: %s\n", timestamp, level, message))
  }
}

# =============================================================================
# INITIALIZATION MESSAGE
# =============================================================================
if (CONFIG$verbose) {
  cat("JULES Configuration loaded successfully\n")
  cat("Run print_config_summary() to see current settings\n")
  cat("Run validate_config() to validate configuration\n")
} 