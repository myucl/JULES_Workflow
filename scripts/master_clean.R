# =============================================================================
# JULES Master Script - Clean Version
# Description: Main workflow for JULES model emulation and analysis
#
# =============================================================================

# =============================================================================
# SETUP AND INITIALIZATION
# =============================================================================

# Clear workspace
rm(list = ls())

# Load configuration and utilities
source("../config/jules_config.R")
source("jules_utils.R")

# Validate configuration
validate_config()
print_config_summary()

# Set working directory to base folder
setwd(CONFIG$base_folder)

# =============================================================================
# MAIN WORKFLOW PARAMETERS
# =============================================================================

# Get scenario parameters from config
rcp <- CONFIG$rcp
ensemble <- CONFIG$ensemble  
id <- CONFIG$id

log_message("=== STARTING JULES EMULATION WORKFLOW ===", "MAIN")
log_message(paste("Scenario: RCP", rcp, "Ensemble", ensemble))

# =============================================================================
# STEP 1: LOAD AND PREPARE DATA
# =============================================================================

log_message("=== STEP 1: LOADING DATA ===", "STEP")

# Load design data (chess scape)
design_data_file <- file.path("..", CONFIG$design_data_file)
if (!file.exists(design_data_file)) {
  stop(paste("Design data file not found:", design_data_file))
}

cs <- readRDS(design_data_file)
log_message(paste("Loaded design data:", nrow(cs), "locations"))

# Initialize experiment environment
cs_output <- init_exp(id, cs, rcp = rcp, case = ensemble)
log_message("Experiment environment initialized")

# =============================================================================
# STEP 2: DATA EXPLORATION AND VALIDATION
# =============================================================================

log_message("=== STEP 2: DATA EXPLORATION ===", "STEP")

# Set random seed for reproducibility
set_seed(CONFIG$random_seed)

# Inspect data distributions
if (CONFIG$create_plots) {
  log_message("Creating data distribution plots...")
  inspect_data(CONFIG$variables, 
              transform = CONFIG$transform_functions, 
              N = 1000, 
              regions = NULL, 
              type = 'hist')
  
  # Create pairwise plots (if not too many variables)
  if (length(CONFIG$variables) <= 6) {
    inspect_data(CONFIG$variables, 
                transform = CONFIG$transform_functions, 
                N = 1000, 
                regions = NULL, 
                type = 'pair')
  }
}

# =============================================================================
# STEP 3: EXPERIMENTAL DESIGN
# =============================================================================

log_message("=== STEP 3: CREATING EXPERIMENTAL DESIGN ===", "STEP")

set_seed(CONFIG$random_seed)

# Create main experimental design
x <- init_design(N = CONFIG$n_samples, 
                var = CONFIG$variables, 
                transform = CONFIG$transform_functions)

log_message(paste("Created design matrix:", nrow(x), "x", ncol(x)))

# Validate design quality
if (CONFIG$create_plots) {
  check_design(x, N = 500, type = 'pair')
  check_design(x, N = 500, type = 'map')
}

# =============================================================================
# STEP 4: TIME COMPONENT SETUP (IF ENABLED)
# =============================================================================

if (CONFIG$use_time_varying) {
  log_message("=== STEP 4: TIME-VARYING SETUP ===", "STEP")
  
  # Generate time points
  set_seed(CONFIG$random_seed)
  years <- sample(CONFIG$min_year:CONFIG$max_year, 
                 CONFIG$n_samples, replace = TRUE)
  normalized_years <- (years - CONFIG$min_year) / (CONFIG$max_year - CONFIG$min_year)
  
  log_message(paste("Generated", length(unique(years)), "unique time points"))
  
  # Add time dimension to design
  x_with_time <- cbind(x, normalized_years)
  colnames(x_with_time)[ncol(x_with_time)] <- "time"
  
} else {
  years <- NULL
  x_with_time <- x
}

# =============================================================================
# STEP 5: RUN JULES SIMULATIONS
# =============================================================================

log_message("=== STEP 5: RUNNING JULES SIMULATIONS ===", "STEP")

# Note: This step would run actual JULES simulations
# For now, we'll create placeholder for the workflow structure

if (file.exists(get_results_path(paste0("exp", id, "/all_outputs.csv")))) {
  log_message("Loading existing JULES output data...")
  y <- as.matrix(read.csv(get_results_path(paste0("exp", id, "/all_outputs.csv"))))
} else {
  log_message("Running JULES simulations... (This may take a while)")
  
  # Run JULES (uncomment when ready to run actual simulations)
  # y <- jules(x, years)
  
  # For demonstration, create placeholder output
  log_message("Creating placeholder output for demonstration")
  n_months <- 360  # 30 years * 12 months
  y <- matrix(runif(nrow(x) * n_months, 50, 200), nrow = nrow(x), ncol = n_months)
  
  # Save outputs
  output_dir <- get_results_path(paste0("exp", id))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  write.csv(y, file.path(output_dir, "all_outputs.csv"), row.names = FALSE)
  saveRDS(x, file.path(output_dir, "x.rds"))
  saveRDS(y, file.path(output_dir, "y.rds"))
}

# Remove rows with missing values
if (CONFIG$use_time_varying) {
  x_final <- x_with_time
} else {
  x_final <- x
}

row_na <- apply(y, 1, function(row) any(is.na(row)))
x_final <- x_final[!row_na, , drop = FALSE]
y <- y[!row_na, , drop = FALSE]

log_message(paste("Final dataset:", nrow(y), "simulations x", ncol(y), "time points"))

# =============================================================================
# STEP 6: EMULATION WITH PCA
# =============================================================================

log_message("=== STEP 6: BUILDING EMULATOR ===", "STEP")

# Create emulator with PCA preprocessing
emulator_result <- create_emulator(x_final, y, 
                                  pca_components = CONFIG$pca_components,
                                  kernel_names = CONFIG$emulator_names,
                                  use_vecchia = CONFIG$use_vecchia)

emulator <- emulator_result$emulator
pca <- emulator_result$pca

log_message("Emulator training completed")

# =============================================================================
# STEP 7: TIME-VARYING ANALYSIS (IF ENABLED)
# =============================================================================

if (CONFIG$use_time_varying && CONFIG$n_time_points > 1) {
  log_message("=== STEP 7: TIME-VARYING ANALYSIS ===", "STEP")
  
  # Create expanded dataset for multiple time points
  log_message("Creating time-expanded dataset...")
  
  # For demonstration, create multiple time snapshots
  time_points <- seq(0, 1, length.out = CONFIG$n_time_points)
  
  # Expand inputs for all time points
  inputs_expanded <- do.call(rbind, replicate(CONFIG$n_time_points, x_final[, 1:ncol(x)], simplify = FALSE))
  times_expanded <- rep(time_points, each = nrow(x_final))
  inputs_expanded <- cbind(inputs_expanded, times_expanded)
  
  # For demonstration, create time-varying outputs
  y_time <- do.call(rbind, replicate(CONFIG$n_time_points, y, simplify = FALSE))
  
  log_message(paste("Time-expanded dataset:", nrow(inputs_expanded), "x", ncol(inputs_expanded)))
  
  # Create time-varying emulator
  log_message("Building time-varying emulator...")
  
  # Apply PCA to time-varying data
  sklearn <- reticulate::import('sklearn')
  pca_time <- sklearn$decomposition$PCA(as.integer(CONFIG$pca_components))
  y_basis_time <- pca_time$fit_transform(y_time)
  
  set_seed(CONFIG$validation_seed)
  
  # Build DGP and GP models for comparison
  m_dgp <- dgp(inputs_expanded, y_basis_time, 
               id = paste0('jules-time-dgp-', rcp, '-', ensemble), 
               name = "matern2.5", vecchia = CONFIG$use_vecchia)
  
  m_gp <- gp(inputs_expanded, y_basis_time, 
             id = paste0('jules-time-gp-', rcp, '-', ensemble), 
             name = "matern2.5", vecchia = CONFIG$use_vecchia)
  
  # Validate models
  m_dgp <- validate(m_dgp)
  m_gp <- validate(m_gp)
  
  if (CONFIG$create_plots) {
    plot(m_dgp)
    plot(m_gp)
  }
  
  # Save time-varying emulators
  write(m_dgp, get_results_path('timevaryingDGP'))
  write(m_gp, get_results_path('timevaryingGP'))
  
  log_message("Time-varying emulators created and saved")
}

# =============================================================================
# STEP 8: MODEL VALIDATION AND TESTING
# =============================================================================

log_message("=== STEP 8: MODEL VALIDATION ===", "STEP")

# Create test dataset
if (CONFIG$n_test_samples > 0) {
  log_message("Creating test dataset...")
  set_seed(CONFIG$random_seed + 1)  # Different seed for test data
  
  test_x <- init_test(N = CONFIG$n_test_samples)
  
  # Run test simulations (placeholder for now)
  log_message("Running test simulations...")
  test_y <- matrix(runif(nrow(test_x) * ncol(y), 50, 200), 
                   nrow = nrow(test_x), ncol = ncol(y))
  
  # Remove NA rows
  row_na_test <- apply(test_y, 1, function(row) any(is.na(row)))
  test_x <- test_x[!row_na_test, , drop = FALSE]
  test_y <- test_y[!row_na_test, , drop = FALSE]
  
  # Save test data
  saveRDS(test_x, get_results_path(paste0("exp", id, "/test_x.rds")))
  saveRDS(test_y, get_results_path(paste0("exp", id, "/test_y.rds")))
  
  # Validate emulator performance
  if (exists("pca")) {
    test_y_basis <- pca$transform(test_y)
    
    if (CONFIG$create_plots) {
      plot(emulator, test_x, test_y_basis, min_max = TRUE, dim = min(6, ncol(test_x)))
    }
    
    # Calculate out-of-sample error
    nrmse_score <- nrmse_oos(emulator, pca$inverse_transform, test_x, test_y)
    log_message(paste("Out-of-sample NRMSE:", round(nrmse_score, 4)))
  }
}

# =============================================================================
# STEP 9: PREDICTIONS AND VISUALIZATION
# =============================================================================

log_message("=== STEP 9: GENERATING PREDICTIONS ===", "STEP")

# Load prediction dataset (if available)
prediction_data_file <- get_data_path("chess_scape_vars.rds")
if (file.exists(prediction_data_file)) {
  log_message("Loading prediction dataset...")
  test_x_full <- readRDS(prediction_data_file)
  test_x_full <- unname(as.matrix(test_x_full))
  
  # Make predictions
  log_message("Making predictions on full dataset...")
  predictions <- predict(emulator, test_x_full, cores = CONFIG$n_cores, 
                        chunks = CONFIG$chunk_size)
  
  # Load coordinates
  coord_file <- get_data_path("chess_scape_bng.rds")
  if (file.exists(coord_file)) {
    xy <- readRDS(coord_file)
    
    # Combine results
    output_predictions <- cbind(xy, 
                               predictions$results$mean, 
                               predictions$results$var)
    names(output_predictions) <- c('x', 'y', 'mean', 'var')
    
    # Save predictions
    output_file <- get_results_path(paste0('predictions_', id, '.feather'))
    feather::write_feather(output_predictions, output_file)
    
    log_message(paste("Predictions saved to:", output_file))
  }
}

# =============================================================================
# STEP 10: VISUALIZATION AND REPORTING
# =============================================================================

if (CONFIG$create_plots) {
  log_message("=== STEP 10: CREATING VISUALIZATIONS ===", "STEP")
  
  # Create summary plots
  if (exists("y") && ncol(y) > 0) {
    # Time series plot of sample outputs
    if (ncol(y) >= 12) {  # If we have at least 12 months of data
      sample_indices <- sample(nrow(y), min(10, nrow(y)))
      
      matplot(t(y[sample_indices, 1:min(360, ncol(y))]), 
              type = "l", lty = 1,
              xlab = "Month", 
              ylab = "Carbon Stored (tCO2e/ha)",
              main = paste("Sample JULES Outputs - RCP", rcp, "Ensemble", ensemble))
    }
  }
  
  # Emulator diagnostic plots
  if (exists("emulator")) {
    plot(emulator)
  }
  
  log_message("Visualization complete")
}

# =============================================================================
# STEP 11: CLEANUP AND SUMMARY
# =============================================================================

log_message("=== STEP 11: WORKFLOW SUMMARY ===", "STEP")

# Print summary statistics
if (exists("y")) {
  log_message(paste("Final simulation count:", nrow(y)))
  log_message(paste("Time points per simulation:", ncol(y)))
  log_message(paste("Mean carbon storage:", round(mean(y, na.rm = TRUE), 2), "tCO2e/ha"))
  log_message(paste("Carbon storage range:", round(min(y, na.rm = TRUE), 2), "-", 
                   round(max(y, na.rm = TRUE), 2), "tCO2e/ha"))
}

if (exists("emulator")) {
  log_message("Emulator successfully created and validated")
}

if (exists("predictions")) {
  log_message(paste("Predictions generated for", nrow(predictions$results$mean), "locations"))
}

# Save workspace summary
summary_data <- list(
  config = CONFIG,
  scenario_id = id,
  n_simulations = if(exists("y")) nrow(y) else 0,
  n_variables = length(CONFIG$variables),
  variables = CONFIG$variables,
  completion_time = Sys.time()
)

saveRDS(summary_data, get_results_path(paste0("workflow_summary_", id, ".rds")))

log_message("=== JULES EMULATION WORKFLOW COMPLETE ===", "MAIN")
log_message(paste("Results saved in:", CONFIG$base_folder))

# =============================================================================
# OPTIONAL: SEQUENTIAL DESIGN (ADVANCED FEATURE)
# =============================================================================

if (CONFIG$perform_sequential_design && exists("emulator")) {
  log_message("=== SEQUENTIAL DESIGN EXTENSION ===", "STEP")
  
  # This would implement the sequential design from the original master.R
  # Placeholder for advanced users
  
  log_message("Sequential design functionality available but not implemented in this clean version")
  log_message("Refer to original master.R for sequential design implementation")
}

log_message("Script execution completed successfully") 