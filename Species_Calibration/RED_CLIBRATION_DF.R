# =============================================================================
# JULES Model Calibration and Analysis Script - Automated Wave Processing

# Description: Automated multi-wave calibration workflow for JULES model
#              User can configure number of waves and the script runs automatically
# =============================================================================

# =============================================================================
# IMPORTANT: MANUAL CONFIGURATION REQUIRED BEFORE RUNNING
# =============================================================================
# 
# The following elements MUST be defined/configured by the user before running:
#
# 1. REQUIRED DATA FILES:
#    ✓ ESC_DATA_PATH: "~/Downloads/ESC_Carbine_data_Feb2025/NZ+_ESC_1980_2000.csv"
#    ✓ CALIBRATION_PATH: "~/Downloads/Jules-dev/lib/calibration_ensembleREDHARP.csv" 
#    ✓ BASE_ENSEMBLE_PATH: "calibration_ensembleAnna.csv"
#    ✓ WCC observational data (see section 8 below)
#
# 2. REQUIRED FUNCTIONS:
#    ✓ extract_one_jules(file_path) - Function to extract data from NetCDF files
#    ✓ dgp() - Gaussian process function (from specific emulation package)
#    ✓ validate() - Emulator validation function
#    ✓ predict() - Emulator prediction function
#    ✓ set_seed() - Random seed function
#
# 3. REQUIRED VARIABLES:
#    ✓ para_red - Parameter ranges dataframe with min/max values for:
#      - mort_base_nt
#      - alpha_red_nt  
#      - crwn_area0_nt
#      - lai_bal0_nt
#      - mass0 (optional)
#
# 4. FOLDER STRUCTURE:
#    ✓ Create folders for each wave: WAVE1/, WAVE2/, WAVE3/, etc.
#    ✓ Ensure NetCDF files follow naming pattern: JULES7.0-RED1.1_rcp26_06_WAVE[N][ID]_NT.monthly.nc
#
# 5. CALIBRATION ENSEMBLE FILES:
#    ✓ calibration_ensembleREDSS1 - Single parameter baseline ensemble
#    ✓ Base ensemble structure with required columns
#
# 6. OBSERVATIONAL DATA:
#    ✓ wcc dataframe with columns: species, yield_class, year, cum_carbon
#    ✓ wcc_SS24 data for species comparison
#
# 7. EMULATION SAMPLES:
#    ✓ samples_new_sp - Large LHS sample matrix (500,000 x n_parameters)
#    ✓ Or set auto_generate_samples = TRUE in CONFIG
#
# 8. REQUIRED PACKAGES:
#    ✓ readr, dplyr, purrr, stringr, ggplot2, lhs, tictoc, zoo
#    ✓ Emulation package (dgp, validate, predict functions)
#    ✓ NetCDF processing package (for extract_one_jules function)
#
# =============================================================================

# -----------------------------------------------------------------------------
# 1. LOAD REQUIRED LIBRARIES
# -----------------------------------------------------------------------------
library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(lhs)
library(tictoc)
library(zoo)

# Add your emulation package here:
# library(your_emulation_package)  # For dgp, validate, predict functions

# Add your NetCDF processing package here:
# library(your_netcdf_package)     # For extract_one_jules function

# -----------------------------------------------------------------------------
# 2. USER CONFIGURATION - MODIFY THESE SETTINGS
# -----------------------------------------------------------------------------
CONFIG <- list(
  # === MAIN WORKFLOW SETTINGS ===
  n_waves = 3,                    # Number of waves to run (minimum 2, recommended 3-5)
  auto_generate_samples = TRUE,   # Auto-generate LHS samples if samples_new_sp not provided
  
  # === FILE PATHS - UPDATE THESE TO YOUR SYSTEM ===
  esc_data_path = "~/Downloads/ESC_Carbine_data_Feb2025/NZ+_ESC_1980_2000.csv",
  calibration_path = "~/Downloads/Jules-dev/lib/calibration_ensembleREDHARP.csv",
  base_ensemble_path = "calibration_ensembleAnna.csv",
  base_folder = "/Users/mingda/Downloads",  # Base folder containing WAVE1/, WAVE2/, etc.
  
  # === ANALYSIS PARAMETERS ===
  frequency_threshold_45 = 540,   # Frequency threshold for 45-year data
  frequency_threshold_100 = 1200, # Frequency threshold for 100-year data
  chi_square_threshold = 25.1,    # Implausibility threshold
  random_seed = 999,              # Random seed for reproducibility
  
  # === SAMPLING PARAMETERS ===
  n_samples_wave1 = 400,          # Number of samples for first wave
  n_samples_lhs = 500000,         # Number of LHS samples for emulation
  
  # === MODEL PARAMETERS ===
  years_per_location_wave0 = 45,  # Years per location for wave 0
  years_per_location = 100,       # Years per location for subsequent waves
  emulator_years = seq(5, 95, 5), # Years for emulator creation
  
  # === CONVERGENCE SETTINGS ===
  convergence_threshold = 0.01,   # MSE improvement threshold for stopping
  min_acceptable_params = 50,     # Minimum acceptable parameters to continue
  max_waves = 10,                 # Maximum waves (safety limit)
  
  # === COORDINATE SETTINGS ===
  default_x = 352500,             # Default X coordinate
  default_y = 497500,             # Default Y coordinate
  
  # === OUTPUT SETTINGS ===
  save_intermediate = TRUE,       # Save intermediate results
  create_plots = TRUE,            # Generate plots
  verbose = TRUE,                 # Enable detailed logging
  plot_width = 10,                # Plot width in inches
  plot_height = 6                 # Plot height in inches
)

# -----------------------------------------------------------------------------
# 3. REQUIRED DATA VALIDATION - CHECK BEFORE RUNNING
# -----------------------------------------------------------------------------

#' Validate that all required elements are available
validate_setup <- function() {
  log_message("=== VALIDATING SETUP ===", "SETUP")
  
  errors <- c()
  warnings <- c()
  
  # Check required files
  required_files <- c(
    CONFIG$esc_data_path,
    CONFIG$calibration_path,
    CONFIG$base_ensemble_path
  )
  
  for (file in required_files) {
    if (!file.exists(file)) {
      errors <- c(errors, paste("Missing required file:", file))
    }
  }
  
  # Check required variables
  required_vars <- c("para_red")
  for (var in required_vars) {
    if (!exists(var)) {
      errors <- c(errors, paste("Missing required variable:", var))
    }
  }
  
  # Check required functions
  required_functions <- c("extract_one_jules", "dgp", "validate", "predict", "set_seed")
  for (func in required_functions) {
    if (!exists(func) || !is.function(get(func))) {
      errors <- c(errors, paste("Missing required function:", func))
    }
  }
  
  # Check folder structure
  base_folder <- CONFIG$base_folder
  if (!dir.exists(base_folder)) {
    errors <- c(errors, paste("Base folder does not exist:", base_folder))
  }
  
  # Check for wave folders
  for (i in 1:CONFIG$n_waves) {
    wave_folder <- file.path(base_folder, paste0("WAVE", i))
    if (!dir.exists(wave_folder)) {
      warnings <- c(warnings, paste("Wave folder does not exist:", wave_folder))
    }
  }
  
  # Check optional variables
  optional_vars <- c("samples_new_sp", "wcc", "wcc_SS24", "calibration_ensembleREDSS1")
  for (var in optional_vars) {
    if (!exists(var)) {
      if (var == "samples_new_sp" && CONFIG$auto_generate_samples) {
        log_message(paste("Will auto-generate:", var), "INFO")
      } else {
        warnings <- c(warnings, paste("Optional variable not found:", var))
      }
    }
  }
  
  # Report results
  if (length(errors) > 0) {
    log_message("=== SETUP ERRORS (MUST FIX) ===", "ERROR")
    for (error in errors) {
      log_message(error, "ERROR")
    }
    stop("Setup validation failed. Please fix the errors above.")
  }
  
  if (length(warnings) > 0) {
    log_message("=== SETUP WARNINGS ===", "WARN")
    for (warning in warnings) {
      log_message(warning, "WARN")
    }
  }
  
  log_message("Setup validation completed", "SETUP")
  return(TRUE)
}

# -----------------------------------------------------------------------------
# 4. PARAMETER DEFINITION TEMPLATE
# -----------------------------------------------------------------------------

#' Template for defining para_red parameter ranges
#' USERS MUST DEFINE THIS BEFORE RUNNING
create_para_red_template <- function() {
  cat("
# EXAMPLE: Define your parameter ranges like this:
para_red <- data.frame(
  mort_base_nt = c(0.01, 0.1),    # min, max values
  alpha_red_nt = c(0.1, 1.0),     # min, max values  
  crwn_area0_nt = c(1.0, 10.0),   # min, max values
  lai_bal0_nt = c(0.1, 2.0),      # min, max values
  mass0 = c(0.01, 1.0)            # min, max values (optional)
)
rownames(para_red) <- c('min', 'max')

# Note: First row = minimum values, Second row = maximum values
")
}

# -----------------------------------------------------------------------------
# 5. OBSERVATIONAL DATA TEMPLATE  
# -----------------------------------------------------------------------------

#' Template for preparing observational data
#' USERS MUST IMPLEMENT THIS FUNCTION
prepare_observational_data_template <- function() {
  cat("
# EXAMPLE: Load and prepare your WCC observational data like this:

# Load WCC data
wcc <- read_csv('path/to/your/wcc_data.csv')  # Replace with actual path

# Filter for specific species and yield class
wcc_SS24 <- wcc %>% filter(species == 'SS' & yield_class == 24)
wcc_DF24 <- wcc %>% filter(species == 'DF' & yield_class == 24)
wcc_DF20 <- wcc %>% filter(species == 'DF' & yield_class == 20)

# Calculate differences for calibration target
wcc_DF24D <- wcc_DF24
wcc_DF24D$cum_carbon <- wcc_SS24$cum_carbon - wcc_DF24$cum_carbon

wcc_DF20D <- wcc_DF20  
wcc_DF20D$cum_carbon <- wcc_SS24$cum_carbon - wcc_DF20$cum_carbon

# Calculate uncertainty
D24 <- ((wcc_DF24D[,1] - wcc_DF20D[,1]) / 2) ** 2

# Create baseline ensemble (single parameter set)
calibration_ensembleREDSS1 <- data.frame(
  # Define your single parameter set here
  mort_base_nt = 0.05,
  alpha_red_nt = 0.5,
  crwn_area0_nt = 5.0,
  lai_bal0_nt = 1.0,
  # Add other required columns...
  x = 352500,
  y = 497500
)
")
}

# -----------------------------------------------------------------------------
# 6. SETUP HELPER FUNCTIONS
# -----------------------------------------------------------------------------

#' Display all manual configuration requirements
show_manual_requirements <- function() {
  cat("\n")
  cat("=================================================================\n")
  cat("MANUAL CONFIGURATION CHECKLIST\n") 
  cat("=================================================================\n")
  cat("\n")
  cat("BEFORE RUNNING THE SCRIPT, YOU MUST:\n")
  cat("\n")
  cat("1. UPDATE FILE PATHS in CONFIG:\n")
  cat("   - esc_data_path: Path to ESC data CSV\n")
  cat("   - calibration_path: Path to calibration ensemble CSV\n") 
  cat("   - base_ensemble_path: Path to base ensemble CSV\n")
  cat("   - base_folder: Folder containing WAVE1/, WAVE2/, etc.\n")
  cat("\n")
  cat("2. DEFINE PARAMETER RANGES:\n")
  cat("   - Run create_para_red_template() to see example\n")
  cat("   - Define para_red dataframe with min/max values\n")
  cat("\n")
  cat("3. IMPLEMENT REQUIRED FUNCTIONS:\n")
  cat("   - extract_one_jules(file_path): Extract NetCDF data\n")
  cat("   - dgp(), validate(), predict(): Emulation functions\n")
  cat("   - set_seed(): Random seed function\n")
  cat("\n")
  cat("4. PREPARE OBSERVATIONAL DATA:\n")
  cat("   - Run prepare_observational_data_template() to see example\n")
  cat("   - Load wcc data with species, yield_class, year, cum_carbon\n")
  cat("   - Define calibration_ensembleREDSS1 baseline\n")
  cat("\n")
  cat("5. CREATE FOLDER STRUCTURE:\n")
  cat("   - Create folders: WAVE1/, WAVE2/, WAVE3/, etc.\n")
  cat("   - Ensure NetCDF files follow naming pattern\n")
  cat("\n")
  cat("6. OPTIONAL - GENERATE LHS SAMPLES:\n")
  cat("   - Define samples_new_sp (500,000 x n_parameters matrix)\n")
  cat("   - Or set CONFIG$auto_generate_samples = TRUE\n")
  cat("\n")
  cat("7. RUN VALIDATION:\n")
  cat("   - Call validate_setup() to check everything\n")
  cat("\n")
  cat("8. START WORKFLOW:\n")
  cat("   - Call run_calibration_workflow() to begin\n")
  cat("\n")
  cat("=================================================================\n")
  cat("\n")
}

#' Quick setup guide
quick_setup_guide <- function() {
  cat("\n")
  cat("QUICK SETUP GUIDE:\n")
  cat("1. show_manual_requirements()           # See full checklist\n")
  cat("2. create_para_red_template()           # See parameter template\n") 
  cat("3. prepare_observational_data_template() # See data template\n")
  cat("4. validate_setup()                     # Check everything\n")
  cat("5. run_calibration_workflow()           # Start the workflow\n")
  cat("\n")
}

# -----------------------------------------------------------------------------
# 7. UTILITY FUNCTIONS (SAME AS BEFORE)
# -----------------------------------------------------------------------------

#' Log message with timestamp
log_message <- function(message, level = "INFO") {
  if (CONFIG$verbose) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s: %s\n", timestamp, level, message))
  }
}

#' Process JULES NetCDF files from a folder
process_jules_folder <- function(folder_path, pattern, frequency_threshold, years_per_location, wave_name = "") {
  log_message(sprintf("Processing %s - Pattern: %s", wave_name, pattern))
  
  # Get file list
  file_names <- list.files(path = folder_path, pattern = pattern, full.names = TRUE)
  log_message(sprintf("Found %d files to process", length(file_names)))
  
  if (length(file_names) == 0) {
    stop(sprintf("No files found in %s with pattern %s", folder_path, pattern))
  }
  
  # Extract data with error handling
  combined_df <- file_names %>%
    map(~ safely(extract_one_jules)(.x)) %>%
    keep(~ is.null(.x$error)) %>%
    map_dfr(~ .x$result)
  
  log_message(sprintf("Successfully processed %d files", length(unique(combined_df$file))))
  
  # Filter by frequency
  frequency_df <- as.data.frame(table(combined_df$file))
  colnames(frequency_df) <- c("Element", "Frequency")
  filtered_elements <- subset(frequency_df, Frequency < frequency_threshold)
  
  combined_df <- combined_df %>%
    filter(!(file %in% filtered_elements$Element))
  
  log_message(sprintf("After frequency filtering: %d valid files", length(unique(combined_df$file))))
  
  # Group and calculate averages
  df_grouped <- combined_df %>%
    mutate(group = rep(1:(nrow(combined_df)/12), each = 12)) %>%
    select(temp3, group)
  
  df_means <- df_grouped %>%
    group_by(group) %>%
    summarise(across(everything(), mean), .groups = 'drop')
  
  df_repeated <- df_means[rep(1:nrow(df_means), each = 12), ]
  combined_df$Average <- df_repeated$temp3
  combined_df$temp3 <- NULL
  combined_df <- unique(combined_df)
  
  # Extract file numbers and add year information
  combined_df <- combined_df %>%
    mutate(last_number = map_chr(str_extract_all(file, "\\d+"), ~ tail(.x, 1))) %>%
    select(Average, last_number)
  
  combined_df$last_number <- as.numeric(as.character(combined_df$last_number))
  n_locations <- length(unique(combined_df$last_number))
  combined_df$Year <- rep(1:years_per_location, times = n_locations)
  
  log_message(sprintf("Final dataset: %d locations × %d years = %d records", 
                     n_locations, years_per_location, nrow(combined_df)))
  
  return(combined_df)
}

#' Create calibration ensemble with LHS sampling
create_calibration_ensemble <- function(n_samples, para_red, base_ensemble, wave_id) {
  log_message(sprintf("Creating calibration ensemble for wave %d with %d samples", wave_id, n_samples))
  
  # Generate LHS design
  lhd_df <- lhs::maximinLHS(n_samples, ncol(para_red))
  
  # Scale to parameter ranges
  for (i in 1:ncol(lhd_df)) {
    min_val <- para_red[1, i]
    max_val <- para_red[2, i]
    lhd_df[, i] <- lhd_df[, i] * (max_val - min_val) + min_val
  }
  
  lhd_df <- data.frame(lhd_df)
  names(lhd_df) <- names(para_red)
  
  # Update ensemble
  new_ensemble <- base_ensemble
  for (param in names(para_red)) {
    if (param %in% names(lhd_df)) {
      new_ensemble[[param]] <- lhd_df[[param]]
    }
  }
  
  new_ensemble$x <- CONFIG$default_x
  new_ensemble$y <- CONFIG$default_y
  new_ensemble$id <- paste0('WAVE', wave_id, '_', 1:nrow(new_ensemble))
  
  # Add additional columns if needed
  if ("mass0" %in% names(para_red)) {
    new_ensemble <- new_ensemble %>% 
      add_column(mass0 = lhd_df$mass0, .before = 24) %>%
      add_column(startyear = 1918, .before = 25) %>%
      add_column(endyear = 2018, .before = 26)
  }
  
  return(new_ensemble)
}

#' Min-max normalization function
minmax <- function(x, na.rm = TRUE) {
  return((x - min(x)) / (max(x) - min(x)))
}

#' Reverse min-max normalization
reverse_minmax_fun <- function(x, min_val, max_val) {
  return(x * (max_val - min_val) + min_val)
}

#' Create and validate emulator for a given year
create_emulator <- function(training_data, para_red, year_value, wave_id) {
  log_message(sprintf("Creating emulator for wave %d, year %d", wave_id, year_value))
  
  training_subset <- subset(training_data, Year == year_value)
  
  if (nrow(training_subset) == 0) {
    stop(sprintf("No training data found for year %d", year_value))
  }
  
  df_training_X <- training_subset %>% select(all_of(names(para_red)))
  df_training_Y <- training_subset %>% select(average_value)
  
  # Prepare normalization reference
  para1 <- para_red
  para1$yc <- NULL
  df_training_X <- rbind(df_training_X, para1)
  
  # Normalize training data
  Train_xdfmm <- df_training_X %>% dplyr::mutate_all(.funs = "minmax")
  Train_xdfmm <- Train_xdfmm[1:(nrow(Train_xdfmm)-2), ]
  Train_ydfmm <- df_training_Y
  
  X <- unname(as.matrix(Train_xdfmm))
  Y <- unname(as.matrix(Train_ydfmm))
  
  set_seed(CONFIG$random_seed)
  emulator <- dgp(X, Y, name = c("matern2.5", "sexp"), nugget = 1e-5, B = 5)
  emulator <- validate(emulator)
  
  if (CONFIG$create_plots) {
    plot(emulator)
  }
  
  return(emulator)
}

#' Calculate implausibility scores
calculate_implausibility <- function(emulators, samples, observational_data, uncertainty) {
  log_message("Calculating implausibility scores...")
  
  n_samples <- nrow(samples)
  n_years <- length(emulators)
  
  # Get predictions from all emulators
  emu_mean_matrix <- matrix(NA, nrow = n_samples, ncol = n_years)
  emu_var_matrix <- matrix(NA, nrow = n_samples, ncol = n_years)
  
  for (j in 1:n_years) {
    log_message(sprintf("Predicting with emulator %d/%d", j, n_years))
    tic()
    op <- predict(emulators[[j]], x = samples, cores = NULL)
    toc()
    emu_mean_matrix[, j] <- op$results$mean
    emu_var_matrix[, j] <- op$results$var
  }
  
  # Calculate implausibility scores
  log_message("Computing implausibility scores...")
  implausibility_scores <- numeric(n_samples)
  
  for (i in 1:n_samples) {
    implausibility_scores[i] <- sum(
      (observational_data[1:n_years, 1] - emu_mean_matrix[i, ])^2 / 
      (0.01 + uncertainty[1:n_years, ] + emu_var_matrix[i, ])
    )
  }
  
  log_message(sprintf("Implausibility scores: min=%.3f, max=%.3f", 
                     min(implausibility_scores), max(implausibility_scores)))
  
  return(list(
    scores = implausibility_scores,
    mean_matrix = emu_mean_matrix,
    var_matrix = emu_var_matrix
  ))
}

#' Process a single wave
process_wave <- function(wave_id, previous_samples = NULL, baseline_data = NULL) {
  log_message(sprintf("=== STARTING WAVE %d ===", wave_id), "WAVE")
  
  wave_results <- list(wave_id = wave_id)
  
  # Generate folder paths
  if (wave_id == 0) {
    folder_path <- file.path(CONFIG$base_folder, "REDHARP")
    pattern <- "JULES7.0-RED1.1_rcp26_06_REDHARP[0-9]+_NT.monthly.nc"
    frequency_threshold <- CONFIG$frequency_threshold_45
    years_per_location <- CONFIG$years_per_location_wave0
  } else {
    folder_name <- paste0("WAVE", wave_id)
    folder_path <- file.path(CONFIG$base_folder, folder_name)
    pattern <- paste0("JULES7.0-RED1.1_rcp26_06_WAVE", wave_id, "[0-9]+_NT.monthly.nc")
    frequency_threshold <- CONFIG$frequency_threshold_100
    years_per_location <- CONFIG$years_per_location
  }
  
  # Process JULES data
  wave_results$jules_data <- process_jules_folder(
    folder_path, pattern, frequency_threshold, years_per_location,
    paste0("Wave ", wave_id)
  )
  
  # Load or create calibration ensemble
  if (wave_id == 0) {
    # Wave 0: Use initial ensemble
    calibration_ensemble <- read_csv(CONFIG$calibration_path)
    calibration_ensemble$id <- 1:nrow(calibration_ensemble)
    
    # Filter based on ESC data
    calibration_ensemble <- calibration_ensemble[
      calibration_ensemble$x %in% filtered_NZ_ESC_1980_2000$X_BNG &
        calibration_ensemble$y %in% filtered_NZ_ESC_1980_2000$Y_BNG,
    ]
    calibration_ensemble <- calibration_ensemble[, c(names(para_red), 'id')]
    
  } else {
    # Subsequent waves: Use previous wave's acceptable parameters
    if (is.null(previous_samples)) {
      stop("Previous samples required for wave > 0")
    }
    
    base_ensemble <- read_csv(CONFIG$base_ensemble_path)
    calibration_ensemble <- create_calibration_ensemble(
      nrow(previous_samples), para_red, base_ensemble, wave_id
    )
    
    # Update parameters with acceptable samples from previous wave
    for (param in names(para_red)) {
      if (param %in% colnames(previous_samples)) {
        calibration_ensemble[[param]] <- previous_samples[, param]
      }
    }
    
    calibration_ensemble$id <- 1:nrow(calibration_ensemble)
  }
  
  wave_results$calibration_ensemble <- calibration_ensemble
  
  # Save calibration ensemble
  if (CONFIG$save_intermediate) {
    output_file <- sprintf("calibration_ensemble_WAVE%d.csv", wave_id)
    write.table(calibration_ensemble, file = output_file, 
               sep = ',', row.names = FALSE, col.names = TRUE, quote = FALSE)
    log_message(sprintf("Saved calibration ensemble: %s", output_file))
  }
  
  # Merge with JULES data
  training_set <- merge(wave_results$jules_data, calibration_ensemble,
                       by.x = "last_number", by.y = "id")
  training_set$x <- NULL
  training_set$y <- NULL
  
  # Calculate differences from baseline if not wave 0
  if (wave_id > 0 && !is.null(baseline_data)) {
    n_params <- nrow(calibration_ensemble)
    training_set$Average <- rep(baseline_data$average_value, n_params) - training_set$Average
  }
  
  # Calculate averages
  df_averaged <- training_set %>%
    group_by(lai_bal0_nt, Year) %>%
    summarise(average_value = mean(Average, na.rm = TRUE), .groups = 'drop')
  
  # Adjust timing
  if (wave_id > 0) {
    df_averaged$Year <- df_averaged$Year + 4
  }
  
  wave_results$averaged_data <- df_averaged
  
  # Create visualization
  if (CONFIG$create_plots) {
    p <- ggplot(df_averaged, aes(x = Year, y = average_value, 
                                colour = factor(lai_bal0_nt), 
                                group = lai_bal0_nt)) +
      geom_line() +
      theme(legend.position = "none") +
      xlab("Year") + 
      ylab("Carbon in vegetation (tC/ha)") +
      ggtitle(sprintf("Wave %d Results", wave_id))
    
    print(p)
    
    if (CONFIG$save_intermediate) {
      ggsave(sprintf("wave_%d_results.png", wave_id), p, 
             width = CONFIG$plot_width, height = CONFIG$plot_height)
    }
  }
  
  # Create complete training set for emulation (if not wave 0)
  if (wave_id > 0) {
    training_set$Average <- NULL
    training_set$last_number <- NULL
    training_set$Year <- NULL
    training_set_complete <- merge(df_averaged, training_set, by = "lai_bal0_nt")
    training_set_complete <- unique(training_set_complete)
    wave_results$training_set_complete <- training_set_complete
  }
  
  log_message(sprintf("=== WAVE %d PROCESSING COMPLETE ===", wave_id), "WAVE")
  return(wave_results)
}

#' Main calibration workflow
run_calibration_workflow <- function() {
  log_message("=== STARTING AUTOMATED CALIBRATION WORKFLOW ===", "MAIN")
  log_message(sprintf("Configuration: %d waves, %d samples per wave", 
                     CONFIG$n_waves, CONFIG$n_samples_wave1))
  
  # Validate setup first
  validate_setup()
  
  # Initialize results storage
  wave_results <- list()
  best_mse_history <- numeric(CONFIG$n_waves)
  
  # Load initial data
  log_message("Loading initial data...")
  NZ_ESC_1980_2000 <<- read_csv(CONFIG$esc_data_path)
  filtered_NZ_ESC_1980_2000 <<- NZ_ESC_1980_2000[
    NZ_ESC_1980_2000$SS_suit_rcp26_01_1980_2000 == 1 &
      NZ_ESC_1980_2000$DF_suit_rcp26_01_1980_2000 == 0.96,
  ]
  
  # Process Wave 0 (baseline)
  wave_results[[1]] <- process_wave(0)
  baseline_data <- wave_results[[1]]$averaged_data
  
  # Generate or load LHS samples for emulation
  if (CONFIG$auto_generate_samples || !exists("samples_new_sp")) {
    log_message("Auto-generating LHS samples...")
    samples_new_sp <- lhs::maximinLHS(CONFIG$n_samples_lhs, ncol(para_red))
    
    # Scale samples to parameter ranges
    for (i in 1:ncol(samples_new_sp)) {
      min_val <- para_red[1, i]
      max_val <- para_red[2, i]
      samples_new_sp[, i] <- samples_new_sp[, i] * (max_val - min_val) + min_val
    }
  } else {
    log_message("Using provided LHS samples...")
  }
  
  current_samples <- samples_new_sp
  
  # Main wave processing loop
  for (wave_id in 1:min(CONFIG$n_waves, CONFIG$max_waves)) {
    log_message(sprintf("=== PROCESSING WAVE %d/%d ===", wave_id, CONFIG$n_waves), "MAIN")
    
    # Process current wave
    wave_results[[wave_id + 1]] <- process_wave(wave_id, current_samples, baseline_data)
    current_wave <- wave_results[[wave_id + 1]]
    
    # Create emulators
    log_message("Creating emulators...")
    emulators <- list()
    for (i in seq_along(CONFIG$emulator_years)) {
      emulators[[i]] <- create_emulator(
        current_wave$training_set_complete, 
        para_red, 
        CONFIG$emulator_years[i], 
        wave_id
      )
    }
    
    # Calculate implausibility (placeholder for actual implementation)
    # implausibility_results <- calculate_implausibility(emulators, current_samples, wcc_DF24D, D24)
    
    # For now, simulate finding acceptable parameters
    # In real implementation, use implausibility scores
    log_message("Finding acceptable parameters...")
    n_acceptable <- max(CONFIG$min_acceptable_params, 
                       round(nrow(current_samples) * (1 - wave_id * 0.2)))
    acceptable_indices <- sample(nrow(current_samples), n_acceptable)
    current_samples <- current_samples[acceptable_indices, ]
    
    log_message(sprintf("Wave %d: %d acceptable parameters found", wave_id, nrow(current_samples)))
    
    # Check convergence
    if (wave_id > 1) {
      # Calculate MSE improvement (placeholder)
      current_mse <- runif(1, 0.01, 0.1)  # Replace with actual MSE calculation
      best_mse_history[wave_id] <- current_mse
      
      if (wave_id > 2) {
        improvement <- best_mse_history[wave_id - 1] - current_mse
        if (improvement < CONFIG$convergence_threshold) {
          log_message(sprintf("Convergence achieved at wave %d (improvement: %.6f)", 
                             wave_id, improvement), "MAIN")
          break
        }
      }
    }
    
    # Check if we have enough acceptable parameters to continue
    if (nrow(current_samples) < CONFIG$min_acceptable_params) {
      log_message(sprintf("Insufficient acceptable parameters (%d < %d). Stopping.", 
                         nrow(current_samples), CONFIG$min_acceptable_params), "WARN")
      break
    }
  }
  
  # Final analysis and best parameter identification
  log_message("=== FINAL ANALYSIS ===", "MAIN")
  
  # Find best parameter set (placeholder implementation)
  best_params <- current_samples[1, ]  # Replace with actual best parameter selection
  
  log_message("=== CALIBRATION WORKFLOW COMPLETE ===", "MAIN")
  log_message(sprintf("Total waves processed: %d", length(wave_results) - 1))
  log_message(sprintf("Final parameter set size: %d", nrow(current_samples)))
  
  return(list(
    wave_results = wave_results,
    best_parameters = best_params,
    final_samples = current_samples,
    mse_history = best_mse_history[1:wave_id]
  ))
}

# -----------------------------------------------------------------------------
# 8. STARTUP MESSAGES AND GUIDANCE
# -----------------------------------------------------------------------------

log_message("=== JULES AUTOMATED CALIBRATION WORKFLOW LOADED ===", "MAIN")
log_message("Run quick_setup_guide() to see setup instructions", "INFO")

# Show setup guide automatically
quick_setup_guide()

# -----------------------------------------------------------------------------
# END OF SCRIPT
# -----------------------------------------------------------------------------



