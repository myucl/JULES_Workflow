# =============================================================================
# JULES Utilities - Essential Functions Only
# Description: Core utility functions for JULES model emulation
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(readr)
  library(tidync)
  library(tidyr)
  library(tidyverse)
  library(dgpsi)
  library(reticulate)
  library(clhs)
})

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

# Create experiment environment to store global variables
exp.env <- new.env()

#' Initialize experiment environment
#' @param id Experiment ID (e.g., "85_15")
#' @param cs Chess scape data
#' @param rcp RCP scenario
#' @param case Ensemble case
init_exp <- function(id, cs, rcp = "85", case = "15") {
  log_message("Initializing experiment environment...")
  
  # Store in experiment environment
  assign('exp_id', id, exp.env)
  assign('rcp', rcp, exp.env)
  assign('case', case, exp.env)
  assign('chess_scape', cs, exp.env)
  
  # Create lookup table for efficiency
  cs_lookup <- cs %>% 
    select(all_of(CONFIG$variables)) %>%
    apply(1, paste, collapse = ", ")
  assign('chess_scape_lookup_table', cs_lookup, exp.env)
  
  # Set paths (these would need to be configured based on your system)
  assign('pd_path', "lib/add_ssh_key.exp.sh", exp.env)
  assign('sh_path', "lib/ai4nz-vn2.sh", exp.env)
  
  log_message(paste("Experiment initialized:", id))
  return(list(
    exp_id = id,
    n_locations = nrow(cs),
    n_variables = length(CONFIG$variables)
  ))
}

# =============================================================================
# DESIGN FUNCTIONS
# =============================================================================

#' Initialize experimental design using LHS sampling
#' @param N Number of design points
#' @param var Variable names
#' @param transform Transformation functions
#' @param minmax Whether to apply min-max normalization
init_design <- function(N = 400, var = NULL, transform = NULL, minmax = TRUE) {
  log_message(paste("Creating experimental design with", N, "points"))
  
  # Setup data folder
  data_folder_name <- get_results_path(paste0("exp", exp.env$exp_id, "/driving_data/"))
  if (!dir.exists(data_folder_name)) {
    dir.create(data_folder_name, recursive = TRUE)
  }
  
  chess_scape <- exp.env$chess_scape
  chess_scape <- chess_scape %>% drop_na()
  
  # Save coordinates
  coord <- chess_scape[, c('x', 'y')]
  saveRDS(coord, file.path(data_folder_name, "chess_scape_bng.rds"))
  assign('chess_scape_xy', coord, exp.env)
  
  # Setup variables
  if (is.null(var)) {
    design_var <- chess_scape %>% select(-any_of(c("x", "y", "z")))
    var <- names(design_var)
  } else {
    design_var <- chess_scape %>% select(all_of(var))
  }
  
  n_var <- length(var)
  
  # Setup transformations
  if (is.null(transform)) {
    transform <- rep(list(identity), n_var)
  }
  
  if (n_var != length(transform)) {
    stop("The supplied 'transform' has different length to the variables considered.")
  }
  
  # Apply transformations
  design_var_transformed <- list()
  for (i in 1:n_var) {
    design_var_transformed[[var[i]]] <- do.call(transform[[i]], list(design_var[[var[i]]]))
  }
  
  # Apply min-max normalization if requested
  if (minmax) {
    design_var_transformed_minmax <- list()
    for (i in 1:n_var) {
      vals <- design_var_transformed[[var[i]]]
      design_var_transformed_minmax[[var[i]]] <- (vals - min(vals)) / (max(vals) - min(vals))
    }
  } else {
    design_var_transformed_minmax <- design_var_transformed
  }
  
  design_var_transformed_minmax <- t(do.call(rbind, design_var_transformed_minmax))
  
  # Use cLHS for space-filling design
  indices <- clhs::clhs(as.data.frame(design_var_transformed_minmax), 
                       size = N, progress = FALSE, simple = TRUE)
  
  selected_design <- design_var_transformed_minmax[indices, , drop = FALSE]
  
  # Save design variables for later use
  saveRDS(design_var, file.path(data_folder_name, "chess_scape_vars.rds"))
  saveRDS(chess_scape, file.path(data_folder_name, "chess_scape.rds"))
  
  log_message(paste("Design created with", nrow(selected_design), "points"))
  return(selected_design)
}

#' Generate coordinate information for design points
#' @param x Design matrix
generate_xy <- function(x) {
  x_s <- apply(x, 1, paste, collapse = ", ")
  idx <- unlist(sapply(x_s, function(item) which(exp.env$chess_scape_lookup_table == item)))
  
  run_coord <- exp.env$chess_scape_xy[idx, , drop = FALSE]
  
  run_coord_labelled <- run_coord %>% 
    mutate(
      scenario = exp.env$rcp,
      ensemble = exp.env$case,
      id = paste("JULESAPPNEW", x, y, sep = "_")
    )
  
  return(run_coord_labelled)
}

# =============================================================================
# DATA INSPECTION FUNCTIONS  
# =============================================================================

#' Inspect data distribution
#' @param var Variable names
#' @param transform Transformation functions
#' @param N Number of samples for inspection
#' @param regions Region filter (optional)
#' @param type Type of inspection ('hist' or 'pair')
inspect_data <- function(var, transform = NULL, N = 1000, regions = NULL, type = 'hist') {
  log_message(paste("Inspecting data:", type, "plot"))
  
  chess_scape <- exp.env$chess_scape
  
  # Apply region filter if specified
  if (!is.null(regions)) {
    chess_scape <- chess_scape %>% filter(region %in% regions)
  }
  
  # Sample data
  if (N < nrow(chess_scape)) {
    chess_scape <- chess_scape %>% sample_n(N)
  }
  
  design_var <- chess_scape %>% select(all_of(var))
  
  # Apply transformations
  if (!is.null(transform)) {
    for (i in 1:length(var)) {
      design_var[[var[i]]] <- do.call(transform[[i]], list(design_var[[var[i]]]))
    }
  }
  
  if (CONFIG$create_plots) {
    if (type == 'hist') {
      # Create histograms
      p <- design_var %>%
        gather(key = "variable", value = "value") %>%
        ggplot(aes(x = value)) +
        geom_histogram(bins = 30, alpha = 0.7) +
        facet_wrap(~variable, scales = "free") +
        theme_minimal() +
        ggtitle("Variable Distributions")
      print(p)
    } else if (type == 'pair') {
      # Create pairs plot (simplified version)
      if (ncol(design_var) <= 6) {  # Only for manageable number of variables
        pairs(design_var, main = "Pairwise Relationships")
      } else {
        log_message("Too many variables for pairs plot, skipping")
      }
    }
  }
  
  return(design_var)
}

#' Check design quality
#' @param x Design matrix
#' @param N Number of points for comparison
#' @param type Type of check ('pair' or 'map')
#' @param regions Region filter (optional)
check_design <- function(x, N = 500, type = 'pair', regions = NULL) {
  log_message(paste("Checking design quality:", type))
  
  if (CONFIG$create_plots) {
    if (type == 'pair' && ncol(x) <= 6) {
      pairs(x, main = "Design Points")
    } else if (type == 'map') {
      # Create map visualization if coordinates available
      coords <- generate_xy(x)
      if (nrow(coords) > 0) {
        plot(coords$x, coords$y, main = "Spatial Distribution", 
             xlab = "X Coordinate", ylab = "Y Coordinate", pch = 16)
      }
    }
  }
  
  return(invisible(x))
}

# =============================================================================
# JULES MODEL INTERFACE
# =============================================================================

#' Run JULES model simulations
#' @param x Design matrix
#' @param years Year values for time-varying simulations
#' @param dr Dimensionality reduction object (optional)
#' @param overwrite Whether to overwrite existing files
#' @param output_dim Output dimensions
#' @param training Whether this is for training (vs testing)
jules <- function(x, years = NULL, dr = NULL, overwrite = FALSE, 
                 output_dim = NULL, training = TRUE, ...) {
  
  log_message("Running JULES simulations...")
  
  master_folder_name <- getwd()
  on.exit(setwd(master_folder_name))
  
  # Setup folders
  exp_folder_name <- get_results_path(paste0("exp", exp.env$exp_id, "/"))
  if (!dir.exists(exp_folder_name)) {
    dir.create(exp_folder_name, recursive = TRUE)
  }
  
  if (training) {
    ensemble_folder_name <- file.path(exp_folder_name, "ensemble/")
  } else {
    ensemble_folder_name <- file.path(exp_folder_name, "ensemble_test/")
  }
  
  if (!dir.exists(ensemble_folder_name)) {
    dir.create(ensemble_folder_name, recursive = TRUE)
  }
  
  setwd(ensemble_folder_name)
  
  if (overwrite) {
    system("rm -f ./JULES7.0-RED1.1*")
  }
  
  # Generate coordinates and run simulations
  run_coord_labelled <- generate_xy(x)
  ids <- run_coord_labelled$id
  
  # Split into manageable chunks
  run_coord_labelled_split <- split(run_coord_labelled, 
                                   (seq(nrow(run_coord_labelled)) - 1) %/% 400)
  
  for (i in 1:length(run_coord_labelled_split)) {
    chunk_data <- run_coord_labelled_split[[i]]
    csv_name <- paste0("single_coordinates_", i, ".csv")
    write.csv(chunk_data, file = csv_name, row.names = FALSE, quote = FALSE)
    
    log_message(paste("Running JULES for chunk", i, "of", length(run_coord_labelled_split)))
    
    # Run JULES (this would call your actual JULES executable)
    # run_jules(csv_name)  # Uncomment when ready to run actual simulations
  }
  
  # Extract results
  log_message("Extracting JULES results...")
  output <- extract_jules(ensemble_folder_name, ids)
  
  # Convert to matrix format
  output_matrix <- do.call(rbind, output)
  output_matrix <- output_matrix[!is.na(rowSums(output_matrix)), , drop = FALSE]
  
  log_message(paste("JULES simulation complete:", nrow(output_matrix), "valid results"))
  
  return(output_matrix)
}

#' Extract JULES results from NetCDF files
#' @param folder_path Path to folder containing NetCDF files
#' @param ids List of simulation IDs to extract
extract_jules <- function(folder_path, ids) {
  log_message("Extracting JULES output data...")
  
  files <- list.files(folder_path)
  nc_files <- files[grepl("JULES7.0-RED1.1", files, fixed = TRUE)]
  
  y_train <- vector('list', length(ids))
  
  for (i in 1:length(ids)) {
    matching_idx <- grep(paste0(ids[i]), nc_files, fixed = TRUE)
    
    if (length(matching_idx) > 0) {
      matching_file <- nc_files[matching_idx[1]]  # Take first match
      
      tryCatch({
        output_dat <- extract_one_jules(file.path(folder_path, matching_file))
        if (!is.null(output_dat)) {
          y_train[[i]] <- output_dat$carbon
        }
      }, error = function(e) {
        log_message(paste("Error extracting", matching_file, ":", e$message), "WARN")
      })
    }
  }
  
  # Remove NULL entries
  y_train <- y_train[!sapply(y_train, is.null)]
  
  log_message(paste("Successfully extracted", length(y_train), "files"))
  return(y_train)
}

#' Extract data from a single JULES NetCDF file
#' @param file_path Path to NetCDF file
extract_one_jules <- function(file_path) {
  if (!file.exists(file_path)) {
    return(NULL)
  }
  
  tryCatch({
    df <- tidync(file_path)
    
    # Extract vegetation carbon
    temp1 <- activate(df, "D0,D1,D2,D6", select_var = "c_veg") %>%
      hyper_tibble() %>%
      filter(pft == 2) %>%
      select(-pft)
    
    # Extract fractional coverage
    temp2 <- activate(df, "D0,D1,D2,D6", select_var = "frac_CA") %>%
      hyper_tibble() %>%
      filter(pft == 2) %>%
      select(-pft)
    
    # Calculate carbon storage (tCO2e/ha)
    carbon <- temp1$c_veg * temp2$frac_CA * 44 * 10 / 12
    
    result <- data.frame(
      carbon = carbon,
      id = str_extract(file_path, "JULESAPPNEW_[0-9]+_[0-9]+")
    )
    
    return(result)
    
  }, error = function(e) {
    log_message(paste("Error reading NetCDF file:", e$message), "WARN")
    return(NULL)
  })
}

# =============================================================================
# EMULATION FUNCTIONS
# =============================================================================

#' Create and validate emulator
#' @param x Input design matrix
#' @param y Output data matrix
#' @param pca_components Number of PCA components
#' @param kernel_names Kernel names for DGP
#' @param use_vecchia Whether to use Vecchia approximation
create_emulator <- function(x, y, pca_components = 1, 
                           kernel_names = c("matern2.5", "sexp"), 
                           use_vecchia = TRUE) {
  
  log_message("Creating emulator with PCA preprocessing...")
  
  # Initialize Python sklearn
  sklearn <- reticulate::import('sklearn')
  
  # Apply PCA
  pca <- sklearn$decomposition$PCA(as.integer(pca_components))
  y_basis <- pca$fit_transform(y)
  
  # Save PCA object
  pca_file <- get_results_path(paste0('pca_obj_', exp.env$exp_id, '.pkl'))
  dgpsi_py <- reticulate::import('dgpsi')
  dgpsi_py$write(pca, pca_file)
  
  log_message("Building DGP emulator...")
  set_seed(CONFIG$validation_seed)
  
  emulator_id <- paste0('jules-', exp.env$rcp, '-', exp.env$case)
  
  if (length(kernel_names) == 1) {
    m <- dgp(x, y_basis, id = emulator_id, name = kernel_names[1], vecchia = use_vecchia)
  } else {
    m <- dgp(x, y_basis, id = emulator_id, name = kernel_names, vecchia = use_vecchia)
  }
  
  # Validate emulator
  log_message("Validating emulator...")
  m <- validate(m)
  
  if (CONFIG$create_plots) {
    plot(m)
  }
  
  # Save emulator
  emulator_file <- get_results_path(paste0('emulator_', exp.env$exp_id))
  write(m, emulator_file)
  
  log_message("Emulator created and validated successfully")
  
  return(list(emulator = m, pca = pca))
}

#' Reconstruct data from PCA components
#' @param pca_obj PCA object
#' @param x PCA-transformed data
reconstruct_from_pca <- function(pca_obj, x) {
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  return(x %*% pca_obj$components_ + 
         matrix(rep(as.vector(pca_obj$mean_), nrow(x)), 
                nrow = nrow(x), byrow = TRUE))
}

#' Calculate NRMSE for out-of-sample testing
#' @param emulator Trained emulator object
#' @param inverse_dr Inverse dimensionality reduction function
#' @param testing_x Test input data
#' @param testing_y Test output data
#' @param workers Number of workers for parallel processing
nrmse_oos <- function(emulator, inverse_dr, testing_x, testing_y, workers = 1) {
  log_message("Calculating out-of-sample NRMSE...")
  
  predictions <- predict(emulator, x = testing_x, cores = workers)
  oos_mean_pca <- predictions$results$mean
  oos_mean_original <- do.call(inverse_dr, list(oos_mean_pca))
  
  true_original <- testing_y
  nrmse <- mean(sqrt(rowMeans((oos_mean_original - true_original)^2)) / 
                (apply(true_original, 1, max) - apply(true_original, 1, min)))
  
  log_message(paste("Out-of-sample NRMSE:", round(nrmse, 4)))
  return(nrmse)
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Set random seed consistently across packages
#' @param seed Random seed value
set_seed <- function(seed) {
  set.seed(seed)
  # Also set Python random seed if using reticulate
  if (reticulate::py_available()) {
    reticulate::py_run_string(paste0("import random; random.seed(", seed, ")"))
    reticulate::py_run_string(paste0("import numpy as np; np.random.seed(", seed, ")"))
  }
}

#' Initialize test design
#' @param N Number of test points
init_test <- function(N = 100) {
  log_message(paste("Creating test design with", N, "points"))
  
  chess_scape <- exp.env$chess_scape %>% drop_na()
  design_var <- chess_scape %>% select(all_of(CONFIG$variables))
  
  # Apply transformations and normalization (similar to init_design)
  design_var_transformed <- list()
  for (i in 1:length(CONFIG$variables)) {
    var_name <- CONFIG$variables[i]
    transform_fun <- CONFIG$transform_functions[[i]]
    design_var_transformed[[var_name]] <- do.call(transform_fun, list(design_var[[var_name]]))
  }
  
  # Min-max normalization
  design_var_transformed_minmax <- list()
  for (i in 1:length(CONFIG$variables)) {
    var_name <- CONFIG$variables[i]
    vals <- design_var_transformed[[var_name]]
    design_var_transformed_minmax[[var_name]] <- (vals - min(vals)) / (max(vals) - min(vals))
  }
  
  design_var_transformed_minmax <- t(do.call(rbind, design_var_transformed_minmax))
  
  # Sample test points
  test_indices <- sample(nrow(design_var_transformed_minmax), N)
  test_design <- design_var_transformed_minmax[test_indices, , drop = FALSE]
  
  return(test_design)
}

# =============================================================================
# INITIALIZATION
# =============================================================================

log_message("JULES utilities loaded successfully")
log_message("Essential functions available for JULES emulation workflow") 