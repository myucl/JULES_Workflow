# =============================================================================
# JULES Clean - Setup Script
# Description: Check environment and install required packages
# =============================================================================

cat("=== JULES Clean Setup Script ===\n")
cat("This script will check your R environment and install required packages.\n\n")

# =============================================================================
# PACKAGE INSTALLATION
# =============================================================================

# Required CRAN packages
required_packages <- c(
  "readr",           # Data reading
  "tidync",          # NetCDF handling
  "tidyr",           # Data tidying
  "tidyverse",       # Data manipulation
  "ggplot2",         # Plotting
  "dplyr",           # Data manipulation
  "reticulate",      # Python interface
  "clhs",            # Latin Hypercube Sampling
  "feather",         # Fast data serialization
  "devtools"         # For installing from GitHub
)

# Special packages (GitHub or special installation)
special_packages <- list(
  dgpsi = "mingdeyu/dgpsi"  # Gaussian Process emulation
)

cat("Checking and installing required packages...\n")

# Install CRAN packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("Installing", pkg, "from CRAN...\n"))
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat(paste("✓", pkg, "already installed\n"))
  }
}

# Install special packages
for (pkg_name in names(special_packages)) {
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    cat(paste("Installing", pkg_name, "from GitHub...\n"))
    devtools::install_github(special_packages[[pkg_name]])
  } else {
    cat(paste("✓", pkg_name, "already installed\n"))
  }
}

# =============================================================================
# PYTHON ENVIRONMENT CHECK
# =============================================================================

cat("\nChecking Python environment...\n")

# Check if Python is available
if (reticulate::py_available()) {
  cat("✓ Python is available\n")
  
  # Check for required Python packages
  required_py_packages <- c("numpy", "sklearn", "dgpsi")
  
  for (py_pkg in required_py_packages) {
    if (reticulate::py_module_available(py_pkg)) {
      cat(paste("✓", py_pkg, "available\n"))
    } else {
      cat(paste("⚠", py_pkg, "not available - may need installation\n"))
      if (py_pkg == "dgpsi") {
        cat("  Note: dgpsi Python package will be installed automatically with R dgpsi\n")
      }
    }
  }
} else {
  cat("⚠ Python not available - some features may not work\n")
  cat("  Consider installing Python and reticulate configuration\n")
}

# =============================================================================
# ENVIRONMENT VALIDATION
# =============================================================================

cat("\nValidating environment...\n")

# Test key functionality
tryCatch({
  # Test data manipulation
  library(dplyr)
  test_df <- data.frame(x = 1:5, y = 6:10)
  test_result <- test_df %>% mutate(z = x + y)
  cat("✓ Data manipulation works\n")
}, error = function(e) {
  cat("✗ Data manipulation test failed:", e$message, "\n")
})

# Test plotting
tryCatch({
  library(ggplot2)
  p <- ggplot(data.frame(x = 1:5, y = 1:5), aes(x, y)) + geom_point()
  cat("✓ Plotting functionality works\n")
}, error = function(e) {
  cat("✗ Plotting test failed:", e$message, "\n")
})

# Test Latin Hypercube Sampling
tryCatch({
  library(clhs)
  test_data <- data.frame(x = rnorm(100), y = rnorm(100))
  sample_idx <- clhs(test_data, size = 10, simple = TRUE)
  cat("✓ Latin Hypercube Sampling works\n")
}, error = function(e) {
  cat("✗ LHS test failed:", e$message, "\n")
})

# Test dgpsi if available
tryCatch({
  library(dgpsi)
  cat("✓ dgpsi package loaded successfully\n")
}, error = function(e) {
  cat("⚠ dgpsi test failed:", e$message, "\n")
  cat("  This may be normal if Python environment needs setup\n")
})

# =============================================================================
# CONFIGURATION CHECK
# =============================================================================

cat("\nChecking configuration...\n")

# Check if config file exists
config_file <- "config/jules_config.R"
if (file.exists(config_file)) {
  cat("✓ Configuration file found\n")
  
  # Try to load config
  tryCatch({
    source(config_file)
    cat("✓ Configuration loaded successfully\n")
    
    # Run validation if available
    if (exists("validate_config")) {
      validate_config()
    }
    
  }, error = function(e) {
    cat("⚠ Configuration loading failed:", e$message, "\n")
  })
} else {
  cat("⚠ Configuration file not found at:", config_file, "\n")
  cat("  Make sure you're running this from the project root directory\n")
}

# =============================================================================
# DATA FILE CHECK
# =============================================================================

cat("\nChecking for essential data files...\n")

# Check data directory
if (dir.exists("data")) {
  data_files <- list.files("data", pattern = "\\.rds$")
  if (length(data_files) > 0) {
    cat("✓ Found", length(data_files), "RDS data files\n")
    for (file in data_files[1:min(3, length(data_files))]) {
      cat("  -", file, "\n")
    }
    if (length(data_files) > 3) {
      cat("  ... and", length(data_files) - 3, "more\n")
    }
  } else {
    cat("⚠ No RDS data files found in data/ directory\n")
  }
} else {
  cat("⚠ data/ directory not found\n")
}

# Check lib directory
if (dir.exists("lib")) {
  lib_files <- list.files("lib", pattern = "\\.csv$")
  if (length(lib_files) > 0) {
    cat("✓ Found", length(lib_files), "CSV files in lib/\n")
  } else {
    cat("⚠ No CSV files found in lib/ directory\n")
  }
} else {
  cat("⚠ lib/ directory not found\n")
}

# =============================================================================
# FINAL RECOMMENDATIONS
# =============================================================================

cat("\n=== Setup Summary ===\n")

cat("Next steps:\n")
cat("1. Review and customize config/jules_config.R for your scenario\n")
cat("2. Ensure your data files are in the data/ and lib/ directories\n")
cat("3. Navigate to scripts/ directory: setwd('scripts')\n")
cat("4. Run the main workflow: source('master_clean.R')\n")
cat("\nFor help, see README.md or check the configuration validation output.\n")

cat("\n=== Setup Complete ===\n")

# =============================================================================
# OPTIONAL: CREATE EXAMPLE CONFIGURATION
# =============================================================================

create_example_config <- function() {
  cat("Creating example configuration...\n")
  
  example_config <- '
# Example configuration for JULES Clean
# Copy this to config/jules_config.R and customize

CONFIG <- list(
  # Basic settings
  rcp = "85",
  ensemble = "15", 
  n_samples = 400,
  
  # Variables (customize for your data)
  variables = c("oneoveralpha", "hcon", "vcrit", "slope", "mean_pr1", 
               "spring_tas1", "winter_tasmin1", "mean_rsds_dec1", "vpd1"),
  
  # File paths (update for your system)
  base_folder = getwd(),
  
  # Other settings
  create_plots = TRUE,
  verbose = TRUE
)
'
  
  if (!dir.exists("config")) {
    dir.create("config")
  }
  
  example_file <- "config/jules_config_example.R"
  writeLines(example_config, example_file)
  cat("Example configuration written to:", example_file, "\n")
}

# Uncomment to create example config
# create_example_config() 