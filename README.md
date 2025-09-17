# JULES Clean - Streamlined JULES Model Emulation Workflow

This is a cleaned and restructured version of the JULES model emulation workflow, incorporating improvements from both the original `master.R` and the advanced `RED_CLIBRATION_DF.R` approaches.

## Purpose

This clean version focuses on the essential components for JULES model emulation:
- Experimental design using Latin Hypercube Sampling
- JULES model simulation interface
- Gaussian Process emulation with PCA preprocessing
- Time-varying analysis capabilities
- Model validation and prediction

##  Project Structure

```
jules_clean/
├── config/
│   └── jules_config.R          # Centralized configuration
├── scripts/
│   ├── master_clean.R          # Main workflow script
│   └── jules_utils.R           # Utility functions
├── data/
│   ├── design_data_rcp_85_15.rds    # Design data (example)
│   └── (other data files...)
├── lib/
│   ├── calibration_ensembleREDHARP.csv  # Calibration data
│   └── (other library files...)
├── results/
│   └── (generated results...)
└── README.md                   # This file
```

## Quick Start

### 1. Prerequisites

Ensure you have the following R packages installed:
```r
install.packages(c("readr", "tidync", "tidyr", "tidyverse", "ggplot2", 
                   "reticulate", "clhs", "feather"))

# For emulation (install dgpsi from GitHub)
devtools::install_github("mingdeyu/dgpsi")
```

### 2. Configuration

Edit `config/jules_config.R` to set up your specific scenario:
```r
CONFIG <- list(
  rcp = '85',                    # Your RCP scenario
  ensemble = '15',               # Your ensemble ID
  n_samples = 400,               # Number of design points
  variables = c("oneoveralpha", "hcon", "vcrit", ...),  # Your variables
  # ... other settings
)
```

### 3. Run the Workflow

```r
# Navigate to the scripts directory
setwd("path/to/jules_clean/scripts")

# Run the main workflow
source("master_clean.R")
```

##  Key Features

### Structured Configuration
- Centralized configuration in `jules_config.R`
- Easy parameter modification
- Automatic path management
- Validation functions

### Clean Workflow Steps
1. **Data Loading**: Load design data and initialize environment
2. **Data Exploration**: Inspect variable distributions
3. **Experimental Design**: Create space-filling design using cLHS
4. **Time Setup**: Configure time-varying analysis (optional)
5. **JULES Simulation**: Run model simulations
6. **Emulation**: Build GP emulator with PCA preprocessing
7. **Time-Varying Analysis**: Extended time-series emulation
8. **Validation**: Test emulator performance
9. **Prediction**: Generate predictions on new data
10. **Visualization**: Create diagnostic plots
11. **Summary**: Generate workflow summary

### Essential Functions
- `init_exp()`: Initialize experiment environment
- `init_design()`: Create experimental design
- `jules()`: Interface to JULES model
- `create_emulator()`: Build and validate emulator
- `inspect_data()`: Data exploration tools
- `extract_jules()`: Process JULES NetCDF outputs

## 🔧 Configuration Options

### Basic Settings
```r
CONFIG$rcp = '85'              # RCP scenario
CONFIG$ensemble = '15'         # Ensemble identifier
CONFIG$n_samples = 400         # Design points
CONFIG$random_seed = 50        # Reproducibility
```

### Advanced Settings
```r
CONFIG$use_time_varying = TRUE      # Enable time-varying analysis
CONFIG$pca_components = 1           # PCA dimensions
CONFIG$emulator_names = c("matern2.5", "sexp")  # GP kernels
CONFIG$create_plots = TRUE          # Generate visualizations
```

## 📈 Outputs

The workflow generates:
- **Emulators**: Trained GP models saved as `.pkl` files
- **Predictions**: Spatial predictions in `.feather` format
- **Diagnostics**: Validation plots and metrics
- **Summary**: Workflow summary with key statistics


##  Usage Notes

### For New Users:
1. Start by reviewing `jules_config.R`
2. Ensure your data files are in the correct locations
3. Run the validation functions first
4. Execute the workflow step by step

### For Existing Users:
1. Your existing data files should work with minimal changes
2. Update file paths in the configuration
3. The core functionality remains the same
4. Advanced features (like sequential design) are noted but not implemented in this clean version

### Data Requirements:
- **Design data**: `.rds` file with spatial and environmental variables
- **Calibration data**: `.csv` files with parameter ranges
- **JULES outputs**: NetCDF files (for actual runs)

## ️ Customization

### Adding New Variables:
```r
CONFIG$variables <- c("your_var1", "your_var2", ...)
CONFIG$transforms <- c("identity", "log10", ...)
```

### Changing Emulation Settings:
```r
CONFIG$pca_components <- 2
CONFIG$emulator_names <- c("matern2.5")
CONFIG$use_vecchia <- FALSE
```

### Modifying Paths:
```r
CONFIG$base_folder <- "/your/path/to/jules_clean"
CONFIG$data_folder <- "your_data_folder"
```

##  Troubleshooting

### Common Issues:
1. **Missing data files**: Check file paths in configuration
2. **Python/reticulate errors**: Ensure Python environment is set up correctly
3. **Memory issues**: Reduce `n_samples` or use `use_vecchia = TRUE`
4. **Plot errors**: Set `create_plots = FALSE` if running headless

### Getting Help:
1. Check the configuration validation output
2. Review log messages for specific errors
3. Ensure all required packages are installed
4. Verify data file formats match expectations

