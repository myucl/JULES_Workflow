# JULES Clean - Streamlined JULES Model Emulation Workflow

## JULES Set–up for ADD-TREES (Intro)

Source: `https://code.metoffice.gov.uk/trac/jules/wiki/RoseJULESonJASMIN`

### 1. To get access to JULES code

- Use the JULES code access link. You will get a screen to fill in your information, accept Terms and Conditions, and click the “Send request” button. If you have issues, read carefully the information in the rectangle of the image on that page.
- Once you get access to the JULES code you will get Met Office Science Repository Service (MOSRS) password advice at “Changing your MOSRS password”.

### 2. To get a JASMIN account

- Apply for a new JASMIN account: Sign up for a JASMIN login account (Sign in with a JASMIN account)
- Create an SSH key pair: see linux/mac instructions for `ssh-keygen`
- Upload your public key: sign in to your JASMIN account. On your JASMIN Profile there will be a blue “Update Key” button, where you can add your public key (registering your ssh key)
- Request access to the `uknetzero` Group Workspace: on your JASMIN Profile click on the green “My Services” button where you can request access to the workspace.
- If you have any problems please see more information and login problems, or send an e-mail to the JASMIN Helpdesk (`suppor@jasmin.ac.uk`).

### 3. Software

Install the following software in advance depending on your operating system; you will need one of these to be able to connect to the JASMIN-HPC.

- Windows: install MobaXterm and configure your computer for MobaXterm following the instructions in this document and instructions in section 4.
- Macintosh: install XQuartz.
- Linux: follow the instructions described in section 4.

### 4. Configuration of your account

1) Configure your SSH config file to use the `login-01.jasmin` server to tunnel through to `cylc2.jasmin` where the suite will run.
2) Edit your `~/.ssh/config` file as follows (changing the `<jasmin_userid>` to your own userid).
3) Create the `~/.ssh` folder on your terminal using:

```bash
mkdir -p ~/.ssh
```

4) Edit the config file:

```bash
vim ~/.ssh/config
```

5) Add the following lines to the file (replace `jasmin_userid`):

```bash
Host *
    ServerAliveInterval 30

Host jlogin1
    Hostname login-01.jasmin.ac.uk
    User jasmin_userid
    ForwardAgent yes

Host xfer-vm-0?
    Hostname %h.jasmin.ac.uk
    User jasmin_userid
    ForwardAgent yes

Host sci-vm-0? sci-ph-0? cylc2
    HostName %h.jasmin.ac.uk
    User jasmin_userid
    ForwardAgent yes
    ProxyCommand ssh -Y jlogin1 -W %h:%p
    ProxyJump jasmin_userid@jlogin1
```

6) Start an ssh-agent session on the terminal and add your key:

```bash
eval $(ssh-agent -s); ssh-add ~/.ssh/id_rsa_jasmin
```

The `~/.ssh/id_rsa_jasmin` file is your JASMIN SSH private key; enter your passphrase when prompted.

7) Log into JASMIN with the following command (example userid):

```bash
ssh -A -X kday002@login-01.jasmin.ac.uk
```

Make sure to use the `-X` option for X-forwarding of graphical applications. If it works, you will find yourself logged into `login-01.jasmin`.

8) Modify the following configuration files on your JASMIN account on the `cylc2.jasmin` server

Ensure that you login to `cylc2.jasmin` — do not try from other systems since they do not have the required setup. If you try instead to submit Cylc jobs from the `sci*.jasmin` servers, the Cylc8 GUI might not appear. If you want to do other (interactive) processing besides running Rose/Cylc, then it is better to use the `sci*.jasmin` servers. The old path of `/apps/contrib/metomi/` has been changed to `/apps/jasmin/metomi/` in several places below.

9) Edit your `~/.bash_profile` by adding the following lines to the top of the file:

```bash
# Get the aliases and functions
if [ -f ~/.bashrc ]; then
        . ~/.bashrc
fi

# Provide access to FCM, Rose and Cylc
PATH=$PATH:/apps/jasmin/metomi/bin

# User specific environment and startup programs
export PATH=$PATH:$HOME/bin
HOST=$(hostname)

if [[ $HOST = "cylc2.jasmin.ac.uk" ]]; then
# Rose/cylc on  cylc2.jasmin node
export PATH=$PATH:/apps/jasmin/metomi/bin
fi
```

10) Add the following at the end of your `~/.bashrc` file:

```bash
[[ $- != *i* ]] && return # Stop here if not running interactively
[[ $(hostname) = "cylc2.jasmin.ac.uk" ]] && .mosrs-setup-gpg-agent
# Enable bash completion for Rose commands
[[ -f /apps/jasmin/metomi/rose/etc/rose-bash-completion ]] && . /apps/jasmin/metomi/rose/etc/rose-bash-completion
```

Now, whenever you login to `cylc2.jasmin` you should be prompted for your Met Office Science Repository Service password. Just press return if you do not want access to the repositories and you should not be prompted again during that login session.

11) A further set-up for JASMIN and MOSRS requires an update to your `~/.subversion/servers` file.

- Create the folder `~/.subversion` using:

```bash
mkdir -p ~/.subversion
```

- Edit the file `~/.subversion/servers`:

```bash
vi ~/.subversion/servers
```

- Add the following lines (use your MOSRS username instead of `myusername`):

```ini
[groups]
metofficesharedrepos = code*.metoffice.gov.uk

[metofficesharedrepos]
# Specify your Science Repository Service username here
username = myusername
store-plaintext-passwords = no
```

12) Open `~/.subversion/config` and comment (by adding the `#` symbol) any lines starting with:

```ini
#password-stores =
```

13)

a) Create the folder for Rose config:

```bash
mkdir -p ~/.metomi
```

b) Add the following lines to the `~/.metomi/rose.conf` file (change `myusername` to your MOSRS username):

```ini
[rosie-id]
prefix-default=u
prefix-location.u=https://code.metoffice.gov.uk/svn/roses-u
prefix-username.u=myusername
#username is all in lower case

prefix-ws.u=https://code.metoffice.gov.uk/rosie/u

[rose-stem]
automatic-options=SITE=jasmin
```

14) Check the Rose configuration by running:

```bash
rose config
```

15)

a) Create the following folder:

```bash
mkdir -p ~/.metomi/fcm
```

b) Open the file `~/.metomi/fcm/keyword.cfg` and add the lines:

```ini
location{primary, type:svn}[jules.x] = https://code.metoffice.gov.uk/svn/jules/main
browser.loc-tmpl[jules.x] = https://code.metoffice.gov.uk/trac/{1}/intertrac/source:/{2}{3}
browser.comp-pat[jules.x] = (?msx-i:\A // [^/]+ /svn/ ([^/]+) /*(.*) \z)

location{primary, type:svn}[jules.xm] = https://code.metoffice.gov.uk/svn/jules/main
browser.loc-tmpl[jules.xm] = https://code.metoffice.gov.uk/trac/{1}/intertrac/source:/{2}{3}
browser.comp-pat[jules.xm] = (?msx-i:\A // [^/]+ /svn/ ([^/]+) /*(.*) \z)

location{primary, type:svn}[jules_doc.x] = https://code.metoffice.gov.uk/svn/jules/doc
browser.loc-tmpl[jules_doc.x] = https://code.metoffice.gov.uk/trac/{1}/intertrac/source:/{2}{3}
browser.comp-pat[jules_doc.x] = (?msx-i:\A // [^/]+ /svn/ ([^/]+) /*(.*) \z)
```

16) Test if the connection to `cylc2.jasmin` works (from `login-01.jasmin`):

```bash
ssh -AX kday002@cylc2.jasmin.ac.uk
xterm &
```

An xterm window will appear, which you can close if you want. Exit the `cylc2.jasmin` server using `logout` or Control-D.

17) Test that the direct connection to `cylc2.jasmin` works (from your own laptop or desktop machine) using:

```bash
ssh -AX kday002@cylc2.jasmin.ac.uk
```

Then, login directly from your own machine to `cylc2.jasmin`.

You should only have to enter your MOSRS password when prompted; there is no login password for `cylc2.jasmin`.

```bash
xterm &
```

An xterm window will appear, which you can close. Then exit (logout) the `cylc2.jasmin` server.

18) To run JULES you will need a suite-ID which contains all the parameters to run the model. Contact your tutor to get access to a suite-ID.

### 5. FAIR USE POLICY

We ask that users of these instructions and suite-ID follow the guidelines for fair use.

The guidelines apply to research applications that result in publication (journal articles, thesis, technical reports).

The creators of these instructions and the suite-ID are Carolina Duran Rojas and Mingda Yuan (Ming Deyu, Dany Williamson and Anna Harper).

- Offer co-authorship to the owners and contributors to this document if they have provided technical or scientific support to the analysis.
- Let the authors of this document and the suite-ID know of plans to use the suite for publication at the early stages, so we can provide updates to run the model on JASMIN.

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

