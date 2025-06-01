install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) {
    install.packages(new_packages)
  }
}


cran_packages <- c(
  # Core packages
  "tidyverse", # Data manipulation and visualization
  "here", # Path management

  # Statistical packages
  "Mass", # Negative binomial distribution

  # Data Manipulation
  "data.table", # fast data operations
  "zoo", # Rolling windows

  # Visualization
  "ggplot2",
  "gridExtra", # Multiple plot layout
  "ggpbur", # Publication-ready plots
  "RColorBrewer", # Color palettes
  "virdis", # Color-blind friendly palettes
  "Scales", # scale function for plots

  # Reporting
  "knitr", # Dynamic reports
  "rmarkdown", # Markdown documents
  "DT", # Interactive tables

  # Dev tools
  "devtools", # Read YAML config files
  "testthat", # Unit testing
  "profvis", "Code profiling"

)

# Install CRAN packages
cat("Installing CRAN packages...\n")
install_if_missing(cran_packages)

# Load all packages to verify installation
cat("\nVerifying package installation...\n")
successfully_loaded <- character()
failed_to_load <- character()

for(pkg in cran_packages) {
  if(require(pkg, character.only = TRUE, quietly = TRUE)) {
    successfully_loaded <- c(successfully_loaded, pkg)
  } else {
    failed_to_load <- c(failed_to_load, pkg)
  }
}

  # Report results
  cat("\n=== Installation Summary ===\n")
  cat("Successfully installed:", length(successfully_loaded), "packages\n")
  if(length(failed_to_load) > 0) {
    cat("Failed to load:", paste(failed_to_load, collapse = ", "), "\n")
  } else {
    cat("All packages loaded successfully!\n")
  }


