#Author: Rogelio Aguilar
#Project: Spatial Dietary Project


#Set working directory
setwd("/home/rogaguilar/Spatial/data/Rose_Li_VisiumHD")

#Functions####
#Input: Vector/list of packages to check system for
#Description:Function to check and install packages
install_missing_packages <- function(package_list) {
  # Load BiocManager if not already loaded
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    stop("BiocManager is not installed.")
  }
  
  # Load remotes if it's not already installed (for GitHub installations)
  if (!requireNamespace("remotes", quietly = TRUE)) {
    message("remotes package is not installed. Installing remotes from CRAN...")
    install.packages("remotes")
  }
  
  # Loop through the package list
  for (package in package_list) {
    # Check if the package is installed
    if (!(package %in% rownames(installed.packages()))) {
      message(paste("Package", package, "is not installed. Checking availability..."))
      
      # Special check for SeuratData and SeuratWrappers (install from GitHub)
      if (package %in% c("SeuratData", "SeuratWrappers")) {
        message(paste("Installing", package, "from GitHub..."))
        if (package == "SeuratData") {
          package <- "seurat-data"
        } else if (package == "SeuratWrappers") {
          package <- "seurat-wrappers"
        }
        remotes::install_github(paste("satijalab", package, sep = "/"))
      } else {
        # Check if the package is available in Bioconductor
        if (package %in% BiocManager::available()) {
          message(paste("Installing", package, "from Bioconductor..."))
          BiocManager::install(package)
        } else {
          # Check if the package is available on CRAN
          if (package %in% rownames(available.packages())) {
            message(paste("Package", package, "not found in Bioconductor. Installing from CRAN..."))
            install.packages(package)
          } else {
            # If not available in either Bioconductor or CRAN
            message(paste("Package", package, "is not available from either Bioconductor or CRAN."))
          }
        }
      }
    } else {
      message(paste("Package", package, "is already installed."))
    }
  }
}


#-------------------------------Script Start-----------------------------------
#Load packages and CSVs
packages <- c("R.utils","hdf5r","arrow",
              "future","Seurat","SeuratData","SeuratWrappers","Banksy",
              "tidyverse","patchwork")

install_missing_packages(packages)

lapply(packages,library, character.only=TRUE)
