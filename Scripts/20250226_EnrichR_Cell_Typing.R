#Author: Rogelio Aguilar
#Project: Spatial Dietary Project
#Objective: Plot cell type enrichment from enrichR databases 


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

#Input: Marker CSV, databases to profile as a list, number on slide, sample being profiled
#Output: Statement of plotting completion
#Description:
enrichment.plotting <- function(df, databases, slide_number,sample_name) {
  #Gather number of clusters
  clusters <- unique(df$cluster)
  
  #Check if dataframe is empty
  if (nrow(df) == 0) {
    cat("Dataframe is empty for this sample.\n")
    return(NULL)  # Exit function if df is empty
  }
  #Create a list of plots for each database
  database.plots <- lapply(databases, function(db) {
    #Check if database has markers for sample
    clust.plots <- list()
    for (clust in clusters) {
      cluster.markers <- df %>%
        dplyr::filter(avg_log2FC > 1) %>%
        filter(cluster==clust)
      if (length(cluster.markers$gene) > 0) {
        #Apply enrichment through accessing enrichR API
        enriched <- enrichR::enrichr(genes = cluster.markers$gene, databases = db)
        
        # Check if enriched[[db]] is NULL or empty
        if (is.null(enriched[[db]]) || nrow(enriched[[db]]) == 0) {
          return(NULL)  # Skip plotting if enrichment result is empty or NULL
        }
        # Plot the enrichment for the top 20 terms ordered by p-value
        plot <- plotEnrich(enriched[[db]], showTerms = 20, numChar = 40,
                           y = "Ratio", orderBy = "P.value",
                           title = paste("Clust:", clust, sep = " "))
        clust.plots[[clust+1]] <- plot        
      }
    }
    # Check if enriched list is empty
    if (length(clust.plots) == 0) {
      return(NULL)  # Skip saving if enrichment plot list is empty
    }
    
    # Save the plots as a PDF
    pdf(file = paste("./Analysis/BANKSY_Normalized_QC_Filtered/Plots/Enrichment/",slide_number,"_",
                     sample_name,"_",db,"_enrichment_plots.pdf",sep = ""),
        width = 24, height = 12)
    plot.grid(plots_list = na.omit(clust.plots))
    dev.off()
    
    # Output message
    message <- paste("Saved plots for", sample_name, "\n", sep = "")    
    return(cat(message))     
  })
}

#Input: Marker CSV, databases to profile as a list, number on slide, sample being profiled
#Output: Statement of file storage completion
#Description:
enrichment.excel <- function(df, databases, slide_number,sample_name) {
  #Gather number of clusters
  clusters <- unique(df$cluster)
  for (clust in clusters) {
    cluster.markers <- df %>%
      dplyr::filter(avg_log2FC > 1) %>%
      filter(cluster==clust)
    if (length(cluster.markers$gene) > 0) {
      #Apply enrichment through accessing enrichR API
      enriched <- enrichR::enrichr(genes = cluster.markers$gene, databases = databases)
      
      # Check if enriched list is empty
      if (length(enriched) == 0) {
        return(NULL)  # Skip saving if enrichment result is empty
      }
      #Store enrichment into Excel file; each database profiled will appear as a tab 
      printEnrich(data = enriched, outFile = "excel",
                  prefix = paste("./Analysis/BANKSY_Normalized_QC_Filtered/Cluster_Markers/Enrichment/",
                                 slide_number,"_", sample_name,"_cluster_",clust,
                                 "_enrichment",
                                 sep = "")) #name Excel file will be stored as
    }
  }
  message <- paste("Saved excel for: ", slide_number, "_", sample_name,"_cluster_",clust,"\n",sep = "")
  return(cat(message))
  cat("Excel files have been saved \n")    
}

#Input: List of plots
#Description: Function to concatenate 4 ggplots into 1 plot
# Function to create a grid of 4 plots
plot.grid <- function(plots_list) {
  # Remove NULL values from the plots list
  plots_list <- Filter(Negate(is.null), plots_list)
  # Loop through the plots list and group them into sets of 4
  for (i in seq(1, length(plots_list), by = 4)) {
    # Get the subset of 4 plots
    plot_subset <- plots_list[i:min(i+3, length(plots_list))]
    
    # Combine the 4 plots using wrap_plots() and print them
    combined_plot <- wrap_plots(plot_subset, ncol = 2, nrow = 2)
    
    # Print the combined plot
    print(combined_plot)
  }
}



#-------------------------------Script Start-----------------------------------
#Load packages and CSVs
packages <- c("tidyverse","patchwork", "enrichR","org.Mm.eg.db")

#Check if packages needed are installed
install_missing_packages(packages)

#Load packages
lapply(packages,library, character.only=TRUE)

sample.csv <- read.csv(file = "./Visium_Mouse_Tumor_MycCap_Slide_Samples.csv")

slides <- unique(sample.csv$Slide_Number)

output.dir <- "./Analysis/BANKSY_Normalized_QC_Filtered/Plots/Enrichment"

if (!dir.exists(output.dir)) {
  dir.create(output.dir, recursive = TRUE)
}

#Enrichment####
lapply(slides, function(slide) {
  sample.df <- sample.csv %>% filter(Slide_Number == slide)
  #Enrichment of sample
  lapply(sample.df$Sample, function(sample_name) {
    #Load in marker csv file
    markers.csv <- read.csv(file = paste("./Analysis/BANKSY_Normalized_QC_Filtered/Cluster_Markers/",
                                         slide,"_", sample_name,"_all_markers.csv",
                                         sep = "")) %>%
      mutate(gene = gsub(pattern = "\\.m0$",replacement = "", x = `gene`))
    
    #Databases to search for enrichment in
    databases <- c("Azimuth_Cell_Types_2021", "Tabula_Muris", "Cancer_Cell_Line_Encyclopedia")
    
    
    #Plot enrichment for selected databases
    enrichment.plotting(df = markers.csv, 
                        databases = databases, 
                        slide_number = slide,
                        sample_name = sample_name)
    #Write excel file for enrichment of selected databases
    enrichment.excel(df = markers.csv, 
                     databases = databases, 
                     slide_number = slide,
                     sample_name = sample_name)
    # Print completion message for the slide
    return(print(paste("Slide", slide, "complete for sample", sample_name)))
  })
})