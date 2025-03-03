#Author: Rogelio Aguilar
#Project: Spatial Dietary Project
#Objective: Apply BANKSY normalization to filtered Visium spatial objects


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

#Input: Seurat object for a slide with bin size 8um and sample's name
#Output: Statement of processed object completion
#Description:
banksy.normalization <- function(spatial_object, sample_name) {
  #Filter spatial object for sample
  cat(paste("Filtering main object for:", sample_name,"Running Banksy", "\n",sep = " "))
  obj <- subset(spatial_object,idents = sample_name)
  #Add mitochondrial percentage
  obj[['percent.mt']] <- PercentageFeatureSet(obj, pattern = '^MT-')
  #Normalize data and find variable genes
  obj %<>% NormalizeData(.) %>%
    FindVariableFeatures(., verbose = TRUE) %>% 
    ScaleData(.)
  cat("Applying BANKSY","\n")
  #Run BANKSY
  obj <- SeuratWrappers::RunBanksy(obj,
                                   lambda = 0.2, verbose = TRUE,
                                   assay = "Spatial.008um", slot = "data", 
                                   features = "variable", k_geom = 30)
  #Set BANKSY as default assay
  DefaultAssay(obj) <- "BANKSY"
  
  cat(paste("PCA for",sample_name, "\n",sep = " "))
  obj <- RunPCA(obj, assay = "BANKSY", reduction.name = "pca.banksy", 
                features = rownames(obj), npcs = 30)
  obj <- RunUMAP(obj,dims = 1:15,reduction = "pca.banksy")
  obj <- FindNeighbors(obj, reduction = "pca.banksy", dims = 1:15)
  obj <- FindClusters(obj, resolution = c(0.5,0.8,1))
  Idents(obj) <- "banksy_cluster"
  cat(paste("Normalization complete for",sample_name, "\n",sep = " "))
  return(obj)
}

#Input: Seurat object for a slide with bin size 8um and sample's name
#Output: Statement of processed object completion
#Description:
qc.plotting <- function(spatial_object, sample_name) {
  #PDF output naming
  pdf(file = paste("./Analysis/Plots/",sample_name,"_QC_plots.pdf",sep = ""), 
      width=8.5, height=11,paper = "USr")
  
  #Violin plots of UMIs, genes, and mitochondrial percentage and PCA elbow plot
  qc.plots <- wrap_plots(VlnPlot(spatial_object, 
                                 features = "nCount_Spatial.008um", pt.size = 0),
                         VlnPlot(spatial_object, 
                                 features = "nFeature_Spatial.008um", pt.size = 0),
                         VlnPlot(spatial_object,
                                 features = "percent.mt", pt.size = 0),
                         ElbowPlot(spatial_object, reduction = "pca.banksy"),
                         ncol = 2, guides = "collect")
  print(qc.plots)
  #Resolution UMAP plots
  res.plots <- wrap_plots(DimPlot(spatial_object,
                                  group.by = "BANKSY_snn_res.0.5", reduction = "umap") + 
                           theme(legend.position = "bottom"),
             DimPlot(spatial_object,
                     group.by = "BANKSY_snn_res.0.8", reduction = "umap") + 
               theme(legend.position = "bottom"),
             DimPlot(spatial_object,
                     group.by = "BANKSY_snn_res.1", reduction = "umap") + 
               theme(legend.position = "bottom")
             )
  print(res.plots)
  #Spatial distribution of cluster plots 
  dim.plots <- wrap_plots(
    DimPlot(spatial_object, pt.size = 0.25, label = TRUE, label.size = 3, 
            repel = TRUE, group.by = "BANKSY_snn_res.0.5"),
    SpatialDimPlot(spatial_object, stroke = NA, label = TRUE, label.size = 3, 
                   repel = TRUE, alpha = 0.5, pt.size.factor = 2, 
                   group.by = "BANKSY_snn_res.0.5"), ncol = 2)
  print(dim.plots)  
  dev.off()
  cat("Plots have been saved to PDF \n")
}


#-------------------------------Script Start-----------------------------------
#Load packages and CSVs
packages <- c("R.utils","hdf5r","arrow",
              "future","Seurat","SeuratData","SeuratWrappers","Banksy",
              "tidyverse","patchwork")

install_missing_packages(packages)

lapply(packages,library, character.only=TRUE)

csv <- read.csv(file = "./Visium_Mouse_Tumor_MycCap_Slide_Samples.csv")
slides <- unique(csv$Slide_Number)

#Data Filtering and Processing####
lapply(slides, function(slide) {
  sample.df <- csv %>% filter(Slide_Number == slide)
  #Read in CSV file including the spots annotated with their tissue of origin
  barcodes <- read.csv(file = paste("./Loupe_Spot_CSV/", slide,
                                    "_Treatment_Spot_Annotation.csv",sep = ""))
  
  #Directory contains read count matrix and image data in a sub directory `spatial`.
  slide.obj <- Load10X_Spatial(data.dir = paste("./BANOSSM_SSM0015_1_PR_Whole_C1_VISHD_",
                                        slide,"_22WJCYLT3/outs/",sep = ""), 
                               bin.size=8)
  
  #Filter object for spots within Loupe browser annotations
  filt.obj <- subset(slide.obj,cells = barcodes$Barcode)
  
  #Order spot barcodes by the order in the Seurat object 
  ord.barcodes <- barcodes %>%
    filter(Barcode %in% rownames(filt.obj@meta.data)) %>%
    arrange(match(Barcode,rownames(filt.obj@meta.data)))
  
  #Add annotation data to Seurat object metadata
  filt.obj[["Treatment"]] <- ord.barcodes
  Idents(filt.obj) <- "Treatment"
  
  #Normalize by BANKSY method and store Seurat object
  sapply(sample.df$Sample, function(sample_name) {
    normalized.obj <- banksy.normalization(filt.obj, sample_name)
    cat(paste("Saving RDS","\n",sep = " "))
    saveRDS(normalized.obj, file = paste("./Analysis/BANKSY_Normalized/Seurat_Objects",
                              slide,"_",sample_name,".rds",sep = ""))
    #QC plotting
    qc.plotting(spatial_object = normalized.obj, sample_name = sample_name)
    
    rm(normalized.obj)
    
    message <-paste("Saved seurat objects and plots for",sample_name,"\n",sep = "")
    return(print(message))
  })
  rm(filt.obj)
  return(cat(paste("Saved seurat objects and plots for",slide,"\n",sep = "")))
})




