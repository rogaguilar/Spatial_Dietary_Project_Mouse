#Author: Rogelio Aguilar
#Project: Spatial Dietary Project
#Objective: Find variably expressed genes per cluster


#Set working directory
setwd("/home/rogaguilar/Spatial/data/Rose_Li_VisiumHD")

#Load packages and CSVs
packages <- c("R.utils","hdf5r","arrow",
              "future","Seurat","SeuratData","SeuratWrappers","Banksy",
              "tidyverse","patchwork")

lapply(packages,library, character.only=TRUE)

csv <- read.csv(file = "./Visium_Mouse_Tumor_MycCap_Slide_Samples.csv")
slides <- unique(csv$Slide_Number)

#Data Filtering and Processing####
lapply(slides, function(slide) {
  sample.df <- csv %>% filter(Slide_Number == slide)
  #Load BANKSY normalized object and find all markers for each cluster
  sapply(sample.df$Sample, function(sample_name) {
    cat(paste("Loading object for ", sample_name,"\n",sep = ""))
    normalized.obj <- readRDS(file = paste("./Analysis/BANKSY_Normalized/",
                                           slide,"_",sample_name,".rds",sep = ""))
    Idents(normalized.obj) <- "BANKSY_snn_res.0.5"
    #Find markers 
    markers <- FindAllMarkers(object = normalized.obj, assay = "BANKSY",
                              only.pos = TRUE, min.pct = 0.25, 
                              logfc.threshold = 0.25)
    #Save as CSV file
    write_csv(x = markers, file = paste("./Analysis/BANKSY_Normalized/",slide,"_",
                                        sample_name,"_all_markers.csv",sep = ""))
    #Remove object and markers data frame to prevent memory overage on server
    rm(normalized.obj,markers)
    
    message <- paste("Saved marker CSV for ",sample_name,"\n",sep = "")
    return(print(message))
  })
  return(cat(paste("Saved marker CSVs for ",slide,"\n",sep = "")))
})