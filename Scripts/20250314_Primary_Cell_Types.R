#Author: Rogelio Aguilar
#Project: Spatial Dietary Project
#Objective: Review spot filtering impact on clusters

#Set working directory
setwd("/Volumes/workl/Rogelio/Analysis/")

#Load packages and CSVs
packages <- c("R.utils","hdf5r","arrow",
              "future","Seurat","SeuratData","SeuratWrappers","Banksy",
              "tidyverse","patchwork")

lapply(packages,library, character.only=TRUE)

csv <- read.csv(file = "./Visium_Mouse_Tumor_MycCap_Slide_Samples.csv")
slides <- unique(csv$Slide_Number)

obj <- readRDS(file = "./BANKSY_Normalized/Seurat_Objects/F07833_5_RT.rds")

Idents(obj) <- "BANKSY_snn_res.0.5"

DimPlot(obj)
sp
gene.list <- c("Ar", #androgen receptor
               "Pecam1", #endothelial 
               "Ptprc", #immune
               "Ccsp")

FeaturePlot(obj, features = gene.list, label = TRUE)

FeaturePlot(obj, features = "Pecam1")
