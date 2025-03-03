#Author: Rogelio Aguilar
#Project: Spatial Dietary Project
#Objective: To BANKSY normalized Visium spatial objects, visualize gene, feature,
# and dimensional reduction plots.


#Set working directory
setwd("/home/rogaguilar/Spatial/data/Rose_Li_VisiumHD")

#Functions####
#Input: Seurat object for a slide with bin size 8um and sample's name
#Output: Statement of processed object completion
#Description:
qc.plotting <- function(spatial_object, slide_number, sample_name) {
  #PDF output naming
  pdf(file = paste("./Analysis/Plots/",slide_number,"_",sample_name,"_QC_plots.pdf",sep = ""), 
      width=17, height=22,paper = "USr")
  
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
                                  group.by = "BANKSY_snn_res.0.5", reduction = "umap",
                                  label = TRUE) + 
                            theme(legend.position = "bottom"),
                          DimPlot(spatial_object,
                                  group.by = "BANKSY_snn_res.0.8", reduction = "umap",
                                  label = TRUE) + 
                            theme(legend.position = "bottom"),
                          DimPlot(spatial_object,
                                  group.by = "BANKSY_snn_res.1", reduction = "umap",
                                  label = TRUE) + 
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

lapply(packages,library, character.only=TRUE)

csv <- read.csv(file = "./Visium_Mouse_Tumor_MycCap_Slide_Samples.csv")
slides <- unique(csv$Slide_Number)

#Data Plotting####
lapply(slides, function(slide) {
  sample.df <- csv %>% filter(Slide_Number == slide)
  #Normalize by BANKSY method and store Seurat object
  sapply(sample.df$Sample, function(sample_name) {
    normalized.obj <- readRDS(file = paste("./Analysis/BANKSY_Normalized/",
                                           slide,"_",sample_name,".rds",sep = ""))
    #QC plotting
    qc.plotting(spatial_object = normalized.obj, slide_number = slide, sample_name = sample_name)
    
    #Remove object from local memory to prevent memory overage
    rm(normalized.obj)
    
    message <-paste("Saved seurat objects and plots for",sample_name,"\n",sep = "")
    return(print(message))
  })
  return(cat(paste("Saved seurat objects and plots for",slide,"\n",sep = "")))
})




