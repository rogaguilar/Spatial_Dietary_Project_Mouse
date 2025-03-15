#Author: Rogelio Aguilar
#Project: Spatial Dietary Project
#Objective: Review UMI mapping in clusters

#Set working directory
setwd("/home/rogaguilar/Spatial/data/Rose_Li_VisiumHD")

#Load packages and CSVs
packages <- c("R.utils","hdf5r","arrow",
              "future","Seurat","SeuratData","SeuratWrappers","Banksy",
              "tidyverse","patchwork","png","grid")

lapply(packages,library, character.only=TRUE)

csv <- read.csv(file = "./Visium_Mouse_Tumor_MycCap_Slide_Samples.csv")
slides <- unique(csv$Slide_Number)

#Data Filtering and Processing####
lapply(slides, function(slide) {
  sample.df <- csv %>% filter(Slide_Number == slide)
  #Load BANKSY normalized object and overlay UMI and Gene counts on UMAP coordinates
  sample.plot.list <- lapply(sample.df$Sample, function(sample_name) {
    cat(paste("Loading object for ", sample_name,"\n",sep = ""))
    spatial_object <- readRDS(file = paste("./Analysis/BANKSY_Normalized_QC_Filtered/",
                                           slide,"_",sample_name,".rds",sep = ""))
    Idents(spatial_object) <- "BANKSY_snn_res.0.5"
    #Violin and UMAP plots of UMIs and genes
    qc.plot <- wrap_plots(FeaturePlot(spatial_object, features = "nCount_Spatial.008um",
                                      label = TRUE,alpha = 0.5, raster = FALSE) + 
                            labs(title = "UMIs"),
                          VlnPlot(spatial_object, group.by = "BANKSY_snn_res.0.5",
                                  features = "nCount_Spatial.008um", pt.size = 0) + 
                            labs(title = "UMIs"),
                          FeaturePlot(spatial_object, features = "nFeature_Spatial.008um",
                                      label = TRUE,alpha = 0.5, raster = FALSE) + 
                            labs(title = "Genes"),
                          VlnPlot(spatial_object, group.by = "BANKSY_snn_res.0.5",
                                  features = "nFeature_Spatial.008um", pt.size = 0) + 
                            labs(title = "Genes"),
                          ncol = 2) + 
      plot_annotation(title = paste(slide,sample_name,sep = "_"))
    ggsave(path = "./Analysis/BANKSY_Normalized_QC_Filtered/Plots/",
           filename = paste(slide,sample_name,"UMI_Genes_UMAP.png",sep = "_"),
           plot = qc.plot, device = "png", width = 24, height = 12, units = "in")
    #Remove object to prevent memory overage on server
    rm(spatial_object, qc.plot)
    message <- paste(slide,sample_name,"UMI/Genes UMAP Saved",sep = " ")
    return(cat(message))
  })
  return(paste(slide,"UMI/Genes UMAP Saved",sep = " "))
})

plot.path.list <-list.files(path = "./Analysis/BANKSY_Normalized_QC_Filtered/Plots",
                            pattern = "UMI_Genes_UMAP.png",full.names = TRUE)

#Plot into PDF
pdf(file = "./Analysis/BANKSY_Normalized_QC_Filtered/Plots/All_Slides_QC_UMI_Gene_plots.pdf", 
    width=24, height=12,paper = "USr")
for (plot.path in plot.path.list) {
  #Read in image 
  img <- readPNG(plot.path)
  
  # Start a new page for each image
  grid.newpage()
  # Plot the image on a new page in the PDF
  grid.raster(img)
}
dev.off()