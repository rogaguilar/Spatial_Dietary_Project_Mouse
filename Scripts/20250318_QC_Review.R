
csv %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)%>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  View(.)

csv <- read.csv(file = "./Analysis/Visium_Mouse_Tumor_MycCap_Slide_Samples.csv")
slides <- unique(csv$Slide_Number)

#Data Filtering and Processing####
lapply(slides[1], function(slide) {
  sample.df <- csv %>% filter(Slide_Number == slide)
  #Read in CSV file including the spots annotated with their tissue of origin
  barcodes <- read.csv(file = paste("./BANOSSM_SSM0015_1_PR_Whole_C1_VISHD_",
                                    slide,"_22WJCYLT3/",slide,
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
    saveRDS(normalized.obj, file = paste("./Analysis/BANKSY_Normalized_10MT_50UMI_50Gene/",
                                         slide,"_",sample_name,".rds",sep = ""))
    #QC plotting
    # qc.plotting(spatial_object = normalized.obj, 
    #             sample_name = sample_name,
    #             slide_number = slide)
    
    rm(normalized.obj)
    
    message <-paste("Saved seurat objects and plots for",sample_name,"\n",sep = "")
    return(print(message))
  })
  rm(filt.obj)
  return(cat(paste("Saved seurat objects and plots for",slide,"\n",sep = "")))
})

#Find Markers####
lapply(slides, function(slide) {
  sample.df <- csv %>% filter(Slide_Number == slide)
  #Load BANKSY normalized object and find all markers for each cluster
  sapply(sample.df$Sample, function(sample_name) {
    cat(paste("Loading object for ", sample_name,"\n",sep = ""))
    normalized.obj <- readRDS(file = paste("./Analysis/BANKSY_Normalized_10MT_50UMI_50Gene/",
                                           slide,"_",sample_name,".rds",sep = ""))
    Idents(normalized.obj) <- "BANKSY_snn_res.0.5"
    #Find markers 
    markers <- FindAllMarkers(object = normalized.obj, assay = "BANKSY",
                              only.pos = TRUE, min.pct = 0.25, 
                              logfc.threshold = 0.25)
    # Specify the directory path
    dir_path <- "./Analysis/BANKSY_Normalized_10MT_50UMI_50Gene/Cluster_Markers"
    
    # Check if the directory exists
    if (!dir.exists(dir_path)) {
      # Create the directory if it doesn't exist
      dir.create(dir_path)
      cat("Directory created:", dir_path, "\n")
    } else {
      cat("Directory already exists:", dir_path, "\n")
    }
    #Save as CSV file
    write_csv(x = markers, file = paste(dir_path,slide,"_",
                                        sample_name,"_all_markers.csv",sep = ""))
    #Remove object and markers data frame to prevent memory overage on server
    rm(normalized.obj,markers)
    
    message <- paste("Saved marker CSV for ",sample_name,"\n",sep = "")
    return(print(message))
  })
  return(cat(paste("Saved marker CSVs for ",slide,"\n",sep = "")))
})


obj <- readRDS(file = "./Analysis/BANKSY_Normalized_10MT_50UMI_50Gene/F07833_5_RT.rds")
Idents(spatial_object) <- "BANKSY_snn_res.0.5"


(FeaturePlot(spatial_object, features = "Ar",label = TRUE)+
  (FeaturePlot(spatial_object, features = "Cd3e",label = TRUE))

FeaturePlot(spatial_object, features = "Ar",label = TRUE)










