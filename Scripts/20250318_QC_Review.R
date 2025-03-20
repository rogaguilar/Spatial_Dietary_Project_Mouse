#Author: Rogelio Aguilar
#Project: Spatial Dietary Project
#Objective: Perform pipeline using >50UMIs,>50genes,and <5%MT


#Set working directory
setwd("/Volumes/workl/Rogelio/")


csv <- read.csv(file = "./Analysis/Visium_Mouse_Tumor_MycCap_Slide_Samples.csv")
slides <- unique(csv$Slide_Number)

#Functions####
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
    pdf(file = paste("./Analysis/BANKSY_Normalized_10MT_50UMI_50Gene/Plots/Enrichment/",slide_number,"_",
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
                  prefix = paste("./Analysis/BANKSY_Normalized_10MT_50UMI_50Gene/Cluster_Markers/Enrichment/",
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
lapply(slides[1], function(slide) {
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

#Enrichment####
output.dir <- "./Analysis/BANKSY_Normalized_10MT_50UMI_50Gene/Plots/Enrichment"

if (!dir.exists(output.dir)) {
  dir.create(output.dir, recursive = TRUE)
}

lapply(slides[1], function(slide) {
  sample.df <- csv %>% filter(Slide_Number == slide)
  #Enrichment of sample
  lapply(sample.df$Sample, function(sample_name) {
    #Load in marker csv file
    markers.csv <- read.csv(file = paste("./Analysis/BANKSY_Normalized_10MT_50UMI_50Gene/Cluster_Markers/",
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


spatial_object <- readRDS(file = "./Analysis/BANKSY_Normalized_10MT_50UMI_50Gene/F07833_5_RT.rds")
Idents(spatial_object) <- "BANKSY_snn_res.0.5"




FeaturePlot(spatial_object, features = "Ar",label = TRUE) + 
              DotPlot(spatial_object, features = c("Cd3d","Cd3e","Cd3g","Cd4","Cd8a"),assay = "Spatial.008um")

DotPlot(spatial_object, features = c("Cd19"),assay = "Spatial.008um")+ 
  DotPlot(spatial_object, features = c("Cd3d","Cd3e","Cd3g","Cd4","Cd8a"),assay = "Spatial.008um")

DotPlot(spatial_object, features = c("Siglec1"),assay = "Spatial.008um")
SpatialFeaturePlot(spatial_object, features = "Psca",image.alpha = 0.02)

SpatialDimPlot(spatial_object, group.by = "BANKSY_snn_res.0.5")
  

DimPlot(spatial_object,group.by = "BANKSY_snn_res.0.5", label = TRUE)
SpatialDimPlot(spatial_object,group.by = "BANKSY_snn_res.0.5")

obj <- readRDS(file = "./Analysis/BANKSY_Normalized_QC_Filtered_minUMI_25/F07833_5_RT.rds")
SpatialDimPlot(obj,group.by = "BANKSY_snn_res.0.5")


marker.csv <- read.csv(file = "D:/Rogelio/Analysis/BANKSY_Normalized_10MT_50UMI_50Gene/Cluster_Markers/F07833_5_RT_all_markers.csv")

top10.marker.csv <- marker.csv %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)%>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  mutate(gene = gsub(pattern = "\\.m0$",replacement = "", x = `gene`))

write.csv(top10.marker.csv,file = "D:/Rogelio/Analysis/BANKSY_Normalized_10MT_50UMI_50Gene/Cluster_Markers/F07833_5_RT_top10_markers.csv")



ImageFeaturePlot(obj, features = c("Ar"))

SpatialFeaturePlot(obj, features = "Ar",) + SpatialFeaturePlot(obj,features = "Psca")

wrap_plots(SpatialDimPlot(obj,group.by = "BANKSY_snn_res.0.5"),
           DimPlot(obj,group.by = "BANKSY_snn_res.0.5", label = TRUE),
           SpatialFeaturePlot(obj, features = "Ar"),
           SpatialFeaturePlot(obj, features = "Psca"))


SpatialFeaturePlot(obj, features = c("Cd8a","Ptprc","Cd3g"))

SpatialFeaturePlot(obj, features = c("nCount_Spatial.008um","nFeature_Spatial.008um"))
#PCa markers
SpatialFeaturePlot(obj, features = c("Ar","Psca","Myc"))


ImageFeaturePlot(obj, features = c("Ar"))

SpatialFeaturePlot(obj, features = "Ar",pt.size.factor = 5) + SpatialFeaturePlot(obj,features = "Psca",pt.size.factor = 5)

wrap_plots(SpatialDimPlot(obj,group.by = "BANKSY_snn_res.0.5"),
           DimPlot(obj,group.by = "BANKSY_snn_res.0.5", label = TRUE),
           SpatialFeaturePlot(obj, features = "Ar"),
           SpatialFeaturePlot(obj, features = "Psca"))


SpatialFeaturePlot(obj, features = c("Cd3d","Cd3e","Cd3g"))











