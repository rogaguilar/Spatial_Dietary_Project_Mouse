#Author: Rogelio Aguilar
#Project: Spatial Dietary Project
#Objective: Review spot filtering impact on clusters

#Set working directory
setwd("D:/Rogelio/Analysis/")

#Load packages and CSVs
packages <- c("R.utils","hdf5r","arrow",
              "future","Seurat","SeuratData","SeuratWrappers","Banksy",
              "tidyverse","patchwork","readxl")

lapply(packages,library, character.only=TRUE)

#Read in CSV file with sample information
csv <- read.csv(file = "./Visium_Mouse_Tumor_MycCap_Slide_Samples.csv")
slides <- unique(csv$Slide_Number)

####Functions####
# Function to read, filter, and assign sample number based on file name
process_file <- function(file_path, sample_name) {
  # Read the sheet "Tabula_Muris"
  data <- tryCatch({
    read_xlsx(file_path, sheet = "Tabula_Muris")
  }, error = function(e) {
    message("Error reading file: ", file_path)
    return(NULL)
  })
  
  if (is.null(data)) return(NULL)
  
  # Ensure consistent column types (convert "Term", "Overlap", and "Genes" columns to character)
  if("Term" %in% colnames(data)) {
    data$Term <- as.character(data$Term)
  }
  if("Overlap" %in% colnames(data)) {
    data$Overlap <- as.character(data$Overlap)
    data$Overlap <- paste0("'", data$Overlap)  # Add single quote to prevent date conversion
  }
  if("Genes" %in% colnames(data)) {
    data$Genes <- as.character(data$Genes)
  }
  
  # Extract the cluster number from the filename
  cluster_number <- str_extract(basename(file_path), "(?<=cluster_)\\d+(?=_enrichment)")
  
  # Filter rows where p-value is less than 0.05
  filter.df <- data %>% filter(P.value < 0.05)
  
  # If the filtered dataframe is empty, add an NA row
  if (nrow(filter.df) == 0) {
    filter.df <- data.frame(Term = NA, Overlap = NA, Genes = NA, 
                            P.value = NA, sample_name = sample_name, 
                            cluster = cluster_number)
  } else {
    # Add sample_name and cluster columns
    filter.df <- filter.df %>%
      mutate(sample_name = sample_name,
             cluster = cluster_number)
  }
  
  return(filter.df)
}

#Data Filtering and Processing####

#Collect cell types enriched within clusters and store in a CSV per sample
lapply(slides, function(slide) {
  sample.df <- csv %>% filter(Slide_Number == slide)
  
  sapply(sample.df$Sample, function(sample_name) {
    # Collect paths for enrichment files 
    folder.path <- "./BANKSY_Normalized_QC_Filtered_minUMI_25/BANKSY_Normalized_QC_Filtered/Cluster_Markers/Enrichment/"
    path.list <- list.files(path = folder.path, 
                            pattern = paste(slide, sample_name, sep = "_"),
                            full.names = TRUE)
    
    # Exclude CSV files
    path.list <- path.list[!grepl("\\.csv$", path.list)]
    cat("Processing: ",slide, sample_name, "\n")
    # Process all files for this sample (all cluster files for this sample)
    combined.data <- path.list %>%
      map_df(~ process_file(.x, sample_name))
    
    if (!is.null(combined.data)) {
      # Order the combined data numerically by the cluster column
      combined.data <- combined.data %>% arrange(as.numeric(cluster))
      # Write the combined data to a CSV file
      write.csv(combined.data, file = paste0(folder.path, slide, "_", sample_name, ".csv"), 
                row.names = FALSE)
    }
  })
})



#-------------------------Testing--------------------------------------------
obj <- readRDS(file = "./BANKSY_Normalized_QC_Filtered_minUMI_25/BANKSY_Normalized_QC_Filtered/F07833_5_RT.rds")

Idents(obj) <- "BANKSY_snn_res.0.5"

DimPlot(obj)

gene.list <- c("Ar", #androgen receptor
               "Pecam1", #endothelial 
               "Ptprc", #immune
               "KRT5")

FeaturePlot(obj, features = gene.list, label = TRUE)

krt <- c("Cd3e","Cd4","Cd8a","Cb8b1","Foxp3")
FeaturePlot(obj, features = c("Krt8", "Krt18", "Psca", "Krt4", "Tacstd2", "Pigr39"), label = TRUE)

FeaturePlot(obj, features = "Psca", label = TRUE)

DotPlot(obj, features = "Cd3e",assay = "BANKSY",)

#Immune
t.cell <- c("Cd3d","Cd3e","Cd3g","Cd4","Cd8a","Cb8b1","Foxp3")
nk.cell <- c("Ncam1","Fcgr3a")
myc <- c("Ar","Myc","Ern1")
FeaturePlot(obj, features =  c("Cd3e","Cd4","Cd8a","Cb8b","Foxp3"), label = TRUE)

FeaturePlot(obj,features =  c("Cd3d","Cd3e","Cd3g"),label = TRUE,raster = FALSE)

DotPlot(obj, features = "Cd3e",assay = "Spatial.008um")

t.cell.list <- FeaturePlot(obj,features =  c("Cd3d","Cd3e","Cd3g"),
                           label = TRUE,raster = FALSE,combine = FALSE)

wrap_plots(FeaturePlot(obj,features =  c("Cd3d","Cd3e","Cd3g"),
                         label = TRUE,raster = FALSE,combine = FALSE),
           DotPlot(obj, features = c("Cd3d","Cd3e","Cd3g"),assay = "Spatial.008um"))

FeaturePlot(obj,features =  t.cell,
            label = TRUE,raster = FALSE,) +
  DotPlot(obj, features = t.cell,assay = "Spatial.008um")

SpatialDimPlot()


FeaturePlot(obj,features =  "Ar",
            label = TRUE,raster = FALSE) +
  DotPlot(obj, features = t.cell,assay = "Spatial.008um")

#------------------------------------------------------------------------------

lapply(slides, function(slide) {
  sample.df <- csv %>% filter(Slide_Number == slide)
  
  sapply(sample.df$Sample, function(sample_name) {
    normalized.obj <- readRDS(file = paste("./BANKSY_Normalized_QC_Filtered_minUMI_25/BANKSY_Normalized_QC_Filtered/",
                                   slide,"_",sample_name,".rds",sep = ""))
    metadata <- normalized.obj@meta.data %>% 
      rownames_to_column("SpotID") %>% 
      select(c("SpotID","nCount_Spatial.008um"))
    plot <- metadata %>% 
      hist(.$nCount_Spatial.008um, breaks = 100, 
           main = paste(slide,sample_name,"UMI_Count",sep = "_"))
    ggsave(filename = paste(slide,sample_name,"UMI_Count.jpg",sep = "_"),
           plot = plot)
  })
})

metadata %>% 
  mutate(rank = rank(-.$nCount_Spatial.008um,ties.method = "first")) %>%
  ggplot(aes(x = rank, y = "nCount_Spatial.008um")) +
  geom_point() +
  geom_histogram(binwidth = 100) +
  theme_minimal() +
  labs(title = "Rank Barcode Plot", x = "Rank", y = "Value")


metadata %>% 
  mutate(rank = rank(-.$nCount_Spatial.008um, ties.method = "first")) %>%
  ggplot(aes(x = rank, y = nCount_Spatial.008um)) +
  geom_bar(stat = "identity") + theme_bw() +
  labs(title = "Rank Barcode Plot for F07833_5_RT", x = "Rank", y = "UMI Count")

metadata %>% 
  mutate(rank = order(nCount_Spatial.008um, decreasing = FALSE)) %>%
  ggplot(aes(x = rank, y = nCount_Spatial.008um)) +
  geom_bar(stat = "identity")

metadata %>% 
  ggplot(aes(x = nCount_Spatial.008um))  +
  geom_histogram(bins=50) + theme_bw()

metadata %>% 
hist(metadata$nCount_Spatial.008um,breaks = 50,main = "F07833_5_RT UMI Count") +
  hist(metadata$nFeature_Spatial.008um,breaks = 50,,main = "F07833_5_RT Gene Count") 


