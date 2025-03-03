#Author: Rogelio Aguilar
#Project: Spatial Dietary Project
#Objective: Plot cell type enrichment from enrichR databases 


#Set working directory
setwd("/home/rogaguilar/Spatial/data/Rose_Li_VisiumHD")

#Functions####
#Input: Marker CSV, databases to profile as a list, number on slide, sample being profiled
#Output: Statement of plotting completion
#Description:
enrichment.plotting <- function(df, databases, slide_number,sample_name) {
  #Gather number of clusters
  clusters <- unique(df$cluster)
  database.plots <- lapply(databases, function(db) {
    clust.plots <- lapply(clusters,function(clust) {
      cluster.markers <- df %>%
        dplyr::filter(avg_log2FC > 1) %>%
        filter(cluster==clust)
    #Apply enrichment through accessing enrichR API
    enriched <- enrichR::enrichr(genes = cluster.markers$gene, databases = db)
    #Plot enrichment for the top 20 terms ordered by p-value
    plot <- plotEnrich(enriched[[db]], showTerms = 20, numChar = 40, 
                       y = "Count", orderBy = "P.value") + 
      ggtitle(label = paste(db,"for cluster", clust, sep = " "))
    return(plot)})
    #Combine the cluster plots
    combined.clust.plots <- wrap_plots(clust.plots)
    return(combined.clust.plots) })
  #PDF output naming
  pdf(file = paste("./Analysis/Plots/",slide_number,"_",sample_name,"_enrichment_plots.pdf",sep = ""), 
      width=17, height=22,paper = "USr")  
  #Print into the PDF
  lapply(database.plots,print)
  
  message <- paste("Saved plots for",sample_name,"\n",sep = "")
  
  dev.off()
  cat("Plots have been saved to PDF \n")
}

#Input: Marker CSV, databases to profile as a list, number on slide, sample being profiled
#Output: Statement of file storage completion
#Description:
enrichment.excel <- function(df, databases, slide_number,sample_name) {
  #Gather number of clusters
  clusters <- unique(df$cluster)
  lapply(clusters, function(clust) {
    #Run enrichment for selected databases and store excel file
    cluster.markers <- df %>%
      dplyr::filter(avg_log2FC > 1) %>%
      filter(cluster==clust)
    #Apply enrichment through accessing enrichR API
    enriched <- enrichR::enrichr(genes = cluster.markers$gene, databases = databases)
    #Store enrichment into Excel file; each database profiled will appear as a tab 
    printEnrich(data = enriched, outFile = "excel",
                prefix = paste("./Analysis/BANKSY_Normalized/Cluster_Markers/",
                               slide_number,"_", sample_name,"_cluster_",clust,
                               sep = "")) #name Excel file will be stored as
    
    message <- paste("Saved excel for: ", slide_number, "_", sample_name,"_cluster_",clust,"\n",sep = "")
    return(cat(message))
  })
  cat("Excel files have been saved \n")
}

#-------------------------------Script Start-----------------------------------
#Load packages and CSVs
packages <- c("tidyverse","patchwork", "enrichR","org.Mm.eg.db")

lapply(packages,library, character.only=TRUE)

sample.csv <- read.csv(file = "./Visium_Mouse_Tumor_MycCap_Slide_Samples.csv")
slides <- unique(sample.csv$Slide_Number)

#Enrichment####
lapply(slides, function(slide) {
  sample.df <- sample.csv %>% filter(Slide_Number == slide)
  #Enrichment of sample
  lapply(sample.df$Sample, function(sample_name) {
    #Load in marker csv file
    markers.csv <- read.csv(file = paste("./Analysis/BANKSY_Normalized/Cluster_Markers/",
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
                        sample_name = sample_name) })
  return(print("Slide complete"))
})