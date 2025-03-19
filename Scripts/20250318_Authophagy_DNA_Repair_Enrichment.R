#Author: Rogelio Aguilar
#Project: Spatial Dietary Project
#Objective: Visualize autophagy related genes and enrichment for DNA repair pathways


#Set working directory
setwd("/Volumes/workl/Rogelio/")

#Load packages and CSVs
packages <- c("R.utils","hdf5r","arrow",
              "future","Seurat","SeuratData","SeuratWrappers","Banksy",
              "tidyverse","patchwork","msigdbr","GSEABase","GSVA","pheatmap")

lapply(packages,library, character.only=TRUE)

#Functions####
#Input: Seurat object for a slide with bin size 8um and sample's name
#Output: Statement of processed object completion
#Description:
gene.plotting <- function(spatial_object, slide_number, sample_name, gene_set) {
  #Plot spatial gene expression of autophagy genes and add to PDF
  autophagy.plot <- SpatialFeaturePlot(spatial_object, features = gene_set) +
    ggtitle(paste(slide_number,sample_name,sep = "_"))
  print(autophagy.plot)
  #Plot UMAP gene expression of autophagy genes and add to PDF
  autophagy.umap.plot <- FeaturePlot(spatial_object, features = gene_set, label = TRUE) +
    ggtitle(paste(slide_number,sample_name,sep = "_"))
  print(autophagy.umap.plot)
  #Plot dot gene expression of autophagy genes and add to PDF
  dot.plot <- DotPlot(spatial_object, features = gene_set) +
    ggtitle(paste(slide_number,sample_name,sep = "_"))
  print(dot.plot)
  dev.off()
}



#Read in sample CSV and filter for control group and low glucose group(s)
csv <- read.csv(file = "./Analysis/Visium_Mouse_Tumor_All_Slide_Samples.csv") %>% 
  mutate(SampleID = paste(Slide_Number,Sample,sep = "_")) %>%
  mutate(SampleID=gsub("CR_RT","CRRT",SampleID)) %>%
  filter(Treatment == "RT"|Treatment == "CRRT"|Treatment == "CR_RT")
slides <- unique(csv$Slide_Number)

output.dir <- "./Analysis/BANKSY_Normalized_QC_Filtered/Autophagy_DNA_Repair"

if (!dir.exists(output.dir)) {
  dir.create(output.dir, recursive = TRUE)
}

#Data Filtering and Processing####

#Plot autophagy gene list and save into a PDF
autophagy.genes <- c("Ulk1", "Atg2a", "Lc3", "Becn1")

#PDF output naming
pdf(file = paste("./Analysis/BANKSY_Normalized_QC_Filtered/Autophagy_DNA_Repair/Autophagy_plots.pdf",sep = ""), 
    width=17, height=22,paper = "USr")
lapply(slides, function(slide) {
  sample.df <- csv %>% filter(Slide_Number == slide)
  #Load BANKSY normalized object and overlay UMI and Gene counts on UMAP coordinates
  sapply(sample.df$Sample, function(sample_name) {
    cat(paste("Loading object for ", sample_name,"\n",sep = ""))
    spatial_object <- readRDS(file = paste("./Analysis/BANKSY_Normalized_QC_Filtered/",
                                           slide,"_",sample_name,".rds",sep = ""))
    Idents(spatial_object) <- "BANKSY_snn_res.0.5"
    #Plot autophagy genes
    gene.plotting(spatial_object = spatial_object, 
                  slide_number = slide, 
                  sample_name = sample_name, 
                  gene_set = autophagy.genes)
    #Remove object to prevent memory overage on server
    rm(spatial_object)
    return(cat(slide,sample_name," Added","\n"))
    })
  })
dev.off()
#Generate pseudobulk matrix
log2.counts.matrix <- csv %>% 
  mutate(SampleID = paste(Slide_Number,Sample,sep = "_"))  %>%
  {lapply(.$SampleID, function(sample_name) {
    cat(paste("Loading object for ", sample_name,"\n",sep = ""))
    spatial_object <- readRDS(file = paste("./Analysis/BANKSY_Normalized_QC_Filtered/",
                                           sample_name,".rds",sep = ""))
    Idents(spatial_object) <- "BANKSY_snn_res.0.5"
    #Pseudobulk of Spatial Object
    counts <- AggregateExpression(spatial_object, group.by = "Treatment",
                                  assays = "Spatial.008um") %>% 
      as.data.frame(.) %>% rename(sample_name = Spatial.008um)
    #Transform to log2 CPM
    cpm <- sweep(x = counts,MARGIN = 2,STATS = colSums(counts)/1E6,FUN = "/")
    log2.counts <- log(cpm+1, base = 2)
    #Remove object to prevent memory overage on server
    rm(spatial_object)
    return(log2.counts)
  })} %>% do.call(what = "cbind")
colnames(log2.counts.matrix) <- csv$SampleID

write_csv(log2.counts.matrix,file = "./Analysis/BANKSY_Normalized_QC_Filtered/Autophagy_DNA_Repair/RT_CRRT_pseudobulk_log2cpm.csv")

#Gene sets####
#Autophagy Genes
autophagy.genes = c("Ulk1", "Ulk2", "Atg2a", "Ar2B", "Atg4a", "Atg4b", "Atg4c", 
                    "Atg4d", "Atg5","Becn1", "Atg7", "Lc3", "Gabarap", "Atg9a", 
                    "Atg9b", "Atg13", "Atg14", "Atg16l1", "Atg16l2", "Acbd5")
#DNA Repair Genes
homologous.genes = c("Wrn", "Rad50", "Nbn", "Brca2", "Atrx", "Brca1","Rad51", "Rpa1", "Rad52", "Gen1")
NHEJ.genes = c("Lig4", "Dclre1c", "Prkdc", "Xrcc4", "Xrcc6", "Atm", "Xrcc1", "Lig3", "Parp1")

#Heatmap####
#Calculate Z-score of genes in log2 cpm matrix 
z.score <- scale(log2.counts.matrix) 

#Dataframe of data for column annotation
annotation.col <- csv %>% 
  select(c(Cell_Line,Treatment_Group))
rownames(annotation.col) <- csv$SampleID

#Dataframe of data for row annotation
gene.df <- data.frame(Pathway=c(rep("Homologous DNA Repair", length(homologous.genes)),
                                       rep("NHEJ DNA Repair", length(NHEJ.genes))),
                             Genes = c(homologous.genes, NHEJ.genes))
rownames(gene.df) <- gene.df$Gene

annotation.row <- select(gene.df,Pathway)

#Plot z-score as heat map
pheatmap(mat=z.score[c(homologous.genes,NHEJ.genes,"Ar"),], 
         annotation_col=annotation.col,
         annotation_row = annotation.row,
         cluster_cols = TRUE,cluster_rows = TRUE,
         main = "Z-Score of DNA Repair Associated Genes")


#Enrichment####
lapply(list(autophagy.genes,homologous.genes,NHEJ.genes), function(gene.set) {
  df <- data.frame(gs_name = as.character(substitute(gene.set)),
                   gene_symbol = I(list(gene.set)))
  df.genes.list <- split(x = df$gene_symbol, 
                         f = df$gs_name)
  gsva.scores <-gsva(param = ssgseaParam(exprData=as.matrix(log2.counts.matrix),
                                         geneSets = df.genes.list))
  
})
interest.gene.df <- data.frame(
  Autophagy = c("Ulk1", "Ulk2", "Atg2a", "Ar2B", "Atg4a", "Atg4b", "Atg4c", 
                "Atg4d", "Atg5","Becn1", "Atg7", "Lc3", "Gabarap", "Atg9a", 
                "Atg9b", "Atg13", "Atg14", "Atg16l1", "Atg16l2", "Acbd5"),
  Homologous_DNA = c("Wrn", "Rad50", "Nbn", "Brca2", "Atrx", "Brca1", 
                     "Rad51", "Rpa1", "Rad52", "Gen1"),
  NHEJ_DNA = c("Lig4", "Dclre1c", "Prkdc", "Xrcc4", "Xrcc6", "Atm", "Xrcc1", "Lig3", "Parp1")
)
df <- data.frame(gs_name = c("Hernandez"),
                 gene_symbol = I(list("PLK3", "SPATA6", "NFIA", "TAF13", "GSTM4",
                                      "KCTD3", "MEIS1", "TMEM87B", "GBE1", "STAG1",
                                      "SCOC", "ICE1", "RAI14", "GDNF", "P4HA2",
                                      "PDLIM4", "B4GALT7", "TSPAN13", "ZNHIT1", "SMO",
                                      "CNTLN", "CHMP5","FAM214B","USP6NL","TRDMT1",
                                      "ASCC1","TOLLIP","CCND1","RHNO1","C2CD5","ARID2",
                                      "DGKA","PDS5B", "UFM1","BCL2L2","SUSD6","KLC1",
                                      "ADPGK","CREBBP","NOL3", "ACADVL", "EFNB3",
                                      "SLC16A3", "ZBTB7A","DDA1","ARHGAP35","ZC3H4",
                                      "PCIF1","POFUT2", "PATZ1","DYNLT3","SPIN4",
                                      "PLXNA3","SLC10A3","MT-CYB")))

df.genes.list <- split(x = df$gene_symbol, 
                       f = df$gs_name)

gsva.scores <-gsva(param = ssgseaParam(exprData=as.matrix(log2.counts),
                                       geneSets = df.genes.list))

















