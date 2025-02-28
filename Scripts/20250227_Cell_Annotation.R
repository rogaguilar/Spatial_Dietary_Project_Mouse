#Load packages and CSVs
packages <- c("tidyverse","patchwork", "enrichR","org.Mm.eg.db")

lapply(packages,library, character.only=TRUE)

#Load in Samples markers csv file
markers <- read_csv(file = "/Volumes/workl/Rogelio/Data_Analysis/Marker_CSV/F07833_5_RT_all_markers.csv") %>%
  mutate(gene = gsub(pattern = "\\.m0$",replacement = "", x = `gene`))

#Filter csv for genes with log fold change over 1
top10 <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() %>% filter(cluster==1)

#Select cluster of interest 
clust0.markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  ungroup() %>% 
  filter(cluster==0)

# Convert Gene Symbols to Entrez IDs for KEGG analysis
entrez_ids <- mapIds(org.Mm.eg.db, 
                     keys = clust0.markers$gene, 
                     column = "ENTREZID", 
                     keytype = "SYMBOL", 
                     multiVals = "first")

#Apply enrichment through accessing enrichR API
enriched <- enrichR::enrichr(genes = clust0.markers$gene,
                 databases = c("Azimuth_Cell_Types_2021", "Tabula_Muris"))

#Store enrichment into Excel file; each database profiled will appear as a tab 
printEnrich(data = enriched, outFile = "excel", 
            prefix = "test") #name Excel file will be stored as







