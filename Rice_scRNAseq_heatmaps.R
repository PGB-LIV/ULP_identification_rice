# Method: To gather the expression of ULP genes at single cell resolution, 
# RDS files were retrieved from NCBI GEO GSE232863(Wang et. al., 2025) and GSE251706(Zhu et. al., 2025). 
# Seurat(v5.2.1) and R version 4.4.0 were used to retrieve the expression levels of ULP genes(MSU/RAPDB annotations) from the rds object and plot as a heatmap using ComplexHeatmap package. 

# Wang X, Huang H, Jiang S, Kang J et al. A single-cell multi-omics atlas of rice. Nature 2025 Aug;644(8077):722-730. PMID: 40634611
# Zhu M, Hsu CW, Peralta Ogorek LL, Taylor IW et al. Single-cell transcriptomics reveal how root tissues adapt to soil stress. Nature 2025 Jun;642(8068):721-729. PMID: 40307555

setwd("./sumo_scrna")


suppressMessages(library(readxl))
suppressMessages(library(svglite))
suppressMessages(library(Seurat)) ### v5.2.1
suppressMessages(library(tidyverse))
suppressMessages(library(edgeR))
suppressMessages(library(ComplexHeatmap))




# Read ULP ids , available as a table in the article
ULPs_Osativa_Japonica_3models_list <- read_excel("ULPs_Osativa_Japonica_3models_list.xlsx")

# get rapdb and MSU ids
rgap_ids_ulp <- ULPs_Osativa_Japonica_3models_list %>% select(`RGAP (MSU)`) %>% filter(!(`RGAP (MSU)` == "None" )) # 42 , used in Zhu et al 2025 (GSE251706)

rapdb_ids_ulp <- ULPs_Osativa_Japonica_3models_list %>% select(`RAP-db`) %>% filter(!(`RAP-db` == "None" )) # 36 , used in Wang et al 2025 (GSE232863)



# 1. Analysis of scRNAseq data from Wang et al 2025 (GSE232863)

#Read RDS
wang_2025_paper_rds <- readRDS("GSE232863.Rds")

# Get RAP-db gene IDs as a vector to match the seurat object data
rapdb_genes <- rapdb_ids_ulp$`RAP-db`


# Aggregate (pseudobulk) expression values and retrieve from seurat object across cell type annotations

wang_pseudo_bulk_paper_scd.seu <- AggregateExpression(wang_2025_paper_rds, return.seurat = TRUE, assays = "RNA", slot = "data" , features = rapdb_genes , group.by = "tissue_cluster_names" ,  normalization.method = "LogNormalize")

log_norm_pseudo_bulk <- GetAssayData(wang_pseudo_bulk_paper_scd.seu, assay = "RNA", layer = "data")[rownames(wang_pseudo_bulk_paper_scd.seu) %in% rapdb_genes, , drop = FALSE] 


# For plotting as Heatmap - add Gene names
log_norm_pseudo_bulk <- log_norm_pseudo_bulk %>% as.data.frame() %>% 
  rownames_to_column(., "rowname") %>% 
  left_join(ULPs_Osativa_Japonica_3models_list[,c(2,5)] , join_by(rowname == `RAP-db`))

log_norm_pseudo_bulk <- log_norm_pseudo_bulk %>% dplyr::rename(gene_name = `Gene name  ` ) %>%
  distinct() %>% 
  column_to_rownames(., "gene_name") %>%
  select(-rowname)


# Convert into a matrix and plot heatmap 
log_norm_pseudo_bulk_matrix <- as.matrix(log_norm_pseudo_bulk)


# Uncomment to export
# svglite("plots/wang_heatmap_log_norm_pseudo_bulk.svg" , width = 11.5,height = 8.5 , bg = "white",pointsize = 12)

#breaks <- seq(min(summed_counts_raw_log_matrix), (max(summed_counts_raw_log_matrix))/4, length.out=100)

# expression level is log normalised and  scaled values as per AggregateExpression function
ComplexHeatmap::pheatmap(log_norm_pseudo_bulk_matrix,  na_col = "grey", cellheight = 8, cellwidth = 10,fontsize_row = 8, cluster_rows = TRUE,
                         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete",
                         cluster_cols = TRUE , heatmap_legend_param = list(
                           title = "Expression level",
                           legend_height = unit(4, "cm"),
                           title_position = "topcenter" , 
                           title_gp = gpar(fontsize = 10, fontface = "bold"),
                           labels_gp = gpar(fontsize = 8),
                           legend_direction = "horizontal",
                           legend_position = "top"))

#dev.off()



# 2. Analysis of scRNAseq data from Zhu et al 2025 (GSE251706 , WT and gel growm samples)

#Read RDS
zhu_2025_paper_rds <- readRDS("GSE251706_WT_atlas_seu3.rds")

# Update to seurat v5
zhu_2025_paper_rds <- UpdateSeuratObject(zhu_2025_paper_rds)

# Gather MSU ids and replace "_" with "-" as per the object data
msu_ids <- unique(rgap_ids_ulp$`RGAP (MSU)`)
msu_ids <- str_replace_all(msu_ids, "_", "-")


# Aggregate (pseudobulk) expression values and retrieve from seurat object across cell type annotations
zhu_pseudo_bulk_paper.seu <- AggregateExpression(zhu_2025_paper_rds, return.seurat = TRUE, assays = "RNA", features = msu_ids , group.by = "celltype.anno" ,  normalization.method = "LogNormalize")

zhu_log_norm_pseudo_bulk <- GetAssayData(zhu_pseudo_bulk_paper.seu, assay = "RNA", layer = "data")[rownames(zhu_pseudo_bulk_paper.seu) %in% msu_ids, , drop = FALSE] 


# For plotting as Heatmap- add Gene names
zhu_log_norm_pseudo_bulk <- zhu_log_norm_pseudo_bulk %>% as.data.frame() %>% 
  rownames_to_column(., "rowname") %>% 
  mutate(rowname = str_replace_all(rowname, "-", "_")) %>% 
  left_join(ULPs_Osativa_Japonica_3models_list[,c(4,5)] , join_by(rowname == `RGAP (MSU)`))


zhu_log_norm_pseudo_bulk <- zhu_log_norm_pseudo_bulk %>% dplyr::rename(gene_name = `Gene name  ` ) %>%
  distinct() %>% 
  column_to_rownames(., "gene_name") %>%
  select(-rowname)



# Convert into a matrix and plot heatmap 
zhu_log_norm_pseudo_bulk_matrix <- as.matrix(zhu_log_norm_pseudo_bulk)

# Uncomment to export
# svglite("plots/zhu_heatmap_log_norm_pseudo_bulk.svg" , width = 8.5,height = 11.5 , bg = "white",pointsize = 12)

breaks <- seq(min(zhu_log_norm_pseudo_bulk_matrix), (max(zhu_log_norm_pseudo_bulk_matrix)), length.out=100)

# expression level is log normalised and scaled values as per AggregateExpression function
ComplexHeatmap::pheatmap(zhu_log_norm_pseudo_bulk_matrix,  na_col = "grey", cellheight = 8, cellwidth = 10,fontsize_row = 8, cluster_rows = TRUE,
                         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete",
                         cluster_cols = TRUE ,  breaks = breaks, heatmap_legend_param = list(
                           title = "Expression level",
                           legend_height = unit(4, "cm"),
                           title_position = "topcenter" , 
                           title_gp = gpar(fontsize = 10, fontface = "bold"),
                           labels_gp = gpar(fontsize = 8),
                           legend_direction = "horizontal",
                           legend_position = "top"))

# dev.off()

