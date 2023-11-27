library(tidyverse)
library(scibet)
library(Seurat)
library(scMCA)

seurat_data <- readRDS("scDART_integrated_seurat.rds")
raw_count_data <- as.data.frame(seurat_data$RNA@counts)
cell_to_sample <- seurat_data@meta.data$orig.ident
# table(cell_to_sample)
cell_to_cluster <- seurat_data@active.ident

result <- scMCA(scdata = raw_count_data, numbers_plot = 3)

cluster_labels <- data.frame(cell = names(cell_to_cluster), 
                             cluster = cell_to_cluster,
                             pred_cell_type = result$scMCA)

pred_celltype_per_cluster <- cluster_labels %>%
  group_by(cluster, pred_cell_type) %>%
  summarize(cells = n()) %>%
  arrange(cluster, desc(cells))

write.table(pred_celltype_per_cluster, file = "scMCA_preds.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)