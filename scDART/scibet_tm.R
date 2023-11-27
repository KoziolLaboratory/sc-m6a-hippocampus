library(tidyverse)
library(scibet)

# TPM divides counts by transcript length but with UMI-tagged data, it isn't necessarily true that longer length -> more counts. I'd recommend not dividing by transcript length for 10X data.
# Just don't divide by transcript lengths. Just take a gene's UMI count and divide it by the total number of UMIs in a cell (this is essentially what TPM is except we're not dividing by transcript length).
library(Seurat)
seurat_data <- readRDS("scDART_integrated_seurat.rds")
raw_count_data <- t(as.data.frame(seurat_data$RNA@counts))
cell_to_sample <- seurat_data@meta.data$orig.ident
# table(cell_to_sample)
cell_to_cluster <- seurat_data@active.ident

model <- readr::read_csv("GSE109774_scibet_core.csv")
model <- pro.core(model)
prd <- LoadModel(model)
length(intersect(colnames(raw_count_data), rownames(model))) #831/1000. 

# Model1: raw counts
pred_label1 <- prd(raw_count_data)

cluster_labels1 <- data.frame(cell = names(cell_to_cluster), 
                              cluster = cell_to_cluster,
                              pred_cell_type = pred_label1)

pred_celltype_per_cluster1 <- cluster_labels1 %>%
  group_by(cluster, pred_cell_type) %>%
  summarize(cells = n()) %>%
  arrange(cluster, desc(cells))

write.table(pred_celltype_per_cluster1, file = "scibet_preds.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

