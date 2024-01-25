# Description:  In this script we are preparing the input for running metaplotR on m/k = 1 variants per cell type. The output is a bed file

library(tidyverse)

variant_filename <- "/GPFS/Magda_lab_temp/maitenat/scs/bullseye/6J-3M-2ws-V2-3-YTH_edit-sites_bc_filtered_all_yth_filt-all.bed"
cell_to_cluster_filename <- "/GPFS/Magda_lab_temp/maitenat/scs/cell_to_cluster.txt"

immune_cells <- as.character(c(0, 1, 4, 8, 9, 10, 14, 15, 17, 18, 19, 20, 21, 24, 25, 27))
oligodendrocytes <- as.character(c(2, 3, 5, 6, 7, 12, 13, 16))
neurons <- as.character(c(11, 23))
endothelial <- "26"
astrocytes <- "22"
cluster_df <- tibble(cell_type = c(rep("IMC", length(immune_cells)), rep("OLG", length(oligodendrocytes)), rep("NEUR", length(neurons)), "EC", "ASC"),
                     cluster_id = c(immune_cells, oligodendrocytes, neurons, endothelial, astrocytes))

variant_data <- read.table(variant_filename, header = FALSE, fill = TRUE)
colnames(variant_data) <- c("chr", "start", "end", "gene", "loc", "change0", "m", "dart_ratio/control_ratio", "mk", "strand", "unk", "unk2", "mk2", "k", "change1", "cell_id", "unk3", "mk_overall", "mk_mut", "mk_ratio")
cell_to_cluster <- read.table(cell_to_cluster_filename, header = TRUE)
variant_data$cluster <- cell_to_cluster$cluster[match(variant_data$cell_id, cell_to_cluster$cell)]
# rm NAs
# los NAs son porque 804 celulas no están asignadas a ningún cluster
variant_data <- variant_data[!is.na(variant_data$cluster),]

cell_type_idx <- match(variant_data$cluster, cluster_df$cluster_id)
variant_data$cell_type <- cluster_df$cell_type[cell_type_idx]

variant_data %>% 
  filter(cell_type == "ASC", mk == 1) %>%
  select(chr, start, end, gene, mk, strand) %>%
write.table("/GPFS/Magda_lab_temp/maitenat/scs/bullseye/6J-3M-2ws-V2-3-YTH_edit-sites_bc_filtered_all_yth_filt-all_mk1_ASC.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

variant_data %>% 
  filter(cell_type == "EC", mk == 1) %>%
  select(chr, start, end, gene, mk, strand) %>%
write.table("/GPFS/Magda_lab_temp/maitenat/scs/bullseye/6J-3M-2ws-V2-3-YTH_edit-sites_bc_filtered_all_yth_filt-all_mk1_EC.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

variant_data %>% 
  filter(cell_type == "IMC", mk == 1) %>%
  select(chr, start, end, gene, mk, strand) %>%
write.table("/GPFS/Magda_lab_temp/maitenat/scs/bullseye/6J-3M-2ws-V2-3-YTH_edit-sites_bc_filtered_all_yth_filt-all_mk1_IMC.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

variant_data %>% 
  filter(cell_type == "NEUR", mk == 1) %>%
  select(chr, start, end, gene, mk, strand) %>%
write.table("/GPFS/Magda_lab_temp/maitenat/scs/bullseye/6J-3M-2ws-V2-3-YTH_edit-sites_bc_filtered_all_yth_filt-all_mk1_NEUR.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

variant_data %>% 
  filter(cell_type == "OLG", mk == 1) %>%
  select(chr, start, end, gene, mk, strand) %>%
write.table("/GPFS/Magda_lab_temp/maitenat/scs/bullseye/6J-3M-2ws-V2-3-YTH_edit-sites_bc_filtered_all_yth_filt-all_mk1_OLG.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## We also want to run metaplotR on all the variants of the cell type, and also on those that have m/k<1, so we need to prepare them. The easiest way for doing so is just take the variants, do the filtering and run metaplotR instead of taking metaplotR results from before. Me da miedo que haya confusión así que lo hago desde este file y punto, aunque no esa lo más eficiente.

#### All the variants
variant_data %>% 
  filter(cell_type == "ASC") %>%
  select(chr, start, end, gene, mk, strand) %>%
  write.table("/GPFS/Magda_lab_temp/maitenat/scs/bullseye/6J-3M-2ws-V2-3-YTH_edit-sites_bc_filtered_all_yth_filt-all_ASC.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

variant_data %>% 
  filter(cell_type == "EC") %>%
  select(chr, start, end, gene, mk, strand) %>%
  write.table("/GPFS/Magda_lab_temp/maitenat/scs/bullseye/6J-3M-2ws-V2-3-YTH_edit-sites_bc_filtered_all_yth_filt-all_EC.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

variant_data %>% 
  filter(cell_type == "IMC") %>%
  select(chr, start, end, gene, mk, strand) %>%
  write.table("/GPFS/Magda_lab_temp/maitenat/scs/bullseye/6J-3M-2ws-V2-3-YTH_edit-sites_bc_filtered_all_yth_filt-all_IMC.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

variant_data %>% 
  filter(cell_type == "NEUR") %>%
  select(chr, start, end, gene, mk, strand) %>%
  write.table("/GPFS/Magda_lab_temp/maitenat/scs/bullseye/6J-3M-2ws-V2-3-YTH_edit-sites_bc_filtered_all_yth_filt-all_NEUR.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

variant_data %>% 
  filter(cell_type == "OLG") %>%
  select(chr, start, end, gene, mk, strand) %>%
  write.table("/GPFS/Magda_lab_temp/maitenat/scs/bullseye/6J-3M-2ws-V2-3-YTH_edit-sites_bc_filtered_all_yth_filt-all_OLG.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#### Variants with m/k<1 in all cells of the cell type
variant_data %>% 
  filter(cell_type == "ASC") %>%
  select(chr, start, end, gene, mk, strand) %>%
  group_by(chr, start, end, gene, strand) %>%
  filter(all(mk < 1)) %>%
  write.table("/GPFS/Magda_lab_temp/maitenat/scs/bullseye/6J-3M-2ws-V2-3-YTH_edit-sites_bc_filtered_all_yth_filt-all_mk_leq1_ASC.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

variant_data %>% 
  filter(cell_type == "EC") %>%
  select(chr, start, end, gene, mk, strand) %>%
  group_by(chr, start, end, gene, strand) %>%
  filter(all(mk < 1)) %>%
  write.table("/GPFS/Magda_lab_temp/maitenat/scs/bullseye/6J-3M-2ws-V2-3-YTH_edit-sites_bc_filtered_all_yth_filt-all_mk_leq1_EC.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

variant_data %>% 
  filter(cell_type == "IMC") %>%
  select(chr, start, end, gene, mk, strand) %>%
  group_by(chr, start, end, gene, strand) %>%
  filter(all(mk < 1)) %>%
  write.table("/GPFS/Magda_lab_temp/maitenat/scs/bullseye/6J-3M-2ws-V2-3-YTH_edit-sites_bc_filtered_all_yth_filt-all_mk_leq1_IMC.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

variant_data %>% 
  filter(cell_type == "NEUR") %>%
  select(chr, start, end, gene, mk, strand) %>%
  group_by(chr, start, end, gene, strand) %>%
  filter(all(mk < 1)) %>%
  write.table("/GPFS/Magda_lab_temp/maitenat/scs/bullseye/6J-3M-2ws-V2-3-YTH_edit-sites_bc_filtered_all_yth_filt-all_mk_leq1_NEUR.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

variant_data %>% 
  filter(cell_type == "OLG") %>%
  select(chr, start, end, gene, mk, strand) %>%
  group_by(chr, start, end, gene, strand) %>%
  filter(all(mk < 1)) %>%
  write.table("/GPFS/Magda_lab_temp/maitenat/scs/bullseye/6J-3M-2ws-V2-3-YTH_edit-sites_bc_filtered_all_yth_filt-all_mk_leq1_OLG.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

