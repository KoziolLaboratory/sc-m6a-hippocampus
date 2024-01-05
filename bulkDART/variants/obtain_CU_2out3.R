# Description:  In this script we obtain the C>U variant lists using replicates. We keep the variants that appear in at least 2/3 replicates and set their freq to the mean of the replicates
# ———————————----------—#

#### Functions ----
read_cims <- function(file) {
  data <- read.table(file, header = FALSE) %>%
    as_tibble()
  colnames(data) <- c("chr", "start", "end", "name", "score", "strand", "k", "m", "fdr", "count", "m/k")
  return(data)
}
#### ---------------

library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
rep1 <- args[1]
rep2 <- args[2]
rep3 <- args[3]
samplename <- args[4]

rep1_data <- read_cims(rep1)
rep2_data <- read_cims(rep2)
rep3_data <- read_cims(rep3)

# Intersection of all 3
rep123 <- rep1_data %>%
  inner_join(rep2_data, by = c("chr", "start", "end", "strand")) %>%
  inner_join(rep3_data, by = c("chr", "start", "end", "strand")) %>%
  select(1, 2, 3, 4, 5, 6, 11, 18, 25) %>%
  magrittr::set_colnames(c("chr", "start", "end", "name", "score", "strand", "ratio1", "ratio2", "ratio3")) %>%
  rowwise() %>%
  mutate(mean_ratio = mean(c(ratio1, ratio2, ratio3))) %>%
  select(-c(ratio1, ratio2, ratio3))

# Pairwise intersections
rep12 <- rep1_data %>%
  inner_join(rep2_data, by = c("chr", "start", "end", "strand")) %>% 
  select(1, 2, 3, 4, 5, 6, 11, 18) %>%
  magrittr::set_colnames(c("chr", "start", "end", "name", "score", "strand", "ratio1", "ratio2")) %>%
  rowwise() %>%
  mutate(mean_ratio = mean(c(ratio1, ratio2))) %>%
  select(-c(ratio1, ratio2))

rep13 <- rep1_data %>%
  inner_join(rep3_data, by = c("chr", "start", "end", "strand")) %>% 
  select(1, 2, 3, 4, 5, 6, 11, 18) %>%
  magrittr::set_colnames(c("chr", "start", "end", "name", "score", "strand", "ratio1", "ratio2")) %>%
  rowwise() %>%
  mutate(mean_ratio = mean(c(ratio1, ratio2))) %>%
  select(-c(ratio1, ratio2))

rep23 <- rep2_data %>%
  inner_join(rep3_data, by = c("chr", "start", "end", "strand")) %>% 
  select(1, 2, 3, 4, 5, 6, 11, 18) %>%
  magrittr::set_colnames(c("chr", "start", "end", "name", "score", "strand", "ratio1", "ratio2")) %>%
  rowwise() %>%
  mutate(mean_ratio = mean(c(ratio1, ratio2))) %>%
  select(-c(ratio1, ratio2))

# Pairwise intersection without intersection of the 3
rep12 %>%
  anti_join(rep123, by = c("chr", "start", "end", "strand")) -> rep12_only

rep13 %>%
  anti_join(rep123, by = c("chr", "start", "end", "strand")) -> rep13_only

rep23 %>%
  anti_join(rep123, by = c("chr", "start", "end", "strand")) -> rep23_only

# Final list
final_df <- bind_rows(rep123, rep12_only, rep13_only, rep23_only)

outfilename <- sprintf("%s_2out3.bed", samplename)
write_tsv(final_df, outfilename, col_names = FALSE)
