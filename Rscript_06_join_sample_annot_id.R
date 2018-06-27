library(tidyverse)
library(stringr)
setwd("/media/nicolas/ce71f585-e4d4-43a3-adea-f042dd352508/TCGA_rnaseq")

sample_annot <- read_tsv("Downloaded_from_GDC/slide.tsv")
sample_annot <- sample_annot %>%
  mutate(tumor_type_factor = factor(str_sub(project_id, 6, 9))) %>%
  rename(submitter_id = sample_submitter_id) %>%
  filter(section_location == "TOP")
  
sample_ids <- read_tsv("corresp_fileUUID_barcode.txt")

full_sample_annot <- left_join(sample_ids, sample_annot, by = "submitter_id")
full_sample_annot <- full_sample_annot %>%
  group_by(submitter_id) %>%
  slice(1) %>%
  ungroup()

write.table(full_sample_annot, file = "full_sample_annot.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
