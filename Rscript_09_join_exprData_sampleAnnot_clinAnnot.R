library(tidyverse)
setwd("/media/nicolas/ce71f585-e4d4-43a3-adea-f042dd352508/TCGA_rnaseq")

aa <- read_tsv("TCGA_rnaseq_FPKM_UQ_hgnc.txt", guess_max = 20000)
bb <- read_tsv("full_sample_annot.txt")
aabb <- left_join(aa,bb,by = "manifest_file_id")

cc <- read_tsv("clinical_annot.txt", guess_max = 640)
aabbcc <- left_join(aabb,cc,by = "submitter_id")

write.table(aabbcc, file = "expr_data_full_clin_annot.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
