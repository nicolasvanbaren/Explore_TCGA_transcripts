library (stringr)
library(tidyverse)
setwd("/media/nicolas/ce71f585-e4d4-43a3-adea-f042dd352508/TCGA_rnaseq")
aaa <- read_tsv("ENSG_biomart_annot.txt")
aaa <- aaa %>%
  arrange(gene_id) %>%
  filter(Gene_type == "protein_coding") %>%
  filter(!is.na(HGNC_symbol)) %>%
  select(gene_id, HGNC_symbol) %>%
  arrange(HGNC_symbol) %>%
  slice(-12808)
bbb <- read_tsv("TCGA_rnaseq_FPKM_UQ.txt")
bbb <- arrange(bbb, gene_id)
aaabbb <- left_join(aaa, bbb, by = "gene_id")
aaabbb <- select(aaabbb, 2:11576)
aaabbb <- arrange(aaabbb, HGNC_symbol)
aaabbb <- aaabbb %>%
  gather(key = hgnc, value = value, 2:ncol(aaabbb)) %>% 
  spread_(key = names(aaabbb)[1],value = 'value')
aaabbb <- rename(aaabbb, manifest_file_id = hgnc)
write.table(aaabbb, file = "TCGA_rnaseq_FPKM_UQ_hgnc.txt", sep = "\t")
ccc <- read_tsv("Full_patient_annot_curated.txt")
ccc <- arrange(ccc, manifest_file_id)
aaabbbccc <- full_join(aaabbb, ccc, by = "manifest_file_id")
write.table(aaabbbccc, file = "TCGA_rnaseq_FPKM_UQ_annot.txt", sep = "\t")
