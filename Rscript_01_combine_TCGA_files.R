library(tidyverse)
library(stringr)
setwd("/media/nicolas/ce71f585-e4d4-43a3-adea-f042dd352508/TCGA_rnaseq2")
FatTibble <- read_tsv("ENSG_gene_id_list.txt", col_names = TRUE)
FileNames <- list.files(path = "TCGA_rnaseq_FPKM_UQ", pattern = ".FPKM-UQ.txt$", all.files = FALSE, full.names = TRUE)
counter <- 0
for (SampleName in FileNames) {
  ForkTibble <- read_tsv(SampleName, col_names = FALSE)
  BufferName <- str_replace(SampleName, "TCGA_rnaseq_FPKM_UQ/", "")
  BufferName <- str_replace(BufferName, ".FPKM-UQ.txt", "")
  ForkTibble <- ForkTibble %>%
    mutate(X3 = stringr::str_sub(X1,1,15)) %>%
    select(X3, X2)
  colnames(ForkTibble) <- c("gene_id", BufferName)
  FatTibble <- full_join(FatTibble, ForkTibble, by = "gene_id")
  counter <- counter + 1
  print(counter)
}
write.table(FatTibble, file = "TCGA_rnaseq_FPKM_UQ.txt", quote = FALSE, sep = "\t", row.names = FALSE)
