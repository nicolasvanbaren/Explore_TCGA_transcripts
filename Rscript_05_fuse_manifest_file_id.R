library(tidyverse)
setwd("/media/nicolas/ce71f585-e4d4-43a3-adea-f042dd352508/TCGA_rnaseq")

a <- read_tsv("fileid_manifestfileid.txt")
b <- read_tsv("corresp_fileUUID_barcode.txt")
ab <- full_join(a,b,by = "file_id")
write.table(ab, file = "manifest_file_corresp.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
