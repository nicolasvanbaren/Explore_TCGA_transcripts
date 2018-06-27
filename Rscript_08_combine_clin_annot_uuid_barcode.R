library(tidyverse)
library(stringr)
setwd("/media/nicolas/ce71f585-e4d4-43a3-adea-f042dd352508/TCGA_rnaseq")

TowerTibble <- tibble::tibble()
FileNames <- list.files(path = "TCGA_xml_clinical_data", pattern = "TCGA-", all.files = FALSE, full.names = TRUE)
for (TumorName in FileNames) {
  print(TumorName)
  BrickTibble <- read.table(TumorName, header = TRUE, sep = "\t", fill = TRUE, quote = "", colClasses = "character")
  TowerTibble <- bind_rows(TowerTibble,BrickTibble, .id = NULL)
}

DungeonTibble <- read_tsv("corresp_fileUUID_barcode.txt")
CastleTibble <- full_join(DungeonTibble, TowerTibble,by ="bcr_patient_barcode")

# next lines of code adds a new factor column with a number that identifies the type of sample (01 = primary tumor, 06 = metastatic tumor, 11 = adjacent normal tissue), and remove rows in which the sample barcode is NA.
CastleTibble <- CastleTibble %>%
  mutate(tissue_type = as.factor(str_sub(submitter_id, 14, 15))) %>%
  filter(!is.na(submitter_id))

write.table(CastleTibble, file = "clinical_annot.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


