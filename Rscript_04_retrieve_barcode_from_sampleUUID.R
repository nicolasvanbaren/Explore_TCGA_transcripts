library(GenomicDataCommons)
library(magrittr)
setwd("/media/nicolas/ce71f585-e4d4-43a3-adea-f042dd352508/TCGA_rnaseq")
manifest <- read.table("Downloaded_from_GDC/gdc_manifest.2018-05-23.txt", header = TRUE)
file_uuids <- manifest$id
TCGAtranslateID = function(file_ids, legacy = FALSE) {
  info <- files(legacy = legacy) %>%
    filter( ~ file_id %in% file_ids) %>%
    select('cases.samples.submitter_id') %>%
    results_all()
  id_list <- lapply(info$cases,function(a) {
    a[[1]][[1]][[1]]})
  barcodes_per_file <- sapply(id_list,length)
  return(data.frame(file_id = rep(ids(info),barcodes_per_file),
                    submitter_id = unlist(id_list)))
}
res <- TCGAtranslateID(file_uuids)
write.table(res, file = "corresp_sampleUUID_barcode2.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
