library(TCGAbiolinks)
library(tidyverse)
FileNames <- c("TCGA-LAML","TCGA-ACC","TCGA-BLCA","TCGA-LGG","TCGA-BRCA","TCGA-CESC","TCGA-CHOL","TCGA-COAD","TCGA-ESCA","TCGA-GBM","TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-DLBC","TCGA-MESO","TCGA-OV","TCGA-PAAD","TCGA-PCPG","TCGA-PRAD","TCGA-READ","TCGA-SARC","TCGA-SKCM","TCGA-STAD","TCGA-TGCT","TCGA-THYM","TCGA-THCA","TCGA-UCS","TCGA-UCEC","TCGA-UVM")
for (TumorType in FileNames) {
  query <- GDCquery(project = TumorType, data.category = "Clinical")
  GDCdownload(query)
  clinical <- GDCprepare_clinic(query, clinical.info = "patient")
  TumorType <- paste("/media/nicolas/ce71f585-e4d4-43a3-adea-f042dd352508/TCGA_rnaseq/TCGA_xml_clinical_data/", TumorType, sep = "")
  write.table(clinical, file = TumorType, sep = "\t", quote = FALSE, row.names = FALSE)
}
