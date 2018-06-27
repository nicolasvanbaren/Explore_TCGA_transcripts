# Explore_TCGA_transcripts
Download, clean, organize RNA-seq and display graphically TCGA data


Downloading and organizing the TCGA RNA-seq data for data analysis in R
NvB 4j18

The aim is to obtain a data frame organized as follows:


Gene 1
Gene 2
Gene 3
...
Gene 11574
Patient annotation 1
Patient annotation 2
Patient annotation n
Sample 1
Expr. Level 1:1
Expr. Level 1:2
Expr. Level 1:3
...
Expr. Level 1:11574
...
...
...
Sample 2
Expr. Level 2:1
Expr. Level 2:2
Expr. Level 2:3
...
Expr. Level 2:11574
...
...
...
...
...
...
...
...
...
...
...
...
Sample 19285
Expr. Level 19285:1
Expr. Level 19285:2
Expr. Level 19285:3
...
Expr. Level 119285:11574
...
...
...

This data frame is a “tibble” (= an optimized data frame format according to the R tidyverse packages), which allows effective data manipulations and graphical outputs.
Genes are identified by their official HGNC symbol (ACTB, MAGEA3, …) and are restricted to protein-coding genes (n= 19285).
Samples are TCGA tumors and adjacent normal tissues for which RNA-seq data are available  (n= 11574) and are identified by a custom code combining the tumor type abbreviation and a rank number (e.g. BRCA1024).
Expression levels are the values downloaded from the TCGA database (available via the Genomic Data Commons), and correspond either to HTSeq-Counts (the number of reads aligned to each gene), HTSeq-FPKM (idem, normalized according to the FPKM method, i.e. for each gene, number of aligned reads divided by 2 (= number of aligned fragments), divided by the summed length of the gene exons, divided by the total number of reads times 1 million), or HTSeq-FPKM-UQ (idem as FPKM, with one additional normalization step: the upper quartile (P75) is set to a constant value (which one?) for each gene).
Patient annotations are clinical annotations corresponding to each sample/patient and downloaded from the TCGA database.

Remark:
Individual genes and samples appear as columns and rows, respectively, and not otherwise, because the tibble format requires the column items to have the same type (here numeric values) except the title, so as to ensure effective data manipulations, and this does not allow to add annotations in lower rows. Thus we have to choose between clinical and gene annotation. The former is more important, which necessitates to put the samples in rows.

Strategy: 


1. Download and organize the gene expression data

A. Download the data from the web
These are publicly available, processed data (gene expression values). The data come as individual files (one per sample).
Source: Genomic Data Commons, GDC Portal (https://portal.gdc.cancer.gov/)
Procedure:
a. Download a manifest file, i.e. a file that lists all the individual files that have to be downloaded.
- in the Portal web page, select “Repository”.
- in the left column, select the tab “Cases”, then select “TCGA”.
- in the left column, select the tab “Files”, choose “RNA-seq” and either “HTSeq-Counts”, “HTSeq-FPKM” or “HTSeq-FPKM-UQ” depending on which type of data you want.

- in the Data Window (right window), select “Files” (you should have 11093 cases).
- choose the “Manifest” tab. A manifest file (gdc_manifest.<date>.txt) is downloaded. Rename it  (here 01_gdc_manifest.2018-05-23.txt) and move it to your desired directory.

b. Download the individual files. To this aim:
- open a Unix Terminal window, and enter the following script (adjust paths/file names):
/apps/gdc-client  download -m  /media/nicolas/ce71f585-e4d4-43a3-adea-f042dd352508/TCGA_rnaseq/01_gdc_manifest.2018-05-23.txt 
gunzip /home/nicolas/*/*.gz 
mv /home/nicolas/*/*.counts /media/nicolas/ce71f585-e4d4-43a3-adea-f042dd352508/TCGA_rnaseq/TCGA_rnaseq_counts
This downloads thousands of files in your home folder, unzip them and transfers them in a dedicated directory.

B. Associate the data in a single table

a. Generate list of genes

Remark:
This list corresponds to the first column of each individual data file (n = 60488), where each gene is identified by its Ensembl code + suffix corresponding to version (e.g. ENSG00000155865.11). This column is the same for every individual data file. As there are often discrepancies related to the version, the best is to remove the version suffix. There are several ways to do that, I report the one I used using a spreadsheet software.
This gene list will also be useful in the next section.

- open one of the individual data files in Excel or LibreOffice Calc using <tab> as field delimiter.
- delete everything but the first column.
- in cell B1 enter “LEFT(A1,15)”, this selects the first 15 left characters of the content of cell A1.
N.B. in Excel-french the function is “GAUCHE”
- copy this function (copy cell B1) and paste it in cells B2 to B60488.
- copy the whole B column into a new document created in a text editor.
- put the cursor just before the first item, type “gene_id” (without quotes) and press “Enter”. This will be used as the name of the column in further steps.
- name and save the file (1. ENSG_gene_id_list.txt).

b. Combine the individual data files into a single data table
- launch RStudio. Open the R script “Rscript_01_combine_TCGA_files.R” or copy the code below in a new window, adjust working directory, file names and paths where appropriate. Select the script and click “Run”.
library(tidyverse)
library(stringr)
setwd("/media/nicolas/ce71f585-e4d4-43a3-adea-f042dd352508/TCGA_rnaseq")
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

This creates a data frame called "TCGA_rnaseq_FPKM_UQ.txt". Its column names are “gene_id” followed by the hexadecimal code of each sample. It contains 60483 rows, one per gene.

Remark:
At this stage there are two major obstacles:
1. The genes are identified by their Ensembl code (i.e. ENSGxxxxxxxxxxxxxxx.yy), where yy stands for the version of the gene data stored in Ensembl. It will be necessary to obtain matched gene annotations, including the HGNC gene symbol, which is often used.
2. The manifest file and the individual data files are identified by a sample_uuid (uuid = randomly generated hexadecimal code), and the clinical annotations by a case_id (= uuid for the patient) and a barcode or submitter_id, which is the previous TCGA patient identifier. There is no direct match between the codes.


2. Download the gene annotations and combine with the expression data

Procedure:
A. Select and download the corresponding gene annotations
Source: Biomart (gene annotation management database maintained by Ensembl, https://www.ensembl.org/biomart)
- choose Database: select “Ensembl genes 92”
- choose Dataset: select “Human genes (GRCh38.p12)”
- in the left window, choose “Filters” (to select the input criteria)
- select “Gene”
- select “Input external references ID list [Max 500 advised]” and then “Gene stable ID(s)”
- select “Browse” and enter the path to your Gene List file (ENSG_gene_id_list.txt, see A).
- in the left window, choose “Attributes” (to select the output items)
- select “Gene” and then	Gene stable ID*
				Gene description 		
				Chromosome/scaffold name 
				Gene start (bp) 
				Gene end (bp) 
				Strand
				Gene name 
				Gene type*
- select “External” and then	HGNC symbol*

Remark:
You may change some of the annotations above according to your needs, however those marked with * are mandatory. Gene stable ID is required for the alignment of the expression and annotation tables, Gene type to select for protein-coding genes, and HGNC symbol is the only identifier that will appear on the final dataframe.

- choose “Results” in the top left tabs. A preview is generated.
- choose “Export  all results to “File” “TSV” and tick Unique results only, click “Go”.
- A file called “mart_export.txt is downloaded. It contains the gene annotations.

B. Combine input gene list and output annotation table

Remark:
For obscure reasons, the number of genes in the Biomart output files can be less than that in the input file. To avoid further problems, we are going to fuse and align the two lists of genes. The R function used here will ensure that genes present in the input but absent in the output are represented in the annotation table, with “NA” replacing empty fields. The fusion requires that the two gene list columns have the same identifier.

- open mart_export.txt and verify that the gene column (first column, with items ENSG...) has gene_id as its name (on top). If necessary, change. The R application that fuses the two files needs a common identifier to perform the matched fusion.
- in Rstudio, open the R script “Rscript_02_join_biomart_annot.R” or copy the code below in a new window, adjust working directory, file names and paths where appropriate. Run.
library(tidyverse) 
setwd("/media/nicolas/ce71f585-e4d4-43a3-adea-f042dd352508/TCGA_rnaseq") 
first <- read_tsv("ENSG_gene_id_list.txt") 
second <- read_tsv("mart_export.txt", na = "NA") 
joined_annot <- left_join(first, second, by = "gene_id") 
write.table(joined_annot, file ="ENSG_biomart_annot.txt", sep = "\t", quote = FALSE, row.names = FALSE)
- a new annotation file (ENSG_biomart_annot.txt) is generated. It contains 60511 rows, one per gene.

C. Fuse gene expression data frame and gene annotation table

Remark:
The annotation file needs additional curation:
- we have to select the protein-coding genes only
- we have to remove potential unidentified gene rows (labeled with “NA”)
- it is expected that each gene is represented by a unique symbol. However, we noticed that the PRAMEF7 gene appeared twice as a protein-coding gene in the annotations, so we had to remove it (in row 12808).

In addition, the fused table needs to be transposed (samples as rows and genes as columns) before the next step.

- in Rstudio, open the R script “Rscript_03_combine_expr_data_gene_annot.R” or copy the code below in a new window, adjust working directory, file names and paths where appropriate. Run.
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
  gather(key = manifest_file_id, value = value, 2:ncol(aaabbb)) %>% 
  spread_(key = names(aaabbb)[1],value = 'value')
write.table(aaabbb, file = "TCGA_rnaseq_FPKM_UQ_hgnc.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
- a new transposed data frame (TCGA_rnaseq_FPKM_UQ_hgnc.txt) is generated. It contains 11574 rows, one per sample. Genes are now identified by their HGNC symbol.


3. Download the sample and clinical annotations and join them to the expression data

This is the most tricky part. There are clinical annotations and sample annotations (one patient can have more than one sample, e.g. a tumor and an adjacent normal tissue). 

A. Download the sample identifiers and sample barcodes, and derive the patient barcode identifiers

a. Download the sample barcodes corresponding to the sample identifiers
The script below is derived from a blog site. It allows to obtain both the sample IDs identical to those in the manifest file, and their correspondence with the TCGA sample barcodes.

Remark:
Example of barcode is “TCGA-JH-5621-01A”, where the items correspond to the project, center, patient ID and sample type, respectively. Barcodes have been replaced since then by hexadecimal uuid identifiers, such as “802f4a27-3cf1-4b3a-b498-9b95djb0b62a”, which are generated randomly and can identify various items such as patients, samples, treatments etc.

- in Rstudio, open the R script “Rscript_04_retrieve_barcode_from_sampleUUID.R” or copy the code below in a new window, adjust working directory, file names and paths where appropriate. Run.
library(GenomicDataCommons)
library(magrittr)
setwd("/media/nicolas/ce71f585-e4d4-43a3-adea-f042dd352508/TCGA_rnaseq")
manifest <- read.table("gdc_manifest.2018-05-23.txt", header = TRUE)
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
write.table(res, file = "corresp_sampleUUID_barcode.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
- a new data frame (corresp_sampleUUID_barcode.txt) is generated. It contains 11093 rows, one per sample.

b. Derive the patient barcode identifiers from the sample barcodes
Open the "corresp_sampleUUID_barcode.txt" file in a spreadsheet software and follow the procedure in 1.B.a: create a new column with the patient barcode (e.g. “TCGA-JH-5621”) adjacent to the truncated source sample barcode (“TCGA-JH-5621-01A”). Name this column “patient_id”. Save the new file as “corresp_fileUUID_barcode.txt”. It contains 11093 rows, one per sample.

c. fuse with the table that matches manifest_id with file_id
- <procedure> where does this fileid_manifestfileid.txt file come from?
- in Rstudio, open the R script “Rscript_05_fuse_manifest_file_id.R” or copy the code below in a new window, adjust working directory, file names and paths where appropriate. Run.
library(tidyverse)
setwd("/media/nicolas/ce71f585-e4d4-43a3-adea-f042dd352508/TCGA_rnaseq")

a <- read_tsv("fileid_manifestfileid.txt")
b <- read_tsv("corresp_fileUUID_barcode.txt")
ab <- full_join(a,b,by = "file_id")
write.table(ab, file = "manifest_file_corresp.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
It contains 11093 rows, one per sample.

B. Download, arrange the sample annotations and fuse them to the sample and patient identifiers from A

a. Download and arrange the sample annotations
- <procedure>
- we are interested only by the slide.tsv, which contains useful sample related information such as the proportion of tumor cells and the percentage of lymphocyte infiltration.
- a quick look at the file shows that individual patients (identified by case_id or barcode) can have different samples (identified by sample_submitter_id/sample barcode), corresponding to e.g. tumor and adjacent normal tissue (codes 01A and 11A, respectively). In addition, some samples are present twice, as they were separated into top and bottom portions that were analyzed separately. As single samples are labeled “top” and top and bottom samples do not seem to differ much, we are going to work with the top samples only (ideally we should calculate mean values).
- 


b. Fuse the sample annotations to the sample and patient identifiers from A

- in Rstudio, open the script “Rscript_06_join_sample_annot_id.R” or copy-paste the following code (adjust paths and file names):
library(tidyverse)
library(stringr)
setwd("/media/nicolas/ce71f585-e4d4-43a3-adea-f042dd352508/TCGA_rnaseq")

sample_annot <- read_tsv("Downloaded_from_GDC/slide.tsv")
sample_annot <- sample_annot %>%
  mutate(tumor_type_factor = factor(str_sub(project_id, 6, 9))) %>%
  rename(submitter_id = sample_submitter_id) %>%
  filter(section_location == "TOP")
  
sample_ids <- read_tsv("manifest_file_corresp.txt")

full_sample_annot <- left_join(sample_ids, sample_annot, by = "submitter_id")
full_sample_annot <- full_sample_annot %>%
  group_by(submitter_id) %>%
  slice(1) %>%
  ungroup()

write.table(full_sample_annot, file = "full_sample_annot.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
- a new completed data frame (full_sample_annot.txt) is generated. It contains 11057 rows, one per sample.

C. Download and combine the clinical annotations

Remarks:
There are two main types of patient annotations available from the GDC portal. The raw or xml  annotations and the indexed annotations. The latter have been standardized across the different types of tumors and contain less detailed information. We have focused on the former.

a. Download the raw clinical annotations
- in Rstudio, open the script “Rscript_07_download_xml_clinical_data.R” or copy-paste the following code (adjust paths and file names):
library(TCGAbiolinks)
library(tidyverse)
setwd("/media/nicolas/ce71f585-e4d4-43a3-adea-f042dd352508/TCGA_rnaseq")
FileNames <- c("TCGA-LAML","TCGA-ACC","TCGA-BLCA","TCGA-LGG","TCGA-BRCA","TCGA-CESC","TCGA-CHOL","TCGA-COAD","TCGA-ESCA","TCGA-GBM","TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-DLBC","TCGA-MESO","TCGA-OV","TCGA-PAAD","TCGA-PCPG","TCGA-PRAD","TCGA-READ","TCGA-SARC","TCGA-SKCM","TCGA-STAD","TCGA-TGCT","TCGA-THYM","TCGA-THCA","TCGA-UCS","TCGA-UCEC","TCGA-UVM")
for (TumorType in FileNames) {
  query <- GDCquery(project = TumorType, data.category = "Clinical")
  GDCdownload(query)
  clinical <- GDCprepare_clinic(query, clinical.info = "patient")
  TumorType <- paste("TCGA_xml_clinical_data/", TumorType, sep = "")
  write.table(clinical, file = TumorType, sep = "\t", quote = FALSE, row.names = FALSE)
}
- this generates, for each tumor type, a table with the clinical annotations for each patient.

b. Combine the per-tumor-type patient annotation files into one single table, fuse with the corresp_fileUUID_barcode.txt table, and add an identifier that distinguishes the type of sample (primary tumor, metastasis (for melanoma only), adjacent normal tissue)

- in Rstudio, open the script “Rscript_08_combine_clin_annot_uuid_barcode.R” or copy-paste the following code (adjust paths and file names):
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
- a new data frame (clinical_annot.txt) with transposed data (samples as rows and genes as columns) is generated. It contains 11622 rows, one per sample.

D. Combine the expression data, the sample and the patient annotations

- in Rstudio, open the script “Rscript_09_join_exprData_sampleAnnot_clinAnnot.R” or copy-paste the following code (adjust paths and file names):
library(tidyverse)
setwd("/media/nicolas/ce71f585-e4d4-43a3-adea-f042dd352508/TCGA_rnaseq")

aa <- read_tsv("TCGA_rnaseq_FPKM_UQ_hgnc.txt", guess_max = 20000)
bb <- read_tsv("full_sample_annot.txt")
aabb <- left_join(aa,bb,by = "manifest_file_id")

cc <- read_tsv("clinical_annot.txt", guess_max = 640)
aabbcc <- left_join(aabb,cc,by = "submitter_id")

write.table(aabbcc, file = "expr_data_full_clin_annot.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
- a new completed data frame (expr_data_full_clin_annot.txt) is generated. It contains 12138 rows, one per sample.

