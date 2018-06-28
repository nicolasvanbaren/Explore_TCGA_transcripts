library(tidyverse)
setwd("/media/nicolas/ce71f585-e4d4-43a3-adea-f042dd352508/TCGA_rnaseq")

optimizedTable <- read_tsv("expr_data_select_clin_annot.txt", col_types = cols(.default = col_double(),bcr_patient_barcode="c",tumor_type_factor="c",submitter_id="c",tissue_type="c",days_to_birth="i",gender="c",race_list="c",ethnicity="c",other_dx="c",history_of_neoadjuvant_treatment="c",person_neoplasm_cancer_status="c",vital_status="c",days_to_last_followup="i",days_to_death="i",stage_event_pathologic_stage="c",stage_event_tnm_categories="c",days_to_last_known_alive="i",tumor_tissue_site="c",histological_type="c",diagnosis_subtype="c",neoplasm_histologic_grade="c",anatomic_neoplasm_subdivision="c",days_to_initial_pathologic_diagnosis="i",age_at_initial_pathologic_diagnosis="i",year_of_initial_pathologic_diagnosis="i",metastatic_site_list="c",other_metastatic_site="c",breast_carcinoma_progesterone_receptor_status="c",breast_carcinoma_estrogen_receptor_status="c",human_papillomavirus_types="c",patient_death_reason="c",death_cause_text="c",treatment="c",microsatellite_instability="c",kras_mutation_found="c",braf_gene_analysis_result="c",loss_expression_of_mismatch_repair_proteins_by_ihc_results="c",h_pylori_infection="c",ldh1_mutation_found="c",p53_gene_analysis="c",egfr_amplication_status="c",hpv_status_by_p16_testing="c",hpv_status_by_ish_testing="c",days_from_date_of_initial_pathologic_diagnosis_to_date_of_birth="i",tumor_type="c",leukemia_french_american_british_morphology_code="c",cytogenetic_abnormalities="c",child_pugh_classification_grade="c",viral_hepatitis_serologies="c",kras_mutation_result="c",egfr_mutation_result="c",eml4_alk_translocation_result="c",breslow_depth_value="d",melanoma_ulceration_indicator="c",malignant_neoplasm_mitotic_count_rate="d",distant_metastasis_anatomic_site="c",percent_tumor_nuclei="d",percent_monocyte_infiltration="d",percent_normal_cells="d",percent_eosinophil_infiltration="d",percent_lymphocyte_infiltration="d",percent_neutrophil_infiltration="d",percent_necrosis="d",percent_granulocyte_infiltration="d",number_proliferating_cells="i",percent_stromal_cells="d",percent_inflam_infiltration="d",percent_tumor_cells="d"))

optimizedTable <- optimizedTable %>%
  mutate(patient_id = as.character(tumor_type_factor)) %>%
  filter(!is.na(submitter_id)) %>%
  filter(!is.na(tumor_type_factor)) %>%
  arrange(patient_id, submitter_id) %>%
  select(submitter_id,A1BG:ZZEF1,tumor_type_factor,tissue_type,days_to_birth,gender,race_list,ethnicity,other_dx,history_of_neoadjuvant_treatment,person_neoplasm_cancer_status,vital_status,days_to_last_followup,days_to_death,stage_event_pathologic_stage,stage_event_tnm_categories,days_to_last_known_alive,tumor_tissue_site,histological_type,diagnosis_subtype,neoplasm_histologic_grade,anatomic_neoplasm_subdivision,days_to_initial_pathologic_diagnosis,age_at_initial_pathologic_diagnosis,year_of_initial_pathologic_diagnosis,metastatic_site_list,other_metastatic_site,breast_carcinoma_progesterone_receptor_status,breast_carcinoma_estrogen_receptor_status,human_papillomavirus_types,patient_death_reason,death_cause_text,treatment,microsatellite_instability,kras_mutation_found,braf_gene_analysis_result,loss_expression_of_mismatch_repair_proteins_by_ihc_results,h_pylori_infection,ldh1_mutation_found,p53_gene_analysis,egfr_amplication_status,hpv_status_by_p16_testing,hpv_status_by_ish_testing,days_from_date_of_initial_pathologic_diagnosis_to_date_of_birth,tumor_type,leukemia_french_american_british_morphology_code,cytogenetic_abnormalities,child_pugh_classification_grade,viral_hepatitis_serologies,kras_mutation_result,egfr_mutation_result,eml4_alk_translocation_result,breslow_depth_value,melanoma_ulceration_indicator,malignant_neoplasm_mitotic_count_rate,distant_metastasis_anatomic_site,percent_tumor_nuclei,percent_monocyte_infiltration,percent_normal_cells,percent_eosinophil_infiltration,percent_lymphocyte_infiltration,percent_neutrophil_infiltration,percent_necrosis,percent_granulocyte_infiltration,number_proliferating_cells,percent_stromal_cells,percent_inflam_infiltration,percent_tumor_cells)

write.table(optimizedTable, file = "expr_data_select_clin_annot.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

group_by(optimizedTable, tumor_type_factor) %>%
  summarize(count = n())
ungroup(optimizedTable)

ggplot(subset(optimizedTable, tissue_type == c("01","11")), aes(x = tumor_type_factor, y = LRRC32, color = tissue_type, fill = tissue_type)) +
  geom_boxplot(outlier.shape = NA, color = "grey", width = 0.4, position = position_dodge(width = 1)) +
  geom_point(size = 0.3, alpha = 0.5, position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1)) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim = c(0, 2500000)) +
  scale_color_manual(values = c("red","blue")) +
  scale_fill_manual(name = "Type of tissue", labels = c("Tumor", "Adjacent non-cancerous"), values = rep(NA, 2))