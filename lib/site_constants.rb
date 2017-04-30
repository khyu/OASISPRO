module SiteConstants
	TUMOR_TYPES = {
		acc:  'Adrenocortical carcinoma',
		blca: 'Bladder urothelial carcinoma',
		lgg:  'Brain lower grade glioma',
		brca: 'Breast invasive carcinoma',
		cesc: 'Cervical squamous cell carcinoma and endocervical adenocarcinoma',
		chol: 'Cholangiocarcinoma',
		coad: 'Colon adenocarcinoma',
		esca: 'Esophageal carcinoma',
		#gbm:  'Glioblastoma multiforme',
		hnsc: 'Head and neck squamous cell carcinoma',
		kich: 'Kidney chromophobe',
		kirc: 'Kidney renal clear cell carcinoma',
		kirp: 'Kidney renal papillary cell carcinoma',
		lihc: 'Liver hepatocellular carcinoma',
		luad: 'Lung adenocarcinoma',
		lusc: 'Lung squamous cell carcinoma',
		laml: 'Lymphoblastic acute myeloid leukemia',
		dlbc: 'Lymphoid neoplasm diffuse large B-cell lymphoma',
		meso: 'Mesothelioma',
		ov:   'Ovarian serous cystadenocarcinoma', 
		paad: 'Pancreatic adenocarcinoma', 
		pcpg: 'Pheochromocytoma and Paraganglioma',
		prad: 'Prostate adenocarcinoma', 
		read: 'Rectum adenocarcinoma',
		sarc: 'Sarcoma',
		skcm: 'Skin cutaneous melanoma',
		stad: 'Stomach adenocarcinoma',
		tgct: 'Testicular germ cell tumors',
		thym: 'Thymoma',
		thca: 'Thyroid carcinoma',
		ucs:  'Uterine carcinosarcoma',
		ucec: 'Uterine corpus endometrial carcinoma',
		uvm:  'Uveal melanoma' 
	}

	DATA_TYPES = {
		rnaseq: 		'Gene expression (RNA-seq)',
		proteomics: 'Proteomics',
		methylation27: 	'DNA methylation (Illumina Infinium Human DNA Methylation 27)',
		#methylation450: 	'DNA Methylation (Illumina Infinium Human DNA Methylation 450)',
		microrna: 	'microRNA'
	}

	TUMOR_TYPES_MSSING_DATA_TYPES = {
		laml: ['proteomics', 'microrna'],
		acc: ['methylation27'],
		blca: ['methylation27'],
		cesc: ['methylation27'],
		chol: ['methylation27'],
		dlbc: ['methylation27'],
		esca: ['methylation27'],
		hnsc: ['methylation27'],
		kich: ['methylation27'],
		lcml: ['methylation27'],
		lgg: ['methylation27'],
		lihc: ['methylation27'],
		meso: ['methylation27'],
		paad: ['methylation27'],
		pcpg: ['methylation27'],
		prad: ['methylation27'],
		sarc: ['methylation27'],
		skcm: ['methylation27'],
		tgct: ['methylation27'],
		thca: ['methylation27'],
		thym: ['methylation27'],
		ucs: ['methylation27'],
		uvm: ['methylation27']
	}

	PREDICTION_TARGTS = {
		stage: 		'Tumor stage (TNM)',
		stage_t: 	'T stage',
		stage_n: 	'N stage',
		stage_m: 	'M stage' 
	}

	FEATURE_SELECTION_METHOD = {
		infog: 		'Information gain',
		gainr: 		'Gain tatio',
		symu: 		'Symmetrical uncertainty',
		randf: 		'Random forest importance',
		_features: 	'Custom' # Prepend SESSION_ID and append TXT when submitting this
	}

	PARTITION_TYPES = ['random', 'batch']

	MEDICAL_CENTERS = []
	medical_center_lines = File.read("public/notes/tissueSourceSite.txt").split("\n")[1..-1]
	medical_center_lines.each do |medical_center_line|
		medical_center_line = medical_center_line.split("\t")
		MEDICAL_CENTERS << {id: medical_center_line[0], name: medical_center_line[1]}
	end

	CLINICAL_VARIABLES = [
		nil,
		'Age',
		'Gender',
		'Race',
		'Ethnicity',
		'History Neoadjuvant Treatment',
		'Vital Status',
		'Last Contact Years To',
		'Death Days To',
		'History Thyroid Disease',
		'Family History Thyroid Cancer',
		'History Radiation Exposure',
		'Anatomic Subdivision',
		'Tumor Focality',
		'Tumor Size Length',
		'Tumor Size Width',
		'Tumor Size Depth',
		'Age At Diagnosis',
		'Extrathyroidal Extension',
		'Ajcc Tumor Pathologic pt'
	]

	CLINICAL_VARIABLE_TYPES = {
		'bcr_patient_uuid' => 'd',
		'bcr_patient_barcode' => 'd',
		'form_completion_date' => 'd',
		'tissue_prospective_collection_indicator' => 'd',
		'tissue_retrospective_collection_indicator' => 'd',
		'gender' => 'd',
		'race' => 'd',
		'ethnicity' => 'd',
		'prior_dx' => 'd',
		'history_of_neoadjuvant_treatment' => 'd',
		'person_neoplasm_cancer_status' => 'd',
		'vital_status' => 'd',
		'radiation_therapy' => 'd',
		'postoperative_rx_tx' => 'd',
		'mitotane_therapy' => 'd',
		'mitotane_therapy_adjuvant_setting' => 'd',
		'therapeutic_mitotane_levels_achieved' => 'd',
		'therapeutic_mitotane_lvl_recurrence' => 'd',
		'mitotane_therapy_for_macroscopic_residual_disease' => 'd',
		'therapeutic_mitotane_lvl_macroscopic_residual' => 'd',
		'therapeutic_mitotane_lvl_progression' => 'd',
		'post_surgical_procedure_assessment_thyroid_gland_carcinoma_status' => 'd',
		'primary_therapy_outcome_success' => 'd',
		'laterality' => 'd',
		'histological_type' => 'd',
		'year_of_initial_pathologic_diagnosis' => 'c',
		'ct_scan' => 'd',
		'ct_scan_findings' => 'd',
		'primary_lymph_node_presentation_assessment' => 'd',
		'lymph_node_examined_count' => 'c',
		'number_of_lymphnodes_positive_by_he' => 'c',
		'weiss_score' => 'c',
		'mitoses_count' => 'c',
		'pathologic_stage' => 'd',
		'residual_tumor' => 'd',
		'metastatic_neoplasm_confirmed_diagnosis_method_name' => 'd',
		'metastatic_neoplasm_confirmed_diagnosis_method_text' => 'd',
		'metastatic_neoplasm_initial_diagnosis_anatomic_site' => 'd',
		'distant_metastasis_anatomic_site' => 'd',
		'excess_adrenal_hormone_history_type' => 'd',
		'excess_adrenal_hormone_diagnosis_method_type' => 'd',
		'molecular_analysis_performed_indicator' => 'd',
		'new_tumor_event_after_initial_treatment' => 'd',
		'age_at_initial_pathologic_diagnosis' => 'c',
		'atypical_mitotic_figures' => 'd',
		'clinical_M' => 'd',
		'clinical_N' => 'd',
		'clinical_T' => 'd',
		'clinical_stage' => 'd',
		'cytoplasm_presence_less_than_equal_25_percent' => 'd',
		'days_to_birth' => 'c',
		'days_to_death' => 'c',
		'days_to_initial_pathologic_diagnosis' => 'c',
		'days_to_last_followup' => 'd',
		'diffuse_architecture' => 'd',
		'disease_code' => 'd',
		'extranodal_involvement' => 'd',
		'icd_10' => 'd',
		'icd_o_3_histology' => 'd',
		'icd_o_3_site' => 'd',
		'informed_consent_verified' => 'd',
		'invasion_of_tumor_capsule' => 'd',
		'mitotic_rate' => 'd',
		'necrosis' => 'd',
		'nuclear_grade_III_IV' => 'd',
		'pathologic_M' => 'd',
		'pathologic_N' => 'd',
		'pathologic_T' => 'd',
		'patient_id' => 'd',
		'project_code' => 'd',
		'ret' => 'd',
		'sdha' => 'd',
		'sdhaf2_sdh5' => 'd',
		'sdhb' => 'd',
		'sdhc' => 'd',
		'sdhd' => 'd',
		'sinusoid_invasion' => 'd',
		'stage_other' => 'd',
		'system_version' => 'd',
		'tissue_source_site' => 'd',
		'tmem127' => 'd',
		'tumor_tissue_site' => 'd',
		'vhl' => 'd',
		'weiss_venous_invasion' => 'd',
		'height' => 'c',
		'weight' => 'c',
		'hist_of_non_mibc' => 'd',
		'non_mibc_tx' => 'd',
		'mibc_90day_post_resection_bcg' => 'd',
		'complete_response_observed' => 'd',
		'induction_course_complete' => 'd',
		'maint_therapy_course_complete' => 'd',
		'resp_maint_from_bcg_admin_month_dur' => 'd',
		'person_occupation_description_text' => 'd',
		'occupation_primary_job' => 'd',
		'chemical_exposure_text' => 'd',
		'person_primary_industry_text' => 'd',
		'person_occupation_years_number' => 'c',
		'tobacco_smoking_history' => 'd',
		'age_began_smoking_in_years' => 'c',
		'stopped_smoking_year' => 'c',
		'number_pack_years_smoked' => 'c',
		'family_medical_history_relative_family_member_relationship_type' => 'd',
		'cancer_diagnosis_cancer_type_icd9_text_name' => 'd',
		'karnofsky_performance_score' => 'c',
		'eastern_cancer_oncology_group' => 'd',
		'diagnosis_subtype' => 'd',
		'diagnosis_age' => 'c',
		'initial_pathologic_diagnosis_method' => 'd',
		'init_pathology_dx_method_other' => 'd',
		'lymphovascular_invasion_present' => 'd',
		'disease_extracapsular_extension_ind-3' => 'd',
		'bladder_carcinoma_extracapsular_extension_status' => 'd',
		'malignant_neoplasm_metastatic_involvement_site' => 'd',
		'other_metastatic_involvement_anatomic_site' => 'd',
		'person_concomitant_prostate_carcinoma_occurrence_indicator' => 'd',
		'person_concomitant_prostate_carcinoma_pathologic_t_stage' => 'd',
		'anatomic_neoplasm_subdivision' => 'd',
		'neoplasm_histologic_grade' => 'd',
		'menopause_status' => 'd',
		'histological_type_other' => 'd',
		'breast_carcinoma_surgical_procedure_name' => 'd',
		'surgical_procedure_purpose_other_text' => 'd',
		'margin_status' => 'd',
		'breast_carcinoma_primary_surgical_procedure_name' => 'd',
		'breast_neoplasm_other_surgical_procedure_descriptive_text' => 'd',
		'breast_cancer_surgery_margin_status' => 'd',
		'axillary_lymph_node_stage_method_type' => 'd',
		'axillary_lymph_node_stage_other_method_descriptive_text' => 'd',
		'cytokeratin_immunohistochemistry_staining_method_micrometastasis_indicator' => 'd',
		'number_of_lymphnodes_positive_by_ihc' => 'd',
		'first_nonlymph_node_metastasis_anatomic_site' => 'd',
		'first_recurrent_non_nodal_metastatic_anatomic_site_descriptive_text' => 'd',
		'breast_carcinoma_estrogen_receptor_status' => 'd',
		'er_level_cell_percentage_category' => 'd',
		'breast_carcinoma_immunohistochemistry_er_pos_finding_scale' => 'd',
		'immunohistochemistry_positive_cell_score' => 'd',
		'positive_finding_estrogen_receptor_other_measurement_scale_text' => 'd',
		'er_detection_method_text' => 'd',
		'breast_carcinoma_progesterone_receptor_status' => 'd',
		'progesterone_receptor_level_cell_percent_category' => 'd',
		'breast_carcinoma_immunohistochemistry_progesterone_receptor_pos_finding_scale' => 'd',
		'breast_carcinoma_immunohistochemistry_pos_cell_score' => 'd',
		'pos_finding_progesterone_receptor_other_measurement_scale_text' => 'd',
		'pgr_detection_method_text' => 'd',
		'lab_proc_her2_neu_immunohistochemistry_receptor_status' => 'd',
		'her2_erbb_pos_finding_cell_percent_category' => 'd',
		'her2_immunohistochemistry_level_result' => 'd',
		'pos_finding_her2_erbb2_other_measurement_scale_text' => 'd',
		'her2_erbb_method_calculation_method_text' => 'd',
		'lab_procedure_her2_neu_in_situ_hybrid_outcome_type' => 'd',
		'her2_neu_breast_carcinoma_copy_analysis_input_total_number' => 'd',
		'fluorescence_in_situ_hybridization_diagnostic_procedure_chromosome_17_signal_result_range' => 'd',
		'her2_neu_and_centromere_17_copy_number_analysis_input_total_number_count' => 'd',
		'her2_neu_chromosone_17_signal_ratio_value' => 'd',
		'her2_and_centromere_17_positive_finding_other_measurement_scale_text' => 'd',
		'her2_erbb_pos_finding_fluorescence_in_situ_hybridization_calculation_method_text' => 'd',
		'metastatic_breast_carcinoma_estrogen_receptor_status' => 'd',
		'metastatic_breast_carcinoma_estrogen_receptor_level_cell_percent_category' => 'd',
		'metastatic_breast_carcinoma_immunohistochemistry_er_pos_cell_score' => 'd',
		'pos_finding_metastatic_breast_carcinoma_estrogen_receptor_other_measuremenet_scale_text' => 'd',
		'metastatic_breast_carcinoma_estrogen_receptor_detection_method_text' => 'd',
		'metastatic_breast_carcinoma_progesterone_receptor_status' => 'd',
		'metastatic_breast_carcinoma_progesterone_receptor_level_cell_percent_category' => 'd',
		'metastatic_breast_carcinoma_immunohistochemistry_pr_pos_cell_score' => 'd',
		'metastatic_breast_carcinoma_pos_finding_progesterone_receptor_other_measure_scale_text' => 'd',
		'metastatic_breast_carcinoma_progesterone_receptor_detection_method_text' => 'd',
		'metastatic_breast_carcinoma_lab_proc_her2_neu_immunohistochemistry_receptor_status' => 'd',
		'metastatic_breast_carcinoma_her2_erbb_pos_finding_cell_percent_category' => 'd',
		'metastatic_breast_carcinoma_erbb2_immunohistochemistry_level_result' => 'd',
		'metastatic_breast_carcinoma_pos_finding_her2_erbb2_other_measure_scale_text' => 'd',
		'metastatic_breast_carcinoma_her2_erbb_method_calculation_method_text' => 'd',
		'metastatic_breast_carcinoma_lab_proc_her2_neu_in_situ_hybridization_outcome_type' => 'd',
		'her2_neu_metastatic_breast_carcinoma_copy_analysis_input_total_number' => 'd',
		'metastatic_breast_carcinoma_fluorescence_in_situ_hybridization_diagnostic_proc_centromere_17_signal_result_range' => 'd',
		'her2_neu_and_centromere_17_copy_number_metastatic_breast_carcinoma_analysis_input_total_number_count' => 'd',
		'metastatic_breast_carcinoma_her2_neu_chromosone_17_signal_ratio_value' => 'd',
		'metastatic_breast_carcinoma_pos_finding_other_scale_measurement_text' => 'd',
		'metastatic_breast_carcinoma_her2_erbb_pos_finding_fluorescence_in_situ_hybridization_calculation_method_text' => 'd',
		'days_to_first_complete_response' => 'c',
		'days_to_first_partial_response' => 'c',
		'days_to_first_response' => 'c',
		'days_to_patient_progression_free' => 'c',
		'days_to_tumor_progression' => 'c',
		'er_disease_extent_prior_er_treatment' => 'd',
		'er_estimated_duration_response' => 'd',
		'er_response_type' => 'd',
		'er_solid_tumor_response_documented_type' => 'd',
		'er_solid_tumor_response_documented_type_other' => 'd',
		'field' => 'd',
		'history_of_radiation_metastatic_site' => 'd',
		'history_of_radiation_primary_site' => 'd',
		'history_prior_surgery_indicator' => 'd',
		'history_prior_surgery_type' => 'd',
		'history_prior_surgery_type_other' => 'd',
		'distant_metastasis_present_ind2' => 'd',
		'molecular_abnormality_results' => 'd',
		'molecular_abnormality_results_other' => 'd',
		'patient_progression_status' => 'd',
		'tumor_tissue_site_other' => 'd',
		'patient_death_reason' => 'd',
		'death_cause_text' => 'd',
		'birth_control_pill_history_usage_category' => 'd',
		'total_number_of_pregnancies' => 'c',
		'number_of_successful_pregnancies_which_resulted_in_at_least_1_live_birth' => 'c',
		'patient_pregnancy_spontaneous_abortion_count' => 'c',
		'patient_pregnancy_therapeutic_abortion_count' => 'c',
		'ectopic_pregnancy_count' => 'c',
		'pregnancy_stillbirth_count' => 'c',
		'pregnant_at_diagnosis' => 'd',
		'patient_history_immune_system_and_related_disorders_name' => 'd',
		'patient_history_immune_system_and_related_disorders_text' => 'd',
		'performance_status_scale_timing' => 'd',
		'performance_status_assessment_timepoint_category_other_text' => 'd',
		'keratinizing_squamous_cell_carcinoma_present_indicator' => 'd',
		'hysterectomy_performed_type' => 'd',
		'hysterectomy_performed_text' => 'd',
		'cervical_neoplasm_pathologic_margin_involved_type' => 'd',
		'cervical_neoplasm_pathologic_margin_involved_text' => 'd',
		'cervical_carcinoma_pelvic_extension_text' => 'd',
		'lymphovascular_invasion_indicator' => 'd',
		'cervical_carcinoma_corpus_uteri_involvement_indicator' => 'd',
		'lymph_node_location_positive_pathology_name' => 'd',
		'lymph_node_location_positive_pathology_text' => 'd',
		'standardized_uptake_value_cervix_uteri_assessment_measurement' => 'd',
		'fdg_or_ct_pet_performed_outcome' => 'd',
		'diagnostic_mri_result_outcome' => 'd',
		'diagnostic_ct_result_outcome' => 'd',
		'human_papillomavirus_type' => 'd',
		'human_papillomavirus_other_type_text' => 'd',
		'human_papillomavirus_laboratory_procedure_performed_name' => 'd',
		'human_papillomavirus_laboratory_procedure_performed_text' => 'd',
		'oligonucleotide_primer_pair_laboratory_procedure_performed_name' => 'd',
		'oligonucleotide_primer_pair_laboratory_procedure_performed_text' => 'd',
		'laboratory_procedure_tumor_marker_squamous_cell_carcinoma_antigen_result_value' => 'd',
		'radiation_therapy_not_administered_reason' => 'd',
		'radiation_therapy_not_administered_specify' => 'd',
		'brachytherapy_method_type' => 'd',
		'brachytherapy_method_other_specify_text' => 'd',
		'brachytherapy_first_reference_point_administered_total_dose' => 'c',
		'rt_administered_type' => 'd',
		'radiation_type_notes' => 'd',
		'rt_pelvis_administered_total_dose' => 'c',
		'external_beam_radiation_therapy_administered_paraaortic_region_lymph_node_dose' => 'c',
		'chemotherapy_negation_radiation_therapy_concurrent_not_administered_reason' => 'd',
		'chemotherapy_negation_radiation_therapy_concurrent_administered_text' => 'd',
		'chemotherapy_regimen_type' => 'd',
		'other_chemotherapy_agent_administration_specify' => 'd',
		'concurrent_chemotherapy_dose' => 'c',
		'dose_frequency_text' => 'd',
		'agent_total_dose_count' => 'c',
		'cd4_counts_at_diagnosis' => 'c',
		'cdc_hiv_risk_group' => 'd',
		'days_to_diagnostic_computed_tomography_performed' => 'c',
		'days_to_diagnostic_mri_performed' => 'c',
		'days_to_fdg_or_ct_pet_performed' => 'c',
		'days_to_hiv_diagnosis' => 'c',
		'days_to_laboratory_procedure_tumor_marker_squamous_cell_carcinoma_antigen_result' => 'c',
		'days_to_performance_status_assessment' => 'c',
		'days_to_sample_procurement' => 'c',
		'hbv_test' => 'd',
		'hcv_test' => 'd',
		'history_immunological_disease' => 'd',
		'history_immunological_disease_other' => 'd',
		'history_immunosuppresive_dx' => 'd',
		'history_immunosuppressive_dx_other' => 'd',
		'history_of_other_malignancy' => 'd',
		'history_relevant_infectious_dx' => 'd',
		'history_relevant_infectious_dx_other' => 'd',
		'hiv_rna_load_at_diagnosis' => 'd',
		'hiv_status' => 'd',
		'hpv_test' => 'd',
		'hysterectomy_Performed_Ind-3' => 'd',
		'kshv_hhv8_test' => 'd',
		'live_birth_number' => 'd',
		'lost_follow_up' => 'd',
		'nadir_cd4_counts' => 'd',
		'on_haart_therapy_at_cancer_diagnosis' => 'd',
		'on_haart_therapy_prior_to_cancer_diagnosis' => 'd',
		'patient_pregnancy_count' => 'd',
		'female_breast_feeding_or_pregnancy_status_indicator' => 'd',
		'prescribed_dose_units' => 'd',
		'prior_aids_conditions' => 'd',
		'targeted_molecular_therapy' => 'd',
		'tumor_response_cdus_type' => 'd',
		'relative_family_cancer_history' => 'd',
		'cancer_first_degree_relative' => 'c',
		'family_member_relationship_type' => 'd',
		'family_cancer_type_txt' => 'd',
		'hist_hepato_carc_fact' => 'd',
		'hist_hepato_carcinoma_risk' => 'd',
		'post_op_ablation_embolization_tx' => 'd',
		'specimen_collection_method_name' => 'd',
		'surgical_procedure_name' => 'd',
		'vascular_tumor_cell_type' => 'd',
		'perineural_invasion_present' => 'd',
		'child_pugh_classification_grade' => 'd',
		'ca_19_9_level' => 'c',
		'ca_19_9_level_lower' => 'c',
		'ca_19_9_level_upper' => 'c',
		'fetoprotein_outcome_value' => 'c',
		'fetoprotein_outcome_lower_limit' => 'c',
		'fetoprotein_outcome_upper_limit' => 'c',
		'platelet_result_count' => 'c',
		'platelet_result_lower_limit' => 'c',
		'platelet_result_upper_limit' => 'c',
		'prothrombin_time_result_value' => 'c',
		'inter_norm_ratio_lower_limit' => 'c',
		'intern_norm_ratio_upper_limit' => 'c',
		'albumin_result_specified_value' => 'c',
		'albumin_result_lower_limit' => 'c',
		'albumin_result_upper_limit' => 'c',
		'bilirubin_upper_limit' => 'c',
		'bilirubin_lower_limit' => 'c',
		'total_bilirubin_upper_limit' => 'c',
		'creatinine_value_in_mg_dl' => 'c',
		'creatinine_lower_level' => 'c',
		'creatinine_upper_limit' => 'c',
		'fibrosis_ishak_score' => 'd',
		'cholangitis_tissue_evidence' => 'd',
		'preoperative_pretreatment_cea_level' => 'c',
		'non_nodal_tumor_deposits' => 'd',
		'circumferential_resection_margin' => 'd',
		'venous_invasion' => 'd',
		'lymphatic_invasion' => 'd',
		'microsatellite_instability' => 'd',
		'number_of_loci_tested' => 'd',
		'number_of_abnormal_loci' => 'd',
		'loss_expression_of_mismatch_repair_proteins_by_ihc' => 'd',
		'loss_expression_of_mismatch_repair_proteins_by_ihc_result' => 'd',
		'kras_gene_analysis_performed' => 'd',
		'kras_mutation_found' => 'd',
		'kras_mutation_codon' => 'd',
		'braf_gene_analysis_performed' => 'd',
		'braf_gene_analysis_result' => 'd',
		'synchronous_colon_cancer_present' => 'd',
		'history_of_colon_polyps' => 'd',
		'colon_polyps_present' => 'd',
		'number_of_first_degree_relatives_with_cancer_diagnosis' => 'c',
		'anatomic_neoplasm_subdivision_other' => 'd',
		'year_of_tobacco_smoking_onset' => 'd',
		'nodal_anatomic_site' => 'd',
		'lymph_node_involvement_site' => 'd',
		'extranodal_disease_involvement_site' => 'd',
		'extranodal_disease_involvement_site_other' => 'd',
		'extranodal_involvment_site_other' => 'd',
		'number_of_involved_extranodal_sites' => 'c',
		'extranodal_sites_involvement_number' => 'd',
		'percentage_of_follicular_component' => 'd',
		'follicular_component_percent' => 'c',
		'hiv_positive_status' => 'd',
		'maximum_tumor_dimension' => 'c',
		'tumor_resected_max_dimension' => 'c',
		'maximum_tumor_bulk_anatomic_location' => 'd',
		'ldh_lab_value' => 'c',
		'ldh_level' => 'c',
		'ldh_upper_limit' => 'c',
		'ldh_norm_range_upper' => 'c',
		'bone_marrow_biopsy_done' => 'd',
		'bone_marrow_involvement' => 'd',
		'histology_of_bone_marrow_sample' => 'd',
		'bone_marrow_sample_histology' => 'd',
		'immunophenotypic_analysis_test' => 'd',
		'immunophenotypic_analysis_tested' => 'd',
		'immunophenotypic_analysis_methodology' => 'd',
		'immunophenotypic_analysis_method' => 'd',
		'result_of_immunophenotypic_analysis' => 'd',
		'immunophenotypic_analysis_results' => 'd',
		'mib1_positive_percentage_range' => 'd',
		'b_lymphocyte_genotyping_method' => 'd',
		'igh_genotype_results' => 'd',
		'igk_genotype_results' => 'd',
		'genetic_abnormality_tested' => 'd',
		'other_genetic_abnormality_tested' => 'd',
		'genetic_abnormality_tested_other' => 'd',
		'abnormality_tested_methodology' => 'd',
		'genetic_abnormality_method' => 'd',
		'abnormality_tested_results' => 'd',
		'genetic_abnormality_results' => 'd',
		'prior_immunologic_disease_type' => 'd',
		'prior_immunologic_disease_other' => 'd',
		'prior_immunosuppressive_therapy_type' => 'd',
		'prior_immunosuppressive_therapy_other' => 'd',
		'prior_infectious_disease' => 'd',
		'prior_infectious_disease_other' => 'd',
		'ebv_antibody_status' => 'd',
		'epstein_barr_viral_status' => 'd',
		'percent_ebv_positive_malignant_cells' => 'c',
		'ebv_positive_malignant_cells_percent' => 'c',
		'ebv_diagnostic_methodology' => 'd',
		'ebv_status_malignant_cells_method' => 'd',
		'pet_scan_results' => 'd',
		'genetic_abnormality_method_other' => 'd',
		'genetic_abnormality_results_other' => 'd',
		'maximum_tumor_bulk_anatomic_site' => 'd',
		'pos_lymph_node_location' => 'd',
		'pos_lymph_node_location_other' => 'd',
		'country_of_birth' => 'd',
		'country_of_procurement' => 'd',
		'state_province_of_procurement' => 'd',
		'city_of_procurement' => 'd',
		'frequency_of_alcohol_consumption' => 'd',
		'amount_of_alcohol_consumption_per_day' => 'd',
		'reflux_history' => 'd',
		'antireflux_treatment_type' => 'd',
		'h_pylori_infection' => 'd',
		'initial_diagnosis_by' => 'd',
		'barretts_esophagus' => 'd',
		'goblet_cells_present' => 'd',
		'history_of_esophageal_cancer' => 'd',
		'number_of_relatives_diagnosed' => 'd',
		'esophageal_tumor_cental_location' => 'd',
		'esophageal_tumor_involvement_site' => 'd',
		'columnar_metaplasia_present' => 'd',
		'columnar_mucosa_goblet_cell_present' => 'd',
		'columnar_mucosa_dysplasia' => 'd',
		'lymph_node_metastasis_radiographic_evidence' => 'd',
		'planned_surgery_status' => 'd',
		'treatment_prior_to_surgery' => 'd',
		'alcohol_history_documented' => 'd',
		'NA' => 'd',
		'lymphnode_neck_dissection' => 'd',
		'lymphnode_dissection_method_left' => 'd',
		'lymphnode_dissection_method_right' => 'd',
		'p53_gene_analysis' => 'd',
		'egfr_amplication_status' => 'd',
		'presence_of_pathological_nodal_extracapsular_spread' => 'd',
		'hpv_status_by_p16_testing' => 'd',
		'hpv_status_by_ish_testing' => 'd',
		'presence_of_sarcomatoid_features' => 'd',
		'percent_tumor_sarcomatoid' => 'd',
		'number_of_lymphnodes_positive' => 'd',
		'lactate_dehydrogenase_result' => 'd',
		'serum_calcium_result' => 'd',
		'hemoglobin_result' => 'd',
		'platelet_qualitative_result' => 'd',
		'white_cell_count_result' => 'd',
		'erythrocyte_sedimentation_rate_result' => 'd',
		'tumor_type' => 'd',
		'prior_hematologic_disorder_diagnosis_indicator' => 'd',
		'hydroxyurea_administration_prior_registration_clinical_study_indicator' => 'd',
		'hydroxyurea_agent_administered_day_count' => 'c',
		'cumulative_agent_total_dose' => 'c',
		'person_history_nonmedical_leukemia_causing_agent_type' => 'd',
		'person_history_leukemogenic_agent_other_exposure_text' => 'd',
		'leukemia_specimen_cell_source_type' => 'd',
		'lab_procedure_blast_cell_outcome_percentage_value' => 'c',
		'leukemia_french_american_british_morphology_code' => 'd',
		'immunophenotype_cytochemistry_testing_result' => 'd',
		'immunophenotype_cytochemistry_percent_positive' => 'c',
		'lab_procedure_bone_marrow_cellularity_outcome_percent_value' => 'c',
		'lab_procedure_leukocyte_result_unspecified_value' => 'c',
		'lab_procedure_hemoglobin_result_specified_value' => 'c',
		'lab_procedure_hematocrit_outcome_percent_value' => 'c',
		'lab_procedure_platelet_result_specified_value' => 'c',
		'lab_procedure_bone_marrow_blast_cell_outcome_percent_value' => 'c',
		'lab_procedure_bone_marrow_promyelocyte_result_percent_value' => 'c',
		'lab_procedure_bone_marrow_myelocyte_result_percent_value' => 'c',
		'lab_procedure_bone_marrow_metamyelocyte_result_value' => 'c',
		'lab_procedure_bone_marrow_band_cell_result_percent_value' => 'c',
		'lab_procedure_bone_marrow_neutrophil_result_percent_value' => 'c',
		'lab_procedure_bone_marrow_lab_eosinophil_result_percent_value' => 'c',
		'lab_procedure_bone_marrow_basophil_result_percent_value' => 'c',
		'lab_procedure_bone_marrow_lymphocyte_outcome_percent_value' => 'c',
		'lab_procedure_monocyte_result_percent_value' => 'c',
		'lab_procedure_bone_marrow_prolymphocyte_result_percent_value' => 'c',
		'lab_procedure_bone_marrow_promonocyte_count_result_percent_value' => 'c',
		'lab_procedure_abnormal_lymphocyte_result_percent_value' => 'c',
		'cytogenetic_analysis_performed_ind' => 'd',
		'fluorescence_in_situ_hybrid_cytogenetics_metaphase_nucleus_result_count' => 'c',
		'acute_myeloid_leukemia_calgb_cytogenetics_risk_category' => 'd',
		'cytogenetic_abnormality' => 'd',
		'cytogenetic_abnormality_other' => 'd',
		'fish_evaluation_performed_ind' => 'd',
		'fluorescence_in_situ_hybridization_abnormal_result_indicator' => 'd',
		'FISH_test_component' => 'd',
		'FISH_test_component_percentage_value' => 'c',
		'disease_detection_molecular_analysis_method_type' => 'd',
		'disease_detection_molecular_analysis_method_type_other_text' => 'd',
		'molecular_analysis_abnormal_result_indicator' => 'd',
		'molecular_analysis_abnormality_testing_result' => 'd',
		'molecular_analysis_abnormality_testing_percentage_value' => 'c',
		'atra_exposure' => 'd',
		'lab_procedure_bone_marrow_diff_not_reported_reason' => 'd',
		'total_dose_units' => 'd',
		'steroid_therapy_administered' => 'd',
		'tumor_location' => 'd',
		'supratentorial_localization' => 'd',
		'history_ionizing_rt_to_head' => 'd',
		'seizure_history' => 'd',
		'headache_history' => 'd',
		'mental_status_changes' => 'd',
		'visual_changes' => 'd',
		'sensory_changes' => 'd',
		'motor_movement_changes' => 'd',
		'first_presenting_symptom' => 'd',
		'first_presenting_symptom_longest_duration' => 'd',
		'asthma_history' => 'd',
		'eczema_history' => 'd',
		'hay_fever_history' => 'd',
		'mold_or_dust_allergy_history' => 'd',
		'first_diagnosis_age_asth_ecz_hay_fev_mold_dust' => 'd',
		'food_allergy_history' => 'd',
		'food_allergy_types' => 'd',
		'first_diagnosis_age_of_food_allergy' => 'd',
		'animal_insect_allergy_history' => 'd',
		'animal_insect_allergy_types' => 'd',
		'first_diagnosis_age_of_animal_insect_allergy' => 'd',
		'preoperative_corticosteroids' => 'd',
		'preoperative_antiseizure_meds' => 'd',
		'family_history_of_cancer' => 'd',
		'family_history_of_primary_brain_tumor' => 'd',
		'ldh1_mutation_tested' => 'd',
		'ldh1_mutation_test_method' => 'd',
		'ldh1_mutation_found' => 'd',
		'inherited_genetic_syndrome_found' => 'd',
		'inherited_genetic_syndrome_result' => 'd',
		'days_to_initial_score_performance_status_scale' => 'c',
		'relative_family_cancer_history_ind_3' => 'd',
		'cancer_diagnosis_first_degree_relative_number' => 'c',
		'history_hepato_carcinoma_risk_factor' => 'd',
		'history_hepato_carcinoma_risk_factor_other' => 'd',
		'viral_hepatitis_serology' => 'd',
		'ablation_embolization_tx_adjuvant' => 'd',
		'surgical_procedure_name_other_specify_text' => 'd',
		'vascular_tumor_cell_invasion_type' => 'd',
		'laboratory_procedure_alpha_fetoprotein_outcome_value' => 'c',
		'laboratory_procedure_alpha_fetoprotein_outcome_lower_limit_of_normal_value' => 'c',
		'laboratory_procedure_alpha_fetoprotein_outcome_upper_limit_of_normal_value' => 'c',
		'laboratory_prcoedure_platelet_result_lower_limit_of_normal_value' => 'c',
		'laboratory_prcoedure_platelet_result_upper_limit_of_normal_value' => 'c',
		'laboratory_procedure_prothrombin_time_result_value' => 'c',
		'laboratory_procedure_international_normalization_ratio_result_lower_limit_of_normal_value' => 'c',
		'laboratory_procedure_international_normalization_ratio_result_upper_limit_of_normal_value' => 'c',
		'laboratory_procedure_albumin_result_specified_value' => 'c',
		'laboratory_procedure_albumin_result_lower_limit_of_normal_value' => 'c',
		'laboratory_procedure_albumin_result_upper_limit_of_normal_value' => 'c',
		'laboratory_procedure_total_bilirubin_result_upper_limit_normal_value' => 'c',
		'laboratory_procedure_total_bilirubin_result_specified_lower_limit_of_normal_value' => 'c',
		'laboratory_procedure_total_bilirubin_result_specified_upper_limit_of_normal_value' => 'c',
		'hematology_serum_creatinine_laboratory_result_value_in_mg_dl' => 'c',
		'laboratory_procedure_creatinine_result_lower_limit_of_normal_value' => 'c',
		'laboratory_procedure_creatinine_result_upper_limit_of_normal_value' => 'c',
		'liver_fibrosis_ishak_score_category' => 'd',
		'adjacent_hepatic_tissue_inflammation_extent_type' => 'd',
		'days_to_definitive_surgical_procedure_performed' => 'c',
		'diagnosis' => 'd',
		'location_in_lung_parenchyma' => 'd',
		'pulmonary_function_test_performed' => 'd',
		'pre_bronchodilator_fev1_percent' => 'd',
		'post_bronchodilator_fev1_percent' => 'd',
		'pre_bronchodilator_fev1_fvc_percent' => 'd',
		'post_bronchodilator_fev1_fvc_percent' => 'd',
		'dlco_predictive_percent' => 'd',
		'kras_mutation_result' => 'd',
		'egfr_mutation_performed' => 'd',
		'egfr_mutation_result' => 'd',
		'eml4_alk_translocation_performed' => 'd',
		'eml4_alk_translocation_method' => 'd',
		'egfr_mutation_identified' => 'd',
		'eml4_alk_translocation_identified' => 'd',
		'eml4_alk_translocation_result' => 'd',
		'pleurodesis_performed_prior' => 'd',
		'pleurodesis_performed_90_days' => 'd',
		'history_asbestos_exposure' => 'd',
		'asbestos_exposure_type' => 'd',
		'asbestos_exposure_source' => 'd',
		'asbestos_exposure_age' => 'c',
		'asbestos_exposure_years' => 'c',
		'asbestos_exposure_age_last' => 'c',
		'primary_occupation' => 'd',
		'primary_occupation_other' => 'd',
		'primary_occupation_years_worked' => 'c',
		'family_history_cancer_type' => 'd',
		'family_history_cancer_type_other' => 'd',
		'serum_mesothelin_prior_tx' => 'd',
		'serum_mesothelin_lower_limit' => 'c',
		'serum_mesothelin_upper_limit' => 'c',
		'creatinine_prior_tx' => 'c',
		'creatinine_norm_range_lower' => 'c',
		'creatinine_norm_range_upper' => 'c',
		'suv_of_pleura_max' => 'c',
		'mesothelioma_detection_method' => 'd',
		'jewish_origin' => 'd',
		'tumor_residual_disease' => 'd',
		'adenocarcinoma_invasion' => 'd',
		'surgery_performed_type' => 'd',
		'surgery_performed_type_other' => 'd',
		'histologic_grading_tier_category' => 'd',
		'source_of_patient_death_reason' => 'd',
		'alcoholic_exposure_category' => 'd',
		'amt_alcohol_consumption_per_day' => 'c',
		'history_of_diabetes' => 'd',
		'days_to_diabetes_onset' => 'c',
		'history_of_chronic_pancreatitis' => 'd',
		'days_to_pancreatitis_onset' => 'c',
		'relative_cancer_type' => 'd',
		'history_pheo_or_para_include_benign' => 'd',
		'history_pheo_or_para_anatomic_site' => 'd',
		'outside_adrenal' => 'd',
		'disease_detected_on_screening' => 'd',
		'zone_of_origin' => 'd',
		'primary_pattern' => 'd',
		'secondary_pattern' => 'd',
		'gleason_score' => 'c',
		'tertiary_pattern' => 'd',
		'tumor_level' => 'd',
		'days_to_bone_scan_performed' => 'c',
		'bone_scan_results' => 'd',
		'diagnostic_ct_abd_pelvis_performed' => 'd',
		'diagnostic_ct_abd_pelvis_result' => 'd',
		'diagnostic_mri_performed' => 'd',
		'diagnostic_mri_result' => 'd',
		'lymphnodes_examined' => 'd',
		'number_of_lymphnodes_examined' => 'c',
		'days_to_psa' => 'c',
		'psa_value' => 'c',
		'biochemical_recurrence' => 'd',
		'days_to_first_biochemical_recurrence' => 'c',
		'leiomyosarcoma_histologic_subtype' => 'd',
		'primary_tumor_lower_uterus_segment' => 'd',
		'leiomyosarcoma_major_vessel_involvement' => 'd',
		'ss18_ssx_fusion_status' => 'd',
		'ss18_ssx_testing_method' => 'd',
		'mpnst_neurofibromatosis' => 'd',
		'mpnst_neurofibromatosis_heredity' => 'd',
		'mpnst_exisiting_plexiform_neurofibroma' => 'd',
		'mpnst_nf1_genetic_testing' => 'd',
		'mpnst_specific_mutations' => 'd',
		'tumor_depth' => 'c',
		'tumor_total_necrosis_percent' => 'd',
		'specific_tumor_total_necrosis_percent' => 'c',
		'mitotic_count' => 'c',
		'tumor_multifocal' => 'd',
		'discontiguous_lesion_count' => 'c',
		'radiologic_tumor_burden' => 'd',
		'pathologic_tumor_burden' => 'd',
		'local_disease_recurrence' => 'd',
		'metastatic_neoplasm_confirmed' => 'd',
		'metastatic_site_at_diagnosis' => 'd',
		'metastatic_site_at_diagnosis_other' => 'd',
		'contiguous_organ_resection_site' => 'd',
		'other_contiguous_organ_resection_site' => 'd',
		'contiguous_organ_invaded' => 'd',
		'well_differentiated_liposarcoma_primary_dx' => 'd',
		'days_to_well_differentiated_liposarcoma_primary_dx' => 'd',
		'days_to_well_differentiated_liposarcoma_resection' => 'd',
		'radiologic_tumor_length' => 'c',
		'radiologic_tumor_width' => 'c',
		'radiologic_tumor_depth' => 'c',
		'pathologic_tumor_length' => 'c',
		'pathologic_tumor_width' => 'c',
		'pathologic_tumor_depth' => 'c',
		'primary_neoplasm_melanoma_dx' => 'd',
		'primary_tumor_multiple_present_ind' => 'd',
		'primary_melanoma_at_diagnosis_count' => 'c',
		'breslow_depth_value' => 'c',
		'melanoma_clark_level_value' => 'd',
		'melanoma_ulceration_indicator' => 'd',
		'malignant_neoplasm_mitotic_count_rate' => 'c',
		'days_to_submitted_specimen_dx' => 'c',
		'melanoma_origin_skin_anatomic_site' => 'd',
		'prior_systemic_therapy_type' => 'd',
		'interferon_90_day_prior_excision_admin_indicator' => 'd',
		'new_tumor_dx_prior_submitted_specimen_dx' => 'd',
		'primary_anatomic_site_count' => 'c',
		'state_province_country_of_procurement' => 'd',
		'antireflux_treatment' => 'd',
		'family_history_of_stomach_cancer' => 'd',
		'number_of_relatives_with_stomach_cancer' => 'd',
		'history_of_undescended_testis' => 'd',
		'level_of_non_descent' => 'd',
		'undescended_testis_corrected' => 'd',
		'undescended_testis_corrected_age' => 'd',
		'undescended_testis_method_left' => 'd',
		'undescended_testis_method_right' => 'd',
		'history_hypospadias' => 'd',
		'history_fertility' => 'd',
		'family_history_testicular_cancer' => 'd',
		'relation_testicular_cancer' => 'd',
		'family_history_other_cancer' => 'd',
		'relative_family_cancer_hx_text' => 'd',
		'postoperative_tx' => 'd',
		'bilateral_diagnosis_timing_type' => 'd',
		'synchronous_tumor_histology_type' => 'd',
		'synchronous_tumor_histology_pct' => 'c',
		'testis_tumor_macroextent' => 'd',
		'testis_tumor_macroextent_other' => 'd',
		'testis_tumor_microextent' => 'd',
		'histological_percentage' => 'c',
		'intratubular_germ_cell_neoplasm' => 'd',
		'serum_markers' => 'd',
		'pre_orchi_hcg' => 'c',
		'pre_orchi_afp' => 'c',
		'pre_orchi_lh' => 'd',
		'pre_orchi_testosterone' => 'd',
		'post_orchi_ldh' => 'c',
		'post_orchi_hcg' => 'c',
		'post_orchi_afp' => 'c',
		'post_orchi_lh' => 'c',
		'post_orchi_testosterone' => 'c',
		'post_orchi_lymph_node_dissection' => 'd',
		'first_treatment_success' => 'd',
		'molecular_test_result' => 'd',
		'igcccg_stage' => 'd',
		'days_to_bilateral_tumor_dx' => 'd',
		'days_to_pre_orchi_serum_test' => 'c',
		'days_to_post_orchi_serum_test' => 'c',
		'pre_orchi_ldh' => 'c',
		'patient_personal_medical_history_thyroid_gland_disorder_name' => 'd',
		'patient_personal_medical_history_thyroid_other_specify_text' => 'd',
		'first_degree_relative_history_thyroid_gland_carcinoma_diagnosis_relationship_type' => 'd',
		'person_lifetime_risk_radiation_exposure_indicator' => 'd',
		'primary_thyroid_gland_neoplasm_location_anatomic_site' => 'd',
		'primary_neoplasm_focus_type' => 'd',
		'neoplasm_length' => 'c',
		'neoplasm_width' => 'c',
		'neoplasm_depth' => 'c',
		'lymph_node_preoperative_scan_indicator' => 'd',
		'lymph_node_preoperative_assessment_diagnostic_imaging_type' => 'd',
		'extrathyroid_carcinoma_present_extension_status' => 'd',
		'genotype_analysis_performed_indicator' => 'd',
		'braf_gene_genotyping_outcome_lab_results_text' => 'd',
		'ras_family_gene_genotyping_outcome_lab_results_text' => 'd',
		'ret_ptc_rearrangement_genotyping_outcome_lab_results_text' => 'd',
		'other_genotyping_outcome_lab_results_text' => 'd',
		'i_131_total_administered_preparation_technique' => 'd',
		'i_131_first_administered_dose' => 'c',
		'i_131_subsequent_administered_dose' => 'c',
		'i_131_total_administered_dose' => 'c',
		'radiation_therapy_administered_preparation_technique_text' => 'd',
		'radiation_therapy_administered_dose_text' => 'd',
		'radiosensitizing_agent_administered_indicator' => 'd',
		'genotyping_results_gene_mutation_not_reported_reason' => 'd',
		'masaoka_stage' => 'd',
		'history_myasthenia_gravis' => 'd',
		'section_myasthenia_gravis' => 'd',
		'horm_ther' => 'd',
		'prior_tamoxifen_administered_usage_category' => 'd',
		'hypertension' => 'd',
		'diabetes' => 'd',
		'pregnancies' => 'c',
		'colorectal_cancer' => 'd',
		'surgical_approach' => 'd',
		'peritoneal_wash' => 'd',
		'pct_tumor_invasion' => 'c',
		'total_pelv_lnr' => 'c',
		'pln_pos_light_micro' => 'c',
		'pln_pos_ihc' => 'c',
		'total_pelv_lnp' => 'c',
		'total_aor_lnr' => 'c',
		'aln_pos_light_micro' => 'c',
		'aln_pos_ihc' => 'c',
		'total_aor-lnp' => 'c',
		'eye_color' => 'd',
		'tumor_morphology_percentage' => 'd',
		'gene_expression_profile' => 'd',
		'pet_ct_standardized_values' => 'c',
		'extravascular_matrix_patterns' => 'd',
		'microvascular_density_mvd' => 'd',
		'tumor_infiltrating_lymphocytes' => 'd',
		'tumor_infiltrating_macrophages' => 'd',
		'tumor_basal_diameter' => 'c',
		'tumor_basal_diameter_mx' => 'd',
		'tumor_thickness' => 'c',
		'tumor_thickness_measurement' => 'd',
		'extrascleral_extension' => 'd',
		'extranocular_nodule_size' => 'c',
		'tumor_shape_pathologic_clinical' => 'd',
		'metastatic_site' => 'd',
		'other_metastatic_site' => 'd'
	}
end

