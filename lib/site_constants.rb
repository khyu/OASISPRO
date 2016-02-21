module SiteConstants
	TUMOR_TYPES = {
		laml: 'Acute Myeloid Leukemia',
		acc:  'Adrenocortical carcinoma',
		blca: 'Bladder Urothelial Carcinoma',
		lgg:  'Brain Lower Grade Glioma',
		brca: 'Breast invasive carcinoma',
		cesc: 'Cervical squamous cell carcinoma and endocervical adenocarcinoma',
		chol: 'Cholangiocarcinoma',
		coad: 'Colon adenocarcinoma',
		esca: 'Esophageal carcinoma',
		gbm:  'Glioblastoma multiforme',
		hnsc: 'Head and Neck squamous cell carcinoma',
		kich: 'Kidney Chromophobe',
		kirc: 'Kidney renal clear cell carcinoma',
		kirp: 'Kidney renal papillary cell carcinoma',
		lihc: 'Liver hepatocellular carcinoma',
		luad: 'Lung adenocarcinoma',
		lusc: 'Lung squamous cell carcinoma',
		dlbc: 'Lymphoid Neoplasm Diffuse Large B-cell Lymphoma',
		meso: 'Mesothelioma',
		ov:   'Ovarian serous cystadenocarcinoma', 
		paad: 'Pancreatic adenocarcinoma', 
		pcpg: 'Pheochromocytoma and Paraganglioma',
		prad: 'Prostate adenocarcinoma', 
		read: 'Rectum adenocarcinoma',
		sarc: 'Sarcoma',
		skcm: 'Skin Cutaneous Melanoma',
		stad: 'Stomach adenocarcinoma',
		tgct: 'Testicular Germ Cell Tumors',
		thym: 'Thymoma',
		thca: 'Thyroid carcinoma',
		ucs:  'Uterine Carcinosarcoma',
		ucec: 'Uterine Corpus Endometrial Carcinoma',
		uvm:  'Uveal Melanoma' 
	}

	DATA_TYPES = {
		maseq: 		'Gene Expression (RNA-seq)',
		proteomics: 'Proteomics',
		meth27: 	'DNA Methylation (Illumina Infinium Human DNA Methylation 27)',
		meth450: 	'DNA Methylation (Illumina Infinium Human DNA Methylation 450)',
		mirnaga: 	'microRNA (Illumina GA platform)',
		mirnahiseq: 'microRNA (Illumina HiSeq platform)'
	}

	PREDICTION_TARGTS = {
		stage: 		'Tumor stage (TNM)',
		stage_t: 	'T stage',
		stage_n: 	'N stage',
		stage_m: 	'M stage' 
	}

	FEATURE_SELECTION_METHOD = {
		infog: 		'Information Gain',
		gainr: 		'Gain Ratio',
		symu: 		'Symmetrical Uncertainty',
		randf: 		'Random Forest Importance',
		_features: 	'Custom' # Prepend SESSION_ID and append TXT when submitting this
	}

	PARTITION_TYPES = ['batch', 'random']

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
end

