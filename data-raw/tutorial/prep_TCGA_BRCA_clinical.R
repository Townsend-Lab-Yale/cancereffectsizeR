library(data.table)
tcga_clinical_file = system.file("extdata/brca_tcga_clinical_data_via_gdac_cbioportal.tsv", package = 'cancereffectsizeR')

tcga_clinical = fread(tcga_clinical_file)
tcga_clinical[`HER2 fish status` == "Positive" | `IHC-HER2` == 'Positive', HER2 := 'HER2+']
tcga_clinical[is.na(HER2) & (`HER2 fish status` == 'Negative' | `IHC-HER2` == 'Negative'), HER2 := 'HER2-']

tcga_clinical[`ER Status By IHC` == 'Positive', ER := 'ER+']
tcga_clinical[`ER Status By IHC` == 'Negative', ER := 'ER-']

tcga_clinical[`PR status by ihc` == 'Positive', PR := 'PR+']
tcga_clinical[`PR status by ihc` == 'Negative', PR := 'PR-']

# HR- if both PR and ER have been tested and are negative
tcga_clinical[PR == 'PR-' & ER == 'ER-', HR := 'HR-']

# HR+ if either PR+ or ER+
tcga_clinical[PR == 'PR+' | ER == 'ER+', HR := 'HR+']


tcga_clinical[! is.na(HR) & ! is.na(HER2), receptor_status := paste(HR, HER2, sep = '/')]
tcga_clinical[receptor_status == 'HR-/HER2-', receptor_status := 'TNBC']


tcga_clinical[`American Joint Committee on Cancer Metastasis Stage Code` == 'M0', pM := 'M0']
tcga_clinical[`American Joint Committee on Cancer Metastasis Stage Code` == 'M1', pM := 'M1']


tcga_clinical = unique(tcga_clinical[, .(patient_id = `Patient ID`, pM, receptor_status)])

output_file = paste0(system.file(package = "cancereffectsizeR"), '/tutorial/TCGA_BRCA_clinical.txt')
fwrite(tcga_clinical, output_file, sep = "\t", na = 'NA')

