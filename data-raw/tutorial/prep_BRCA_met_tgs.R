# Download from cBioPortal on 11/11/21
load("~/cancereffectsizeR/")
maf_file = 'data_mutations.txt'
chain_file = 'hg19ToHg38.over.chain'

# Load data, remove patients with more than one sample, lift to hg38, filter germline variants, save
# minimal necessary MAF data for cancereffectsizeR
maf = fread(maf_file)
maf[, patient_id := gsub('-T.*', '', Tumor_Sample_Barcode)]
multisample_patients = maf[, .(uniqueN(Tumor_Sample_Barcode) > 1), by = "patient_id"][V1 == T, patient_id]
maf = maf[! multisample_patients, on = 'patient_id']
maf = preload_maf(maf, refset = ces.refset.hg38, chain_file = chain_file, sample_col = 'patient_id')
maf = maf[germline_variant_site == FALSE]
maf = maf[, .(Unique_Patient_Identifier, Chromosome, Start_Position, Reference_Allele, Tumor_Allele)]
out_file = paste0(system.file("tutorial", package = "cancereffectsizeR") , '/metastatic_breast_2021_hg38.maf')
fwrite(maf, out_file, sep = "\t")
