# Run from dev directory
library(devtools)
library(data.table)
load_all()


if(! 'cancereffectsizeR' %in% devtools::dev_packages()) {
  stop('Run from development version of package (i.e., load with devtools).')
}

brca_file = system.file('brca.maf.gz', package = 'cancereffectsizeR')
if (! file.exists(brca_file)) {
  get_TCGA_project_MAF(project = 'BRCA', filename = brca_file)
}

brca = fread(brca_file)
brca = brca[, .(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, 
                Tumor_Sample_Barcode, Unique_Patient_Identifier, source_file_id)]
set.seed(2001)
random_patients = sample(unique(brca$Unique_Patient_Identifier), size = 200, replace = FALSE)
random_patients_tiny = sample(unique(brca$Unique_Patient_Identifier), size = 5, replace = FALSE)

tiny_maf = brca[random_patients_tiny, on = 'Unique_Patient_Identifier']
brca = brca[random_patients, on = 'Unique_Patient_Identifier']

fwrite(brca, 'tests/test_data/brca_subset_hg38.maf.gz', sep = "\t")


# Make a very small MAF for additional testing
tiny_maf = tiny_maf[, .(Unique_Patient_Identifier,	Chromosome,	Start_Position,
                        Reference_Allele, Tumor_Seq_Allele2)]
more_records = fread('tests/test_data/invalid_or_unusual_maf_records.maf')

tiny_maf = rbind(tiny_maf, more_records)
fwrite(tiny_maf, 'tests/test_data/tiny.hg38.maf.txt', sep = "\t")

