# Run with care, after loading correct version of cancereffectsizeR
setwd(system.file("tests/generate_test_data", package = "cancereffectsizeR"))
source("generate_tiny_cesa.R")
# source("generate_refcds_test.R") # [ Not typically run because it downloads tons of data from Ensembl]
source("generate_luad_cesa.R")
source("generate_luad_cesa_multi.R")
source("generate_trinuc_data.R")
