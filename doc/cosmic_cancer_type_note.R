## ----eval = F-----------------------------------------------------------------
#    suggest_cosmic_v3_signatures_to_remove(cancer_type = "BRCA", treatment_naive = TRUE)
#    suggest_cosmic_v3_signatures_to_remove(cancer_type = "Kidney-RCC")

## ---- echo=FALSE, warning=FALSE, message=FALSE--------------------------------
  data_source = paste0(system.file("extdata", package = "cancereffectsizeR"), '/pcawg_tcga_cancer_types.txt')
  cancer_type = data.table::fread(data_source)
  cancer_type[is.na(cancer_type)] = "(none)"
  formattable::formattable(cancer_type[, .(PCAWG, Applicable_TCGA, Number_of_tumors, Description)])

