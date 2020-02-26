prev_dir = setwd(system.file("tests/test_data/", package = "cancereffectsizeR"))
tiny_ak = load_maf(cesa = CESAnalysis(genome="hg19"), maf = "tiny.hg19.maf.txt", sample_col = "sample_id")
saveRDS(tiny_ak, "tiny_hg19_maf_loaded.rds")
setwd(prev_dir)
