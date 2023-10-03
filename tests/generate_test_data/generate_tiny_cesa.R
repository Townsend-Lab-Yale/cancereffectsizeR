prev_dir = setwd(system.file("tests/test_data/", package = "cancereffectsizeR"))

# 2 SNVs will trigger a warning for being reference mismatches, and 2 for being chrMT
maf = preload_maf('tiny.hg38.maf.txt', refset = 'ces.refset.hg38')
fwrite(maf, 'tiny_preloaded.txt', sep = "\t", na = 'NA')

tiny_ak = load_maf(cesa = CESAnalysis(refset ="ces.refset.hg38"), maf = maf)
save_cesa(tiny_ak, "tiny_hg38_maf_loaded.rds")

# Also doing an hg19 version for a few tests
maf = preload_maf('tiny.hg19.maf.txt', refset = 'ces.refset.hg19')
tiny_ak = load_maf(cesa = CESAnalysis(refset = 'ces.refset.hg19'), maf = maf)
save_cesa(tiny_ak, 'tiny_hg19_maf_loaded.rds')

setwd(prev_dir)

