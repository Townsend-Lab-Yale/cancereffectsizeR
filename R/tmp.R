# 
# 
# 
# tmp = as.data.table(cesa$reference_data$gene_ranges)
# setkey(tmp, "names")
# 
# microbenchmark::microbenchmark(foverlaps(tmp[names == "CDKN2A.p14arf"], tmp[names == "CDKN2A.p16INK4a"], by.x = c("seqnames", "start", "end"), by.y = c("seqnames", "start", "end"), nomatch = NULL)[, .N > 0],
#                                GenomicRanges::intersect(cesa$reference_data$gene_ranges[cesa$reference_data$gene_ranges$names == "CDKN2A.p14arf"], 
#                          cesa$reference_data$gene_ranges[cesa$reference_data$gene_ranges$names == "CDKN2A.p16INK4a"]), times = 10)
