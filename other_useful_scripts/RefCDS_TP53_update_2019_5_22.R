###
# building new RefCDS_TP53splice data


# old RefCDS
rm(list=ls());load("data/RefCDS_TP53splice.RData")
gr_genes_oursold <- gr_genes
RefCDS_oursold <- RefCDS
names_old_RefCDS <- sapply(X = RefCDS, function(x) x$gene_name)


# new RefCDS
load("../local_work/dndscv/dndscv-master/data/refcds_hg19.rda")
names_new_RefCDS <- sapply(X = RefCDS, function(x) x$gene_name)


# RefCDS_oursold[which(names_old_RefCDS=="TP53")]
# replace TP53 entree with the one that has the splice site correctly assigned
RefCDS[which(names_new_RefCDS=="TP53")] <- RefCDS_oursold[which(names_old_RefCDS=="TP53")]

# resave
save(RefCDS,
     gr_genes,
     file = "data/RefCDS_TP53splice.RData")





