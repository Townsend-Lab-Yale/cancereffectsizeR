###
# rebuilding PCAs with CDKN2A variants

pca_files <- dir(path = "data",pattern = "pca",full.names = T)

# add two rows with new protein names
# p.s. I know what the "a" in "pca" stands for, this is an analysis of that pca!
for(i in 1:length(pca_files)){
  pca_analysis <- get(load(pca_files[i]))
  pca_analysis$rotation <- rbind(pca_analysis$rotation,pca_analysis$rotation["CDKN2A",],pca_analysis$rotation["CDKN2A",])
  rownames(pca_analysis$rotation)[(nrow(pca_analysis$rotation)-1):(nrow(pca_analysis$rotation))] <- c("CDKN2A.p14arf","CDKN2A.p16INK4a")
  save(pca_analysis,file = pca_files[i])
}

# fixing pca names

pca_files <- dir(path = "data",pattern = "pca",full.names = T)
for(i in 1:length(pca_files)){
  pca_analysis <- get(load(pca_files[i]))
  pca_analysis$rotation <- pca_analysis$rotation[-which(rownames(pca_analysis$rotation)=="CDKN2A"),]
  new_name <- strsplit(strsplit(pca_files[i],split="/")[[1]],split = "\\.")[[2]][1]
  assign(new_name, pca_analysis)
  save(list=new_name,file = pca_files[i])
}
