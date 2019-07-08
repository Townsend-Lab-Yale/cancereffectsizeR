## Script to convert new COSMIC signatures for use within deconstructSigs

# Downloaded new cosmic signatues from https://cancer.sanger.ac.uk/cosmic/signatures/ ... v3 (May 2019)

signatures_cosmic_May2019 <- read.csv(file = "~/Documents/slack_downloads/sigProfiler_SBS_signatures_2019_05_22.csv",stringsAsFactors = F)

# look at structure of the data
head(signatures_cosmic_May2019)
head(t(signatures_cosmic_May2019[,3:ncol(signatures_cosmic_May2019)]))

# align like input for deconstructSigs
signatures_cosmic_May2019_t <- t(signatures_cosmic_May2019[,3:ncol(signatures_cosmic_May2019)])

paste(signatures_cosmic_May2019$Type,signatures_cosmic_May2019$SubType)

colnames(signatures.cosmic)

# checking to make sure the data is in the same order as deconstructSigs

upstream_tri <- unlist(strsplit(signatures_cosmic_May2019$SubType,split = ""))[c(T,F,F)]

downstream_tri <- unlist(strsplit(signatures_cosmic_May2019$SubType,split = ""))[c(F,F,T)]

old_up_tri <- unlist(strsplit(colnames(signatures.cosmic),split = ""))[c(T,F,F,F,F,F,F)]
old_down_tri <- unlist(strsplit(colnames(signatures.cosmic),split = ""))[c(F,F,F,F,F,F,T)]

upstream_tri == old_up_tri

downstream_tri == old_down_tri

# all match!

colnames(signatures_cosmic_May2019_t) <- colnames(signatures.cosmic)

head(signatures_cosmic_May2019_t)

signatures_cosmic_May2019 <- as.data.frame(signatures_cosmic_May2019_t, stringsAsFactors=F)

# rates seem to be binned into two categories,
# 1) rates consistent in magnitude with the previous v2 release and
# 2) rates fixed to a floor several orders of
# magnitude less than v2 (several are exactly equal to 2.23e-16).
# changing these to zero, consistent with the v2 release.
signatures_cosmic_May2019[which(signatures_cosmic_May2019 < 1e-15,arr.ind = T)] <- 0



save(signatures_cosmic_May2019, file = "data/signatures_cosmic_May2019.RData")


