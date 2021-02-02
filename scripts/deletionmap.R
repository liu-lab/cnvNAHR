# setwd() to working directory
del.all <- read.table("NAHR_GRCh38.withgenes.txt", sep="\t", stringsAsFactors = F, header = T)
del.all.nochr17alt <- subset(del.all, seqnames != "chr17_GL000258v2_alt")

idx.17q23 <- del.all[which(del.all$seqnames=="chr17_GL000258v2_alt" & del.all$end == 1355937),]
idx.17q23$start <- 45607798 
idx.17q23$end <- 46210001 
idx.17q23$seqnames <-"chr17"

del.all.modified <- rbind(del.all.nochr17alt, idx.17q23)
del.all.modified.curated <- subset(del.all.modified, freqCMA!="NA")
