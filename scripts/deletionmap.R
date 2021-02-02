# setwd() to working directory
del.all <- read.table("NAHR_GRCh38.withgenes.txt", sep="\t", stringsAsFactors = F, header = T)
del.curated <- subset(del.all, freqCMA != "NA", select = c(seqnames, start, end))
