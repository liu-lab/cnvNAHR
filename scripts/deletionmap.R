# setwd() to working directory
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install("karyoploteR")

library(karyoploteR)
library(IRanges)
library(GenomicRanges)

del.all <- read.table("NAHR_GRCh38.withgenes.txt", sep="\t", stringsAsFactors = F, header = T)
del.all.nochr17alt <- subset(del.all, seqnames != "chr17_GL000258v2_alt")

idx.17q23 <- del.all[which(del.all$seqnames=="chr17_GL000258v2_alt" & del.all$end == 1355937),]
idx.17q23$start <- 45607798 
idx.17q23$end <- 46210001 
idx.17q23$seqnames <-"chr17"

del.all.modified <- rbind(del.all.nochr17alt, idx.17q23)
del.all.modified.curated <- subset(del.all.modified, freqCMA!="NA")
del.all.modified.gr <- with(del.all.modified, GRanges(del.all.modified$seqnames, IRanges(start=as.numeric(del.all.modified$start), end=as.numeric(del.all.modified$end)))) #hg38
del.all.modified.curated.gr <- with(del.all.modified.curated, GRanges(del.all.modified.curated$seqnames, IRanges(start=as.numeric(del.all.modified.curated$start), end=as.numeric(del.all.modified.curated$end)))) #hg38

SegDup.hg38 <- read.table("C:/Users/byuan/OneDrive - SCH/Work/Manuscript_ongoing/NAHR/segdup38.txt", sep = "\t", header=F, stringsAsFactors = FALSE)
SegDup.hg38.gr <- with(SegDup.hg38, GRanges(SegDup.hg38$V2, IRanges(start=as.numeric(SegDup.hg38$V3), end=as.numeric(SegDup.hg38$V4)))) #hg38

#all in one plot
kp <- plotKaryotype(plot.type=2)
kpPlotDensity(kp, data=SegDup.hg38.gr)
kpPlotRegions(kp, data=del.all.modified.gr, data.panel=2,col="#008000")

kp <- plotKaryotype(plot.type=2)
kpPlotDensity(kp, data=SegDup.hg38.gr)
kpPlotRegions(kp, data=del.all.modified.curated.gr, data.panel=2,col="#008000")

#split
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = c("chr1","chr2","chr3","chr4","chr5", "chr6","chr7","chr8","chr9","chr10","chr11","chr12"))
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = c("chr13","chr14","chr15","chr16","chr17", "chr18","chr19","chr20","chr21","chr22","chrX","chrY"))
kpPlotDensity(kp, data=SegDup.hg38.gr, window.size = 1000)
kpPlotRegions(kp, data=del.all.modified.gr, data.panel=2, col="#008000")

kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = c("chr1","chr2","chr3","chr4","chr5", "chr6","chr7","chr8","chr9","chr10","chr11","chr12"))
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes = c("chr13","chr14","chr15","chr16","chr17", "chr18","chr19","chr20","chr21","chr22","chrX","chrY"))
kpPlotDensity(kp, data=SegDup.hg38.gr, window.size = 1000)
kpPlotRegions(kp, data=del.all.modified.curated.gr, data.panel=2, col="#008000")
