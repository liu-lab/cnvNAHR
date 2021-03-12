if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install("karyoploteR")

library(karyoploteR)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)

del.all <- read.table("./output/NAHR_GRCh38.withgenes.tsv", sep="\t", stringsAsFactors = F, header = T)
del.all.nochr17alt <- subset(del.all, seqnames != "chr17_GL000258v2_alt")

idx.17q23 <- del.all[which(del.all$seqnames=="chr17_GL000258v2_alt" & del.all$end == 1355937),]
idx.17q23$start <- 45607798 
idx.17q23$end <- 46210001 
idx.17q23$seqnames <-"chr17"

del.all.modified <- rbind(del.all.nochr17alt, idx.17q23)
del.all.modified.curated <- subset(del.all.modified, freqCMA!="NA")
del.all.modified.gr <- with(del.all.modified, GRanges(seqnames, IRanges(start=as.numeric(start), end=as.numeric(end)))) #hg38
del.all.modified.curated.gr <- with(del.all.modified.curated, GRanges(seqnames, IRanges(start=as.numeric(start), end=as.numeric(end)))) #hg38

mySession = browserSession("UCSC")
genome(mySession) <- "hg38"
segdup.raw <- getTable(ucscTableQuery(mySession, track="genomicSuperDups"))

SegDup.hg38 <- read.table("C:/Users/byuan/OneDrive - SCH/Work/Manuscript_ongoing/NAHR/segdup38.txt", sep = "\t", header=F, stringsAsFactors = FALSE)
SegDup.hg38.gr <- with(segdup.raw, GRanges(chrom, IRanges(start=as.numeric(chromStart), end=as.numeric(chromEnd)))) #hg38

#all in one plot
pdf(file = "./output/NAHR_GRCh38.plot.pdf")
kp <- plotKaryotype(plot.type=2)
kpPlotDensity(kp, data=SegDup.hg38.gr, window.size = 1000)
kpPlotRegions(kp, data=del.all.modified.gr, data.panel=2,col="#008000", r0 = 0, r1=3.5)
dev.off()
