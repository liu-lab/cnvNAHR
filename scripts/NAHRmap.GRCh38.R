library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)
library(biomaRt)
library(magrittr)
library(tidyr)
library(VariantAnnotation)

mySession = browserSession("UCSC")
genome(mySession) <- "hg38"

ensembl <- useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", 
                      dataset="hsapiens_gene_ensembl",
                      mirror = "uswest")

attributes <- c("ensembl_gene_id","chromosome_name","start_position","end_position","hgnc_symbol", "hgnc_id", "ensembl_transcript_id", 
                "entrezgene_id", "gene_biotype") 

chr.names <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
               "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chr17_GL000258v2_alt")

## download tables from UCSC
segdup.raw <- getTable(ucscTableQuery(mySession, track="genomicSuperDups"))
cen.raw <- getTable(ucscTableQuery(mySession, track="centromeres"))
gap.raw <- getTable(ucscTableQuery(mySession, track="gap"))

gene.bm <- getBM(attributes=attributes, values="*", mart=ensembl) %>% 
  mutate(chromosome_name = paste0("chr", chromosome_name), hgnc_id = as.integer(gsub("HGNC:", "", hgnc_id)))
gene.bm.mod <- read_tsv("./input/chr17_GL000258v2_alt_genes.txt") %>% bind_rows(gene.bm) %>% filter(hgnc_id !="")

gr.gene.all <- with(gene.bm.mod, GRanges( seqnames = chromosome_name,  IRanges(start=start_position, end=end_position), 
                                          hgnc_symbol = hgnc_symbol, hgnc_id = hgnc_id, 
                                          ensembl_transcript_id=ensembl_transcript_id, 
                                          contig_hgnc=paste0(chromosome_name, "_", hgnc_id))) # created this key to separate groups of genes on differnt contigs

gr.gene.grouped <- gr.gene.all %>% splitAsList(gr.gene.all$contig_hgnc) 
gr.gene.grouped.maxidx <- gr.gene.grouped %>% width %>% which.max
#save only longest transcript for each gene
gr.gene.all <- do.call(c, lapply(1:length(gr.gene.grouped.maxidx), function(x)
{
  gr.gene.grouped[[names(gr.gene.grouped.maxidx[x])]][unname(gr.gene.grouped.maxidx[x])]
}))

coding.gene.id <- gene.bm.mod %>% filter(gene_biotype=="protein_coding") %>% .$hgnc_id
gr.gene.coding <- gr.gene.all[gr.gene.all$hgnc_id %in% coding.gene.id]




# load input data
Constraint <- read_tsv("./input/gnomad.v2.1.1.lof_metrics.by_gene.txt") # can be downloaded from the gnomAD website
DDD <- read_csv("./input/DDG2P_5_1_2021.csv.gz") %>% mutate(hgnc_id=paste0("HGNC:", `hgnc id`)) # can be downloaded from the DDD website
Mim2gene <- read_tsv("./input/mim2gene.txt", skip=4)  # 01/04/2021 version from OMIM
Genemap2 <- read_tsv("./input/genemap2.txt", skip=3)  # 01/04/2021 version from OMIM
# extract relevant information from the OMIM genemap2 file
Genemap2_2 <- Genemap2 %>%
  mutate(Phenotypes = as.character(Phenotypes)) %>%
  dplyr::select(`MIM Number`,`Gene Symbols`,`Approved Symbol`, `Entrez Gene ID`, `Ensembl Gene ID`,Phenotypes) %>%
  filter(Phenotypes!="" & `Entrez Gene ID` !="") %>%
  mutate(Phenotypes = strsplit(Phenotypes, "; ")) %>%
  unnest(Phenotypes) %>%
  mutate(MIM_dz_name = str_extract_all(Phenotypes, '^[^\\{\\[].*(?=\\s\\d{6})') %>% str_c(., sep=";"),
         MIM_dz_name_ass = str_extract_all(Phenotypes, '(^\\{).*(?=\\s\\d{6})') %>% str_c(., sep=";"),
         MIM_dz_id =  str_extract_all(Phenotypes, '\\d{6}' ),
         MIM_dz_anno = str_extract_all(Phenotypes, '\\(\\d{1}\\)')  %>% as.matrix %>% as.character,
         Inheritance = str_extract_all(Phenotypes, '(?<=\\(\\d{1}\\),).*') %>% as.matrix %>% as.character) %>%
  filter(MIM_dz_anno=="(3)" & MIM_dz_name != "character(0)")  # exclude [] and {}, only select (3)

Biallelic.OMIM <- c(" Autosomal recessive"," Digenic recessive" ,
                    " Autosomal recessive, Digenic recessive",
                    " Autosomal recessive, Digenic dominant",
                    " Autosomal recessive, Mitochondrial",                       
                    " Autosomal recessive, ?Autosomal dominant"  ,             
                    " Autosomal recessive, Multifactorial",    
                    " Autosomal recessive, Somatic mosaicism",  
                    " Pseudoautosomal recessive" ,
                    "Autosomal recessive, Autosomal dominant", 
                    " Autosomal recessive, Multifactorial, Autosomal dominant",   
                    " Autosomal recessive, Digenic dominant, Autosomal dominant", 
                    " Autosomal recessive, Autosomal dominant, Isolated cases" , 
                    " Autosomal recessive, Autosomal dominant, Somatic mutation")

## convert all gene names to hgnc gene id

Biallelic.OMIM.gene.value <- Genemap2_2 %>% filter(Inheritance %in% Biallelic.OMIM) %>% .$`Entrez Gene ID` %>% unique
Biallelic.OMIM.gene.bm <- getBM(attributes=attributes, filter="entrezgene_id", values=Biallelic.OMIM.gene.value, mart=ensembl)
Biallelic.OMIM.gene <- Biallelic.OMIM.gene.bm$hgnc_id %>% unique

MonoAndBi.OMIM.gene.bm <- getBM(attributes=attributes, filter="entrezgene_id", values=unique(Genemap2_2$`Entrez Gene ID`), mart=ensembl)
MonoAndBi.OMIM.gene  <- MonoAndBi.OMIM.gene.bm$hgnc_id %>% unique


###### DDD  #########

Biallelic.DDD <- c("biallelic", "biallelic,monoallelic","biallelic,uncertain")

MonoAndBi.DDD.gene <-  DDD$hgnc_id %>% unique

Biallelic.DDD.gene <- DDD %>% filter(`allelic requirement` %in% Biallelic.DDD) %>% .$hgnc_id %>% unique

###### ClinGen #######

clingen <- read_csv("./input/Clingen-Gene-Disease-Summary-2021-01-31.csv", skip = 6, 
                col_names = c("GENE_SYMBOL",	"hgnc_id",	"DISEASE_LABEL",	"DISEASE_MONDO_ID",	"MOI",	"SOP",	"CLASSIFICATION", 
                              "ONLINE_REPORT", "CLASSIFICATION_DATE",	"GCEP"))
Biallelic.clingen.gene <- clingen %>% filter(clingen$MOI=="AR") %>% .$hgnc_id
MonoandBiallelic.clingen.gene <- clingen %>% .$hgnc_id


######## pRec, pLI #######
genes.pRec.value <- Constraint %>% filter(pRec > 0.99) %>% .$transcript
genes.pRec.bm <- getBM(attributes=attributes, filter="ensembl_transcript_id", values=genes.pRec.value, mart=ensembl)
genes.pRec <- genes.pRec.bm$hgnc_id %>% unique

genes.pLI.value <- Constraint %>% filter(pLI > 0.99) %>% .$transcript
genes.pLI.bm <- getBM(attributes=attributes, filter="ensembl_transcript_id", values=genes.pLI.value, mart=ensembl)
genes.pLI <- genes.pLI.bm$hgnc_id %>% unique

#Biallelic.gene.omim.ddd <- unique(c(Biallelic.OMIM.gene, Biallelic.DDD.gene))
Biallelic.gene.omim.ddd.clingen <- setdiff(unique(c(Biallelic.OMIM.gene, Biallelic.DDD.gene, Biallelic.clingen.gene)), "")
MonoAndBi.gene.omim.ddd.clingen <- setdiff(unique(c(MonoAndBi.OMIM.gene, MonoAndBi.DDD.gene, MonoandBiallelic.clingen.gene)), "")
write_lines(Biallelic.gene.omim.ddd.clingen, "./output/biallelicgenes.txt")
write_lines(MonoAndBi.gene.omim.ddd.clingen, "./output/biallelicmonoallelicgenes.txt")

hgncid.to.gene.gr <- function(ids)
{
  ids1 <- gsub("HGNC:", "", ids)
  gr.gene.all[gr.gene.all$hgnc_id %in% ids1]
}

ens.transcript.to.gene.gr <- function(ids)
{
  hgncids <- getBM(attributes=attributes, filter="ensembl_transcript_id", values=ids, mart=ensembl) %>% .$hgnc_id
  hgncid.to.gene.gr(hgncids)
}

gr.genes.pLI999 <- ens.transcript.to.gene.gr(Constraint %>% filter(pLI > 0.999) %>% .$transcript)
gr.genes.pLI99 <- ens.transcript.to.gene.gr(Constraint %>% filter(pLI > 0.99 & pLI <=0.999)  %>% .$transcript)
gr.genes.pLI90 <- ens.transcript.to.gene.gr(Constraint %>% filter(pLI > 0.9 & pLI <= 0.99) %>% .$transcript)


gr.genes.biallelic <- hgncid.to.gene.gr(Biallelic.gene.omim.ddd.clingen)
gr.genes.monobiallelic <- hgncid.to.gene.gr(MonoAndBi.gene.omim.ddd.clingen)
gr.genes.pRec.bm <- ens.transcript.to.gene.gr(Constraint %>% filter(pRec > 0.99) %>% .$transcript)

## prepare segdup file for 17q21.31
# cd soft/
# wget "http://www.littlest.co.uk/software/pub/bioinf/freeold/miropeats_2.02.tgz"
# wget "http://www.littlest.co.uk/software/pub/bioinf/freeold/icatools_2.5.tgz"
# gunzip miropeats_2.02.tgz
# gunzip icatools_2.5.tgz
# tar -xvf miropeats_2.02.tar
# tar -xvf icatools_2.5.tar
# cd icatools_2.5/
# PATH=$PATH:/home/pengfei/soft/icatools_2.5
# cd /mnt/pure/rnd/pengfei/NAHR/17q21.31
# /home/pengfei/soft/miropeats_2.02/miropeats -s 200 -onlyinter -seq /mnt/pure/rnd/pengfei/temp/chr17_GL000258v2_alt_1-1821992.fa -seq /mnt/pure/rnd/pengfei/temp/chr17_45309498-46836265.fa
# /home/pengfei/soft/miropeats_2.02/miropeats -s 1000 -onlyinter -seq /mnt/pure/rnd/pengfei/temp/chr17_GL000258v2_alt_1-1821992.fa -seq /mnt/pure/rnd/pengfei/temp/chr17_45309498-46836265.fa
# /home/pengfei/soft/miropeats_2.02/miropeats -s 1000 -onlyintra -seq /mnt/pure/rnd/pengfei/temp/chr17_GL000258v2_alt_1-1821992.fa
# /home/pengfei/soft/miropeats_2.02/miropeats -s 1000 -onlyinter -seq /mnt/pure/rnd/pengfei/temp/chr17_GL000258v2_alt_450000-80000.fa -seq /mnt/pure/rnd/pengfei/temp/chr17_GL000258v2_alt_1250000-1460000.fa

segdup.chr17_GL000258v2_alt <- read.table("./input/chr17_GL000258v2_alt_segdup.txt", stringsAsFactors = F, header = T)

segdup.raw <- rbind(segdup.raw, segdup.chr17_GL000258v2_alt)

#define the chromosomes used 
## preprocess gaps
gap <- gap.raw[which(gap.raw$chrom %in% chr.names),]
gap$chrom <- factor(gap$chrom)
gap$chrom <- as.character(gap$chrom)
# make ranges that are intervals in between segdups
gr.gap <- with(gap, GRanges( seqnames = chrom,  ranges = IRanges(chromStart, end = chromEnd)))
# heterochromatin
gap.heterochromatin <- gap[gap$type=="heterochromatin",]
gap.heterochromatin$chrom <- factor(gap.heterochromatin$chrom)
gap.heterochromatin$chrom <- as.character(gap.heterochromatin$chrom)
gr.gap.hetochromatin <- with(gap.heterochromatin, GRanges( seqnames = chrom,  ranges = IRanges(chromStart, end = chromEnd)))
## preprocess centromere
cen.merged <- do.call(rbind, lapply(chr.names[1:24], function(x)
{
  Chr.sel <- x
  Start.sel <- min( cen.raw[which(cen.raw$chrom==Chr.sel),"chromStart"] )
  End.sel <- max( cen.raw[which(cen.raw[,2]==Chr.sel),"chromEnd"] )
  data.frame(chrom=Chr.sel, chromStart=as.numeric(Start.sel), chromEnd=as.numeric(End.sel))
}))
# make granges for centromeres
gr.cen.merged <- with(cen.merged, GRanges( seqnames = chrom,  ranges = IRanges(chromStart, end = chromEnd)))

gr.cen.merged <- c(gr.cen.merged, gr.gap.hetochromatin)
gr.cen.merged <- reduce(gr.cen.merged)

## preprocess segdup table
#only in chr1-22, X, Y, and chr17_GL000258v2_alt
Segdup0 <- segdup.raw[which(segdup.raw$chrom %in% chr.names),]
Segdup0$chrom <- factor(Segdup0$chrom)
Segdup0$chrom <- as.character(Segdup0$chrom)

# same chr, same orientation
Segdup1 <- Segdup0[which( (Segdup0$strand == "+") & (Segdup0$chrom == Segdup0$otherChrom) ), ]

# annotate CNVs
Segdup1$key <- paste0(Segdup1$chrom, ":", Segdup1$chromStart, "-", Segdup1$chromEnd)
Segdup1$otherkey <- paste0(Segdup1$otherChrom, ":", Segdup1$otherStart,  "-", Segdup1$otherEnd)
Segdup1$cnvStart <- sapply(1:nrow(Segdup1), function(x)
{
  round(mean(c(Segdup1$chromStart[x], Segdup1$chromEnd[x])))
})

Segdup1$cnvEnd <- sapply(1:nrow(Segdup1), function(x)
{
  round(mean(c(Segdup1$otherStart[x], Segdup1$otherEnd[x])))
})

# remove duplicates, Start < Start.other
Segdup2 <- Segdup1[which(Segdup1$chromStart < Segdup1$otherStart),]

gr.segdup1 <- with(Segdup1, GRanges( seqnames = chrom,  ranges = IRanges(chromStart, end = chromEnd)))
## make granges for NAHR mediated CNV regions. Note this CNV regions, Not segdup regions. 
gr.NAHR <- with(Segdup2, GRanges( seqnames = chrom,  ranges = IRanges(cnvStart, end = cnvEnd)))
# NAHR CNV regions do not span centromere
if(!identical(queryHits(findOverlaps( gr.NAHR, gr.cen.merged)), integer(0)))
{
  gr.NAHR <- gr.NAHR[-queryHits(findOverlaps( gr.NAHR, gr.cen.merged))]
} 
# NAHR CNV regions contain coding genes
gr.NAHR <- gr.NAHR[unique(queryHits(findOverlaps(gr.NAHR, gr.gene.coding)))]

## make granges for segdups. These are the raw segdup data with many overlapping ones. The next few lines of code aim to merge them. 
gr.sd <- with(Segdup0, GRanges( seqnames = chrom,  ranges = IRanges(chromStart, end = chromEnd)))
# add gaps into segdups to be merged
gr.sd <- c(gr.sd, gr.gap)
# first simple merge so that they are reduced to nonoverlapping intervals
gr.sd <- reduce(gr.sd)
# sorting to prepare for the next merging 
seqlevels(gr.sd) <- sort(seqlevels(gr.sd))
gr.sd<- sort(gr.sd)
# break them into chromosomes. make ranges for intervals in between the collapsed segdups. The idea is that if the interval is small, then the elements before and after the interval will be assigned to the same group to be merged. 
gr.sd.interval <- do.call(c,lapply(chr.names, function(x)
{
  gr.sd.chr <- gr.sd[seqnames(gr.sd)==x]
  gr.sd.chrone <- GRanges(seqnames = x, IRanges(min(start(gr.sd.chr)), max(end(gr.sd.chr))))
  bksd <- disjoin(c(gr.sd.chrone, gr.sd.chr))
  bksd[countOverlaps(bksd, gr.sd.chr) == 0]
}))
# sort again
seqlevels(gr.sd.interval) <- sort(seqlevels(gr.sd.interval))
gr.sd.interval<- sort(gr.sd.interval)
# This is a logic list of YES or NO whether element i from the vector of reduced segdups can be merged to the next element i+1
# the logic is based on 1) the interval is < 1kb, or 2) the interval is between 1 to 50 kb but there is no genes in it
# the output is a list splitted by chromosomes
sd.goup.logic.list <- lapply(split(gr.sd.interval, seqnames(gr.sd.interval)), function(x){
  sd.goup.logic <- width(x)<1e3 |
    (width(x)>1e3 & width(x)<5e4 & !(1:length(x) %in% unique(queryHits(findOverlaps(x, gr.gene.coding))) ))
  sd.goup.logic
})

# based on the logic assign group names to each element. Elements with the same group name can be merged. Assign them to the segdup reduced Granges
sd.goup.list <- lapply(sd.goup.logic.list, function(xx){
  i=1
  sd.goup <- c(i, sapply(xx, function(x){
    if(x)
    {i} else {i<<-i+1}
  }))
  sd.goup
})
sd.goup.list <- lapply(names(sd.goup.list), function(x){
  paste0(x, "_", sd.goup.list[[x]])
})

gr.sd$sdGroup = unlist(sd.goup.list)

# merge the sedgups based on the assigned group number
gr.sd.merged <- do.call(c,lapply(unique(gr.sd$sdGroup), function(x)
{
  sd.group.sel <- gr.sd[gr.sd$sdGroup==x]
  GRanges(seqnames = seqnames(sd.group.sel)[1], IRanges(min(start(sd.group.sel)), max(end(sd.group.sel))), strand="+")
}))



# get indexes of self overlap within gr.NAHR
overlap.gr.NAHR <- findOverlaps(gr.NAHR, drop.self=TRUE, drop.redundant=TRUE)

### create temporary list with pairwise overlaps; 
### keep removing rows from temporary list when pairs have been assigned
### to an overlapping group of regions
temp <- as.data.frame(overlap.gr.NAHR)
list.overlaps <- list()
while(nrow(temp) > 0){
  
  qh <- temp$queryHits[1]
  overlaps <- c(qh,temp$subjectHits[temp$queryHits == qh])
  temp <- temp[temp$queryHits != qh,]
  
  ### find all the other regions that overlap when the 'subjectHit' 
  ### becomes the 'queryHit'
  i <- 2
  while(i <= length(overlaps)){
    newHits <- temp$subjectHits[temp$queryHits == overlaps[i]]
    overlaps <- c(overlaps, newHits[!(newHits %in% overlaps)])
    temp <- temp[temp$queryHits != overlaps[i],]
    i <- i + 1
  }
  
  ### store list of overlapping regions when no new paired regions are found
  list.overlaps <- append(list.overlaps, list(overlaps))
  #print(nrow(temp))
}

nonoverlaps <- setdiff(1:length(gr.NAHR), unlist(list.overlaps) %>% unique) %>% sapply(list)
list.overlaps <- append(list.overlaps, nonoverlaps)

## merging NAHR CNVs
gr.CNVmerged <- do.call(c,lapply(1:length(list.overlaps), function(i)  ## for each cluster of CNVs overlapping, do a loop
{
  idx.group <- list.overlaps[[i]]
  
  # logic to merge: if two CNVs both have left and right breakpoints fall into the same pair of merged segdup, the two CNVs can be merged
  # make keys that identify which CNVs can be merged
  # although different CNVs here may be mediated by different segdup pairs, they can be harmonized to the merged segdup regions to form different sub-groups
  keys <- sapply(idx.group, function(x){
    gr.left<- gr.NAHR[x]
    end(gr.left)=start(gr.left)+1
    key.left <- subsetByOverlaps(gr.sd.merged, gr.left) %>% paste0
    
    gr.right<- gr.NAHR[x]
    start(gr.right)=end(gr.right)-1
    key.right <- subsetByOverlaps(gr.sd.merged, gr.right) %>% paste0
    paste(key.left, "__", key.right)  
  })
  
  # merge CNV by key
  gr.CNVmerged.group <- do.call(c,lapply(unique(keys), function(x)
  {
    
    idx.tomerge <- which(keys==x)
    #gr.tomerge <- gr.NAHR[idx.group][idx.tomerge]
    sd.left <-strsplit(x, "__")[[1]][1]
    chr.merged <- strsplit(sd.left, ":")[[1]][1]
    # granges for the left breakpoint region. This region is the reduced region. 
    gr.sd.left <- GRanges(seqnames = chr.merged, 
                          IRanges(start=as.numeric(strsplit(strsplit(sd.left, ":")[[1]][2], "-")[[1]][1]), 
                                  end=as.numeric(strsplit(strsplit(sd.left, ":")[[1]][2], "-")[[1]][2])), strand = "+")
    sd.right <-strsplit(x, "__")[[1]][2] %>% gsub(" ", "", .)
    # granges for the right breakpoint region. This region is the reduced region. 
    gr.sd.right <- GRanges(seqnames = chr.merged, 
                           IRanges(start=as.numeric(strsplit(strsplit(sd.right, ":")[[1]][2], "-")[[1]][1]), 
                                   end=as.numeric(strsplit(strsplit(sd.right, ":")[[1]][2], "-")[[1]][2])), strand = "+")
    
    # identify all the non-reduced segdups in the left breakpoint region;
    # find their pairs based on the raw segdup map;
    # check if these pairs fall in the right breakpoint region;
    # if they are, save this pair as an element of flanking repeats 
    key.left.individual.sds <- paste0(subsetByOverlaps(gr.segdup1, gr.sd.left))
    key.left.individual.sds.pair <- Segdup1$otherkey[which(Segdup1$key %in% key.left.individual.sds)]
    key.right.individual.sds <- paste0(subsetByOverlaps(gr.segdup1, gr.sd.right))
    idx.qh.with.pair <- which(key.left.individual.sds.pair %in% key.right.individual.sds)
    idx.segdup1.left.tomerge <- queryHits(findOverlaps(gr.segdup1, gr.sd.left))[idx.qh.with.pair]
    
    segdup.tomerge <- Segdup1[idx.segdup1.left.tomerge, ]
    segdup.tomerge <- segdup.tomerge[segdup.tomerge$chromStart < segdup.tomerge$otherStart, ]
    length.score <- sapply(1:nrow(segdup.tomerge), function(x)
    {
      if(segdup.tomerge$otherSize[x]>100e3)
      {250}   else
      {
        round(segdup.tomerge$otherSize[x]/100e3*250)
      }
    })
    
    identity.score<- round((segdup.tomerge$fracMatch-0.9)*10*250)
    
    segdup.tomerge$score <- length.score + identity.score
    
    length.merged <- sum(segdup.tomerge$alignB)
    identity.merged <- sum(segdup.tomerge$matchB)/sum(segdup.tomerge$alignB)
    
    start.merged <- weighted.mean(ave(segdup.tomerge$chromStart,segdup.tomerge$chromEnd), w = segdup.tomerge$score)
    end.merged <- weighted.mean(ave(segdup.tomerge$otherStart,segdup.tomerge$otherEnd), w = segdup.tomerge$score)
    
    GRanges(seqnames=chr.merged, IRanges(start=start.merged, end=end.merged), length=length.merged, identity=identity.merged)
  }))
  
  gr.CNVmerged.group
}))

# CNV does not span centromeres
if(!identical(queryHits(findOverlaps( gr.CNVmerged, gr.cen.merged)), integer(0)))
{
  gr.CNVmerged <- gr.CNVmerged[-queryHits(findOverlaps( gr.CNVmerged, gr.cen.merged))]
} 

# CNV does not completely contained within a segdup cluster
if(!identical(queryHits(findOverlaps( gr.CNVmerged, gr.sd.merged, type = "within")), integer(0)))
{
  gr.CNVmerged <- gr.CNVmerged[-queryHits(findOverlaps(gr.CNVmerged, gr.sd.merged, type = "within"))]
}

gr.CNVmerged <- unique(gr.CNVmerged)

# NAHR CNV regions contain genes
gr.CNVmerged <- gr.CNVmerged[unique(queryHits(findOverlaps(gr.CNVmerged, gr.gene.coding)))]

CNVmerged <- as.data.frame(gr.CNVmerged)


length.score2 <- sapply(1:nrow(CNVmerged), function(x)
{
  if(CNVmerged$length[x]>200e3)
  {200}   else
  {
    round(CNVmerged$length[x]/200e3*200)
  }
})

identity.score2<- round((CNVmerged$identity-0.9)*10*200)

width.score <- sapply(1:nrow(CNVmerged), function(x)
{
  if(CNVmerged$width[x]>5e6)
  {0}   else if (CNVmerged$width[x]<500e3)
  {100} else
  {
    round((100-(CNVmerged$width[x]-500e3)/4.5e6*100))
  }
})

gene.score <- sapply(1:nrow(CNVmerged), function(x)
{
  hits90 <- queryHits(findOverlaps(gr.genes.pLI90, gr.CNVmerged[x]))
  count90 <- elementMetadata(gr.genes.pLI90)[hits90,1] %>% unique %>% length
  hits99 <- queryHits(findOverlaps(gr.genes.pLI99, gr.CNVmerged[x]))
  count99 <- elementMetadata(gr.genes.pLI99)[hits99,1] %>% unique %>% length
  hits999 <- queryHits(findOverlaps(gr.genes.pLI999, gr.CNVmerged[x]))
  count999 <- elementMetadata(gr.genes.pLI999)[hits999,1] %>% unique %>% length
  500 - 10*count90 - 30*count99 - 150*count999
})


CNVmerged$score <- length.score2 + identity.score2 + width.score + gene.score
CNVmerged.final <- CNVmerged[which(CNVmerged$score>0), ]
CNVmerged.final$info <- with(CNVmerged.final, paste0(round(length/1e3), "kb_", round(identity*100, 1), "%"))

gr.CNVmerged.final <- with(CNVmerged.final, GRanges(seqnames = seqnames, IRanges(start= start, end= end)))

overlaps2 <- findOverlaps(gr.CNVmerged.final, gr.genes.biallelic)
idx.all.nahr.recessive <- as.data.frame(overlaps2) %>% 
  cbind(gene=elementMetadata(gr.genes.biallelic)[subjectHits(overlaps2), "hgnc_symbol"]) %>%
  group_by(queryHits) %>%
  summarise(genes = paste0(unique(gene), collapse = ",")) %>%
  dplyr::select(queryHits, genes) %>%
  unique

overlaps3 <- findOverlaps(gr.CNVmerged.final, gr.genes.monobiallelic)
idx.all.nahr.monobiallelic <- as.data.frame(overlaps3) %>% 
  cbind(gene=elementMetadata(gr.genes.monobiallelic)[subjectHits(overlaps3), "hgnc_symbol"]) %>%
  group_by(queryHits) %>%
  summarise(genes = paste0(unique(gene), collapse = ",")) %>%
  dplyr::select(queryHits, genes) %>%
  unique

overlaps4 <- findOverlaps(gr.CNVmerged.final, gr.genes.pRec.bm)
idx.genes.pRec.bm <- as.data.frame(overlaps4) %>% 
  cbind(gene=elementMetadata(gr.genes.pRec.bm)[subjectHits(overlaps4), "hgnc_symbol"]) %>%
  group_by(queryHits) %>%
  summarise(genes = paste0(unique(gene), collapse = ",")) %>%
  dplyr::select(queryHits, genes) %>%
  unique

overlaps5 <- findOverlaps(gr.CNVmerged.final, gr.gene.coding)
idx.genes.bm <- as.data.frame(overlaps5) %>% 
  cbind(gene=elementMetadata(gr.gene.coding)[subjectHits(overlaps5), "hgnc_symbol"]) %>%
  group_by(queryHits) %>%
  summarise(genes = paste0(unique(gene), collapse = ",")) %>%
  dplyr::select(queryHits, genes) %>%
  unique

overlaps6 <- findOverlaps(gr.CNVmerged.final, gr.gene.all)
idx.genes.all.bm <- as.data.frame(overlaps6) %>% 
  cbind(gene=elementMetadata(gr.gene.all)[subjectHits(overlaps6), "hgnc_symbol"]) %>%
  group_by(queryHits) %>%
  summarise(genes = paste0(unique(gene), collapse = ",")) %>%
  dplyr::select(queryHits, genes) %>%
  unique

num.genes.all.bm <- as.data.frame(overlaps6) %>% 
  cbind(gene=elementMetadata(gr.gene.all)[subjectHits(overlaps6), "hgnc_symbol"]) %>%
  group_by(queryHits) %>%
  summarise(length(unique(gene)))

num.genes.coding.bm <- as.data.frame(overlaps5) %>% 
  cbind(gene=elementMetadata(gr.gene.coding)[subjectHits(overlaps5), "hgnc_symbol"]) %>%
  group_by(queryHits) %>%
  summarise(length(unique(gene)))


CNVmerged.final$biallelic <- ""
CNVmerged.final$biallelic[idx.all.nahr.recessive$queryHits] <- idx.all.nahr.recessive$genes

CNVmerged.final$omim_ddd <- ""
CNVmerged.final$omim_ddd[idx.all.nahr.monobiallelic$queryHits] <- idx.all.nahr.monobiallelic$genes

CNVmerged.final$pRec <- ""
CNVmerged.final$pRec[idx.genes.pRec.bm$queryHits] <- idx.genes.pRec.bm$genes

CNVmerged.final$allgenescoding <- ""
CNVmerged.final$allgenescoding[idx.genes.bm$queryHits] <- idx.genes.bm$genes

CNVmerged.final$allgenes <- ""
CNVmerged.final$allgenes[idx.genes.all.bm$queryHits] <- idx.genes.all.bm$genes

CNVmerged.final$numallgenes <- num.genes.all.bm$`length(unique(gene))`

CNVmerged.final$numcodinggenes <- num.genes.coding.bm$`length(unique(gene))`

CNVmerged.bed <- CNVmerged.final[, c("seqnames", "start", "end", "info", "score")]
NAHR.final.bed <- file("./output/NAHR_GRCh38.bed", open="wt")
writeLines(paste("track name=NAHR description='NAHR' useScore=1 visibility=3"), NAHR.final.bed)
write.table(unique(CNVmerged.bed), NAHR.final.bed, quote=F, sep="\t", row.names= F, col.names = F)
close(NAHR.final.bed)

### make table for NAHR and frequencies
CNVmerged.final1 <- CNVmerged.final %>% arrange(desc(score) )
Recurrentdel <- read_tsv("./input/cnv_known.freq.GRCh38.txt") 
Recurrentdel0 <- Recurrentdel %>% dplyr::select(end, name, freq, freqCMA)
CNVmerged.final1 <- left_join(unique(CNVmerged.final1), Recurrentdel0) %>% unique

write_tsv(CNVmerged.final1, "./output/NAHR_GRCh38.withgenes.tsv")

