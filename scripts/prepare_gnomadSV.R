library(dplyr)
library(readr)
library(stringr)
library(GenomicRanges)
library(VariantAnnotation)
library(rtracklayer)
library(rbcf)


## load recurrent genomic deletion coordinates. Some of these are represented in gnomadSV. They need to be eliminated because the current calculation is based on small SV in gnomad. 
Recurrentdel.hg19 <- read_tsv("./output/NAHR_GRCh37.withgenes.tsv") %>% mutate(seqnames=gsub("chr", "",seqnames)) %>% filter(freq>0)
idx.17q23.hg19 <- which(Recurrentdel.hg19$seqnames=="17_ctg5_hap1" & Recurrentdel.hg19$end == 1307210)
Recurrentdel.hg19$start[idx.17q23.hg19] <- 43685164	
Recurrentdel.hg19$end[idx.17q23.hg19] <- 44287367
Recurrentdel.hg19$seqnames[idx.17q23.hg19] <-"17"
gr.recurrentdel.hg19 <- with(Recurrentdel.hg19, GRanges(seqnames = seqnames, IRanges(start=start, end=end)))


## download gnomad SV from website

filename <- "./input.bigfiles/gnomad/gnomad_v2.1_sv.sites.vcf.gz"
filenameout <- "./input.bigfiles/gnomad/gnomad_v2.1_sv.sites.Lof.nosegdup.vcf"
fp <- bcf.open(filename,FALSE)
# error on opening (exit 0 for tests)
if(is.null(fp)) quit(save="no",status=0,runLast=FALSE)
# current variant
vc <- NULL
out <- bcf.new.writer(fp,filenameout)

## function for filtering rule: CNV cannot be within segdup
mySession = browserSession("UCSC")
genome(mySession) <- "hg19"
segdup19.raw <- getTable(ucscTableQuery(mySession, track="genomicSuperDups")) %>% mutate(chrom = gsub("chr", "", chrom))
gr.segdup19 <- with(segdup19.raw, GRanges( seqnames = chrom,  ranges = IRanges(chromStart, end = chromEnd))) %>% reduce

logic.segdup = function(vc)
{
  gr.gnomadSV.LOF <- GRanges(seqnames = variant.chrom(vc), IRanges(variant.start(vc), variant.end(vc)))
  gr.overlap <- subsetByOverlaps(gr.gnomadSV.LOF, gr.segdup19)
  ## remove gnomad cnvs that has >80% content located within any segdup region.
  if(length(gr.overlap)<1) {return(TRUE)}
  sum(width(reduce(gr.overlap)))/ width(gr.gnomadSV.LOF) < 0.8
  
}

while(!is.null(vc<-bcf.next(fp))) 
{
  if(paste(variant.alt.alleles(vc), collapse = ",") %in% c("<DEL>", "<DUP>"))
  {
    if((variant.has.attribute(vc,"PROTEIN_CODING__LOF") | 
        variant.has.attribute(vc, "PROTEIN_CODING__DUP_LOF") | 
        variant.id(vc)=="gnomAD-SV_v2.1_DEL_16_152725") &  ## this SV is incorrectly annotated by gnomAD. The annotation missed the ABCC6 gene overlap so did not label it as LOF. 
       !variant.is.filtered(vc) &
       variant.float.attribute(vc = vc, att = "POPMAX_AF") < 0.01 &
       (variant.int.attribute(vc = vc, att = "N_HOMALT") == 0 | variant.id(vc)=="gnomAD-SV_v2.1_DEL_15_146284") &  ## the OCA2 2.7kb exon 7 deletion is the only one with hmz in gnomad but pathogenic
       variant.qual(vc) > 500 & 
       logic.segdup(vc))
    {
      bcf.write.variant(out,vc)
    }
    
  }
}
# dispose the vcf reader
bcf.close(fp)
# show
bcf.close(out)

## VEP annotation is done with the VEP web portal. Below is the equivilant command
## ./vep --appris --biotype --buffer_size 500 --check_existing --distance 5000 --mane --merged --plugin LoFtool,[path_to]/LoFtool_scores.txt --polyphen b --pubmed --regulatory --sift b --species homo_sapiens --symbol --transcript_version --tsl --cache --input_file [input_data] --output_file [output_file] --port 3337
filename <- "./input.bigfiles/gnomad/gnomad_v2.1_sv.sites.Lof.nosegdup.vep.vcf"
fp <- bcf.open(filename,FALSE)
# error on opening (exit 0 for tests)
if(is.null(fp)) quit(save="no",status=0,runLast=FALSE)
# current variant
vc <- NULL
df.allele.hgnc <- NULL
while(!is.null(vc<-bcf.next(fp))) 
{
  ## affect all Refseq transcripts 
  df.vepparse <- variant.vep(vc) %>% filter(startsWith(Feature, "NM_")) %>% mutate(HGNC_ID = paste0("HGNC:", HGNC_ID))
  df.vepparse.alltx <- df.vepparse %>% group_by(HGNC_ID) %>%
    summarise(perc_tx_affected = length(grep("coding_sequence", Consequence))/n()) %>%
    filter(perc_tx_affected == 1) %>% 
    ungroup
  
  ## recurrent large CNVs are present in gnomadSV. In order to not double count them, remove recurrent CNVs from gnomadSV. 
  gr.vc <- GRanges(variant.chrom(vc), IRanges(variant.start(vc), variant.end(vc)))
  hits <- findOverlaps(gr.vc, gr.recurrentdel.hg19)
  if(length(hits)==0) {is.NAHR <- FALSE} else
  {
    is.NAHR <- width(subsetByOverlaps(gr.vc, gr.recurrentdel.hg19))/max(width(gr.vc), width(gr.recurrentdel.hg19[subjectHits(hits)])) > 0.4
  }
  
  if(nrow(df.vepparse.alltx)>0 & 
     !is.NAHR)
  {
    df.allele.hgnc <- rbind(df.allele.hgnc, 
                           data.frame(CSRA = variant.id(vc), 
                                      Chrom = variant.chrom(vc),
                                      Start = variant.start(vc), 
                                      End = variant.end(vc),
                                      HGNC_ID = unique(df.vepparse.alltx$HGNC_ID), 
                                      AF = variant.float.attribute(vc = vc, att = "AF"), 
                                      AFR_AF = variant.float.attribute(vc = vc, att = "AFR_AF"), 
                                      AMR_AF = variant.float.attribute(vc = vc, att = "AMR_AF"), 
                                      EAS_AF = variant.float.attribute(vc = vc, att = "EAS_AF"), 
                                      EUR_AF = variant.float.attribute(vc = vc, att = "EUR_AF"),
                                      MALE_AC = variant.int.attribute(vc = vc, att = "MALE_AC"),
                                      FEMALE_AC = variant.int.attribute(vc = vc, att = "FEMALE_AC"),
                                      AFR_AC = variant.int.attribute(vc = vc, att = "AFR_AC"),
                                      AFR_MALE_AC = variant.int.attribute(vc = vc, att = "AFR_MALE_AC"),
                                      AFR_FEMALE_AC = variant.int.attribute(vc = vc, att = "AFR_FEMALE_AC"),
                                      AMR_AC = variant.int.attribute(vc = vc, att = "AMR_AC"),
                                      AMR_MALE_AC = variant.int.attribute(vc = vc, att = "AMR_MALE_AC"),
                                      AMR_FEMALE_AC = variant.int.attribute(vc = vc, att = "AMR_FEMALE_AC"),
                                      EAS_AC = variant.int.attribute(vc = vc, att = "EAS_AC"),
                                      EAS_MALE_AC = variant.int.attribute(vc = vc, att = "EAS_MALE_AC"),
                                      EAS_FEMALE_AC = variant.int.attribute(vc = vc, att = "EAS_FEMALE_AC"),
                                      EUR_AC = variant.int.attribute(vc = vc, att = "EUR_AC"),
                                      EUR_MALE_AC = variant.int.attribute(vc = vc, att = "EUR_MALE_AC"),
                                      EUR_FEMALE_AC = variant.int.attribute(vc = vc, att = "EUR_FEMALE_AC"),
                                      OTH_AC = variant.int.attribute(vc = vc, att = "OTH_AC"),
                                      OTH_MALE_AC = variant.int.attribute(vc = vc, att = "OTH_MALE_AC"),
                                      OTH_FEMALE_AC = variant.int.attribute(vc = vc, att = "OTH_FEMALE_AC"),
                                      stringsAsFactors = F))
  }
}


## gnomad has certain SVs that are likely representing the same SV. 
## These SVs are characterized by 
## (a) they overlap with each other, and 
## (b) they have different population frequency breakdowns compared with each other
## the following function combined with the following loop identifies these duplicate SVs, removes the extra duplicates, and retains the one with the highest MAF. 
remove.gnomadSV.overlapping.duplicates <- function(df)
{
  df <- df %>% mutate(key = paste(CSRA, HGNC_ID, sep="__"))
  gr.all.genes <- with(df, GRanges(Chrom, IRanges(Start, End), 
                       genes = HGNC_ID, 
                       AF = AF, 
                       AFR_AF = AFR_AF, 
                       AMR_AF = AMR_AF, 
                       EAS_AF = EAS_AF, 
                       EUR_AF = EUR_AF,
                       MALE_AC = MALE_AC,
                       FEMALE_AC = FEMALE_AC,
                       AFR_AC = AFR_AC,
                       AFR_MALE_AC = AFR_MALE_AC,
                       AFR_FEMALE_AC = AFR_FEMALE_AC,
                       AMR_AC = AMR_AC,
                       AMR_MALE_AC = AMR_MALE_AC,
                       AMR_FEMALE_AC = AMR_FEMALE_AC,
                       EAS_AC = EAS_AC,
                       EAS_MALE_AC = EAS_MALE_AC,
                       EAS_FEMALE_AC = EAS_FEMALE_AC,
                       EUR_AC = EUR_AC,
                       EUR_MALE_AC = EUR_MALE_AC,
                       EUR_FEMALE_AC = EUR_FEMALE_AC,
                       OTH_AC = OTH_AC,
                       OTH_MALE_AC = OTH_MALE_AC,
                       OTH_FEMALE_AC = OTH_FEMALE_AC,
                       key = key))

  gr.processed <- do.call(c, lapply(unique(gr.all.genes$genes), function(y) 
  {
    gr.per.gene <- gr.all.genes[gr.all.genes$genes==y]
    overlap <- findOverlaps(gr.per.gene)
    
    do.call(c,unlist(lapply(unique(queryHits(overlap)), function(x)
    {
      pair.idx <- setdiff(subjectHits(overlap[which(queryHits(overlap)==x)]), x)
      cutoff.percentoverlap <- 0.2
      cutoff.percentACmatching <- 17/17
      remove.x.or.not= if(length(pair.idx)==0) { FALSE} else
      {
        remove.x.logic.vector <- do.call(c, lapply(pair.idx, function(z)
        {
          overlapregion <- intersect(gr.per.gene[z], gr.per.gene[x])
          percentOverlap <- width(overlapregion)/max(width(gr.per.gene[z]), width(gr.per.gene[x]))
          z.ACs <- elementMetadata(gr.per.gene[z])[c("MALE_AC", "FEMALE_AC", "AFR_AC", "AFR_MALE_AC", "AFR_FEMALE_AC", 
                                                               "AMR_AC", "AMR_MALE_AC", "AMR_FEMALE_AC", "EAS_AC", "EAS_MALE_AC", 
                                                               "EAS_FEMALE_AC", "EUR_AC", "EUR_MALE_AC", "EUR_FEMALE_AC" , "OTH_AC", 
                                                               "OTH_MALE_AC", "OTH_FEMALE_AC")] 
          x.ACs <- elementMetadata(gr.per.gene[x])[c("MALE_AC", "FEMALE_AC", "AFR_AC", "AFR_MALE_AC", "AFR_FEMALE_AC", 
                                                               "AMR_AC", "AMR_MALE_AC", "AMR_FEMALE_AC", "EAS_AC", "EAS_MALE_AC", 
                                                               "EAS_FEMALE_AC", "EUR_AC", "EUR_MALE_AC", "EUR_FEMALE_AC" , "OTH_AC", 
                                                               "OTH_MALE_AC", "OTH_FEMALE_AC")]
          percentACmatching <- length(unlist(which(z.ACs==x.ACs)))/17
          if(percentACmatching == cutoff.percentACmatching & percentOverlap > cutoff.percentoverlap & (gr.per.gene$AF[x] < gr.per.gene$AF[z]))
          {TRUE} else if (percentACmatching == cutoff.percentACmatching & percentOverlap > cutoff.percentoverlap & (gr.per.gene$AF[x] == gr.per.gene$AF[z]) & (x<z))
          {TRUE} else {FALSE}
        })) 
      remove.x.logic.vector %>% any
      }
      if(!remove.x.or.not) {gr.per.gene[x]} else {NULL}
    })))
  }))
  df[df$key %in% gr.processed$key, ] %>% unique
}

i=1
input <- df.allele.hgnc
df.allele.hgnc.processed <- NULL
repeat{
  print(paste0("iteration:", i))
  input.length <- nrow(input)
  print(paste0("input length: ", input.length))
  df.allele.hgnc.processed <- remove.gnomadSV.overlapping.duplicates(input)
  output.length <- nrow(df.allele.hgnc.processed)
  if(input.length==output.length) { break } 
  input <- df.allele.hgnc.processed
  i=i+1
  print(paste0("output length: ", output.length))
}

gnomadSV.genome.table <- df.allele.hgnc.processed %>% dplyr::select(CSRA, hgnc_id=HGNC_ID, AF, AFR_AF, AMR_AF, EAS_AF, EUR_AF) %>% filter(hgnc_id !="HGNC:")
hgnc_table <- read_tsv("./input/hgnc_table.tsv")
gnomadSV.genome.table <- gnomadSV.genome.table %>% 
  left_join(hgnc_table) %>% 
  arrange(CSRA) %>% 
  dplyr::rename(., GeneSymbol = hgnc_symbol)

## all validation CNV obtained from HGMD
# validationCnv = c("gnomAD-SV_v2.1_DEL_16_152725", "gnomAD-SV_v2.1_DEL_22_181424", "gnomAD-SV_v2.1_DEL_2_21008", "gnomAD-SV_v2.1_DEL_2_21010", "gnomAD-SV_v2.1_DEL_15_146284", "gnomAD-SV_v2.1_DEL_15_146513")
# sapply(validationCnv, function(x){x %in% gnomadSV.genome.table$CSRA})

write_tsv(gnomadSV.genome.table, "./output/gnomadSV.genome.table.tsv")

