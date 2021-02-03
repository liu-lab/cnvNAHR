library(dplyr)
library(readr)
library(stringr)
library(GenomicRanges)
library(biomaRt)
library(VariantAnnotation)
library(rtracklayer)

## download gnomad SV from website
gnomadSV <- readVcf("./input.bigfiles/gnomad_v2.1_sv.sites.vcf.gz")
# system('tabix -p vcf ./input.bigfiles/gnomad_v2.1_sv.sites.vcf.gz')

Recurrentdel <- read_tsv("./output/NAHR_GRCh38.withgenes.tsv")  %>% mutate(seqnames=gsub("chr", "",seqnames))
## genes in the 17_GL000258v2_alt haplotype do not map properly. Therefore the region needs to be converted to the chromosome 17 equivilant. 
idx.17q23 <- which(Recurrentdel$seqnames=="17_GL000258v2_alt" & Recurrentdel$end == 1355937)
Recurrentdel$start[idx.17q23] <- 45607798	
Recurrentdel$end[idx.17q23] <- 46210001
Recurrentdel$seqnames[idx.17q23] <-"17"
Recurrentdel.hg19 <- read_tsv("./output/NAHR_GRCh37.withgenes.tsv") %>% mutate(seqnames=gsub("chr", "",seqnames))
idx.17q23.hg19 <- which(Recurrentdel.hg19$seqnames=="17_ctg5_hap1" & Recurrentdel.hg19$end == 1307210)
Recurrentdel.hg19$start[idx.17q23.hg19] <- 43685164	
Recurrentdel.hg19$end[idx.17q23.hg19] <- 44287367
Recurrentdel.hg19$seqnames[idx.17q23.hg19] <-"17"

Recurrentdel.w.freq <- Recurrentdel %>% filter(freq>0)
Recurrentdel.w.freqCMA <- Recurrentdel %>% filter(freqCMA>0)
gr.Recurrentdel.w.freq <-  with(Recurrentdel.w.freq, GRanges(seqnames=seqnames, IRanges(start = start, end= end )))
gr.Recurrentdel.w.freqCMA <-  with(Recurrentdel.w.freqCMA, GRanges(seqnames=seqnames, IRanges(start = start, end= end )))
## Mb size for total unique intervals with considerable population frequency
sum(width(reduce(gr.Recurrentdel.w.freq)))/1e6
sum(width(reduce(gr.Recurrentdel.w.freqCMA)))/1e6
## Number of genes involved
do.call(c,lapply(Recurrentdel.w.freq$allgenes, function(x){strsplit(x, ",")[[1]]})) %>% unique %>% length()
do.call(c,lapply(Recurrentdel.w.freqCMA$allgenes, function(x){strsplit(x, ",")[[1]]})) %>% unique %>% length()
## aggregate population allele frequency
sum(Recurrentdel.w.freq$freq)/1e6

gr.recurrentdel<- with(Recurrentdel, GRanges(seqnames = seqnames, IRanges(start = start, end = end), regionname = name))
gr.recurrentdel.hg19 <- with(Recurrentdel.hg19, GRanges(seqnames = seqnames, IRanges(start=start, end=end)))

Biallelic.gene.omim.ddd.clingen <- read_lines("./output/biallelicgenes.txt")

## gnomad sv is in GRCh37. So ranges for genes need to be in GRCh37. 
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
gene.raw.know.AR.grch37 <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position","external_gene_name",
                                     "hgnc_symbol", "hgnc_id"), 
                        filter="chromosome_name", values=c(1:22,"X", "Y"), mart=grch37) %>% 
                        mutate(hgnc_id=paste0("HGNC:", hgnc_id)) %>%
                        filter(hgnc_id %in% Biallelic.gene.omim.ddd.clingen)


## make filter indexes for LOF SVs. 
## 1. variant has PROTEIN_CODING__LOF annotation from the original gnomad file
## 2. in AR genes
## 3. no hmz counts in gnomad
## 4. popmax AF < 1%
## 5. pass filter; quality socre > 500
## 6. is deletion CNV
## 7. add exceptions (ABCC6 and OCA2 are identified bny reviewing all deletion pathogenic CNVs in AR genes in NAHR region from HGMD)
idx.sv.LOF0 <- c(which(sapply(info(gnomadSV)[["PROTEIN_CODING__LOF"]], function(x)
                                                                          {
                                                                            any(unlist(x) %in% gene.raw.know.AR.grch37$hgnc_symbol)!= 0
                                                                          }
                               )
                        ),
                  which(names(gnomadSV)=="gnomAD-SV_v2.1_DEL_16_152725")) %>% ## this SV is incorrectly annotated by gnomAD. The annotation missed the ABCC6 gene overlap so did not label it as LOF. 
              unique 

idx.sv.LOF <- intersect(idx.sv.LOF0, which(
    (info(gnomadSV)[["N_HOMALT"]] == 0 | names(gnomadSV)=="gnomAD-SV_v2.1_DEL_15_146284") &  ## the OCA2 2.7kb exon 7 deletion is the only one with hmz in gnomad but pathogenic
    info(gnomadSV)[["POPMAX_AF"]] < 0.01 &
    elementMetadata(rowRanges(gnomadSV))[,"FILTER"]=="PASS" &
    elementMetadata(rowRanges(gnomadSV))[,"QUAL"]>500 &
    paste(elementMetadata(rowRanges(gnomadSV))[,"ALT"],  collapse = ";")=="<DEL>"))

gnomadSV.LOF <- gnomadSV[idx.sv.LOF]
gr.gnomadSV.LOF <- rowRanges(gnomadSV.LOF)
end(gr.gnomadSV.LOF) <- info(gnomadSV.LOF)[["END"]]  ## modify end coordinates
gr.gnomadSV.LOF$AF <- unlist(info(gnomadSV.LOF)[["AF"]])
gr.gnomadSV.LOF$AFR_AF <- unlist(info(gnomadSV.LOF)[["AFR_AF"]])
gr.gnomadSV.LOF$AMR_AF <- unlist(info(gnomadSV.LOF)[["AMR_AF"]])
gr.gnomadSV.LOF$EAS_AF <- unlist(info(gnomadSV.LOF)[["EAS_AF"]])
gr.gnomadSV.LOF$EUR_AF <- unlist(info(gnomadSV.LOF)[["EUR_AF"]])

## additional filtering steps to remove low quality CNVs that are less likely to be pathogenic. 
## rule: CNV cannot be within segdup
mySession = browserSession("UCSC")
genome(mySession) <- "hg19"
segdup19.raw <- getTable(ucscTableQuery(mySession, track="genomicSuperDups")) %>% mutate(chrom = gsub("chr", "", chrom))
gr.segdup19 <- with(segdup19.raw, GRanges( seqnames = chrom,  ranges = IRanges(chromStart, end = chromEnd))) %>% reduce
hits <- findOverlaps(gr.gnomadSV.LOF, gr.segdup19)

## remove gnomad cnvs that has >80% content located within any segdup region.
logical.overlap.w.segdup <- do.call(c,lapply(unique(queryHits(hits)), function(x)
{
  hits.select <- hits[queryHits(hits)==x]
  total.overlap.width.insegdup.select <- pintersect(gr.gnomadSV.LOF[queryHits(hits.select)], gr.segdup19[subjectHits(hits.select)]) %>% reduce %>% width %>% sum
  if(total.overlap.width.insegdup.select/width(gr.gnomadSV.LOF[unique(queryHits(hits.select))]) > 0.8)
  {
    x
  }
}))

if(!is.null(logical.overlap.w.segdup))
{
  gr.gnomadSV.LOF.no.segdup <- gr.gnomadSV.LOF[-logical.overlap.w.segdup]
} else 
{
  gr.gnomadSV.LOF.no.segdup <- gr.gnomadSV.LOF
}

## add gene annotatoin to gr.gnomadSV.LOF.no.segdup
gr.gene.recessive.grch37 <- with(gene.raw.know.AR.grch37, GRanges( seqnames = chromosome_name,  ranges = IRanges(start_position, end = end_position),
                                                                   genename=hgnc_symbol, name=hgnc_id ))
overlaps <- findOverlaps(gr.gnomadSV.LOF.no.segdup, gr.gene.recessive.grch37)
overlaps.genes <- as.data.frame(overlaps) %>%
  mutate(genename = gr.gene.recessive.grch37$genename[subjectHits(overlaps)]) %>%
  group_by(queryHits) %>% mutate(genenames= paste0(unique(genename), collapse = ",")) %>%
  dplyr::select(queryHits, genenames) %>% unique %>% arrange(queryHits)
gr.gene.recessive.grch37.gnomadSV.LOF <- gr.gnomadSV.LOF.no.segdup
gr.gene.recessive.grch37.gnomadSV.LOF$genes <- overlaps.genes$genenames

## recurrent large CNVs are present in gnomadSV. In order to not double count them, remove recurrent CNVs from gnomadSV. 
hits <- findOverlaps(gr.gene.recessive.grch37.gnomadSV.LOF, gr.recurrentdel.hg19)
overlaps <- pintersect(gr.gene.recessive.grch37.gnomadSV.LOF[queryHits(hits)], gr.recurrentdel.hg19[subjectHits(hits)])
percentOverlap1 <- width(overlaps) / width(gr.gene.recessive.grch37.gnomadSV.LOF[queryHits(hits)])
percentOverlap2 <- width(overlaps) / width(gr.recurrentdel.hg19[subjectHits(hits)])
hits <- hits[percentOverlap1 > 0.4 & percentOverlap2 > 0.4] # NPHP1 overlaps 46.3%
## these are the recurrent CNVs in gnomadSV
gr.gene.recessive.grch37.gnomadSV.LOF[unique(queryHits(hits))]
## ruling out recurrent CNVs
gr.gene.recessive.grch37.gnomadSV.LOF1 <- gr.gene.recessive.grch37.gnomadSV.LOF[-unique(queryHits(hits))]

## next, need to further filter these SVs based on their consequence to the genes. 
## first, annotation is down with VEP
bed.gene.recessive.grch37.gnomadSV.LOF1 <- data.frame(seqnames=seqnames(gr.gene.recessive.grch37.gnomadSV.LOF1),
                                                      starts=start(gr.gene.recessive.grch37.gnomadSV.LOF1)+1,
                                                      ends=end(gr.gene.recessive.grch37.gnomadSV.LOF1),
                                                      names=names(gr.gene.recessive.grch37.gnomadSV.LOF1))

write.table(bed.gene.recessive.grch37.gnomadSV.LOF1, file="./output/gnomad_v2.1_sv.sites.AR.LOF.nosegdup.bed", quote=F, sep="\t", row.names=F, col.names=F)

## VEP annotation is done with the bed file on the VEP web portal. Check RefSeq transcripts. Below is the equivilant command
## ./vep --af --appris --biotype --buffer_size 500 --check_existing --distance 5000 --mane --polyphen b --pubmed --refseq --regulatory --sift b --species homo_sapiens --symbol --transcript_version --tsl --cache --input_file [input_data] --output_file [output_file] --port 3337
## load the annotated file
gnomad_v2.1_sv.sites.AR.LOF.nosegdup.vcf.vep <- read_tsv('./input.bigfiles/gnomad_v2.1_sv.sites.AR.LOF.nosegdup.vcf.vep.refseq.txt')
## filter down to only recessive genes and refseq transcript
gnomad_v2.1_sv.sites.AR.LOF.nosegdup.vcf.vep.1 <- gnomad_v2.1_sv.sites.AR.LOF.nosegdup.vcf.vep[gnomad_v2.1_sv.sites.AR.LOF.nosegdup.vcf.vep$SYMBOL %in% gene.raw.know.AR.grch37$hgnc_symbol & 
                                                                                                 startsWith(gnomad_v2.1_sv.sites.AR.LOF.nosegdup.vcf.vep$Feature, "NM_"),]
## detailed filtering logic
## 1. SV has to affect all refseq transcripts altering coding sequence. 
all.tx.coding <- gnomad_v2.1_sv.sites.AR.LOF.nosegdup.vcf.vep.1 %>% 
  group_by(Allele, SYMBOL) %>% 
  summarise(perc_tx_affected = length(grep("coding_sequence", Consequence))/n()) %>%
  filter(perc_tx_affected == 1) %>% ungroup %>% mutate(key= paste0(Allele, "_", SYMBOL))

gr.gene.recessive.grch37.gnomadSV.LOF2 <- gr.gene.recessive.grch37.gnomadSV.LOF1[names(gr.gene.recessive.grch37.gnomadSV.LOF1) %in% 
                                                                                   all.tx.coding$Allele]

## 2. gnomad has certain SVs that are likely representing the same SV. 
## These SVs are characterized by 
## (a) they overlap with each other, and 
## (b) they have different population frequency breakdowns compared with each other
## the following function combined with the following loop identifies these duplicate SVs, removes the extra duplicates, and retains the one with the highest MAF. 
remove.gnomadSV.overlapping.duplicates <- function(gr.all.genes)
{
  do.call(c, lapply(unique(gr.all.genes$genes), function(y) 
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
        
        #print(paste0("   index=", x, "; name is ", names(gr.per.gene[x])))
        remove.x.logic.vector <- do.call(c, lapply(pair.idx, function(z)
        {
          overlapregion <- intersect(gr.per.gene[z], gr.per.gene[x])
          percentOverlap <- width(overlapregion)/max(width(gr.per.gene[z]), width(gr.per.gene[x]))
          z.ACs <- info(gnomadSV.LOF[names(gr.per.gene[z])])[c("MALE_AC", "FEMALE_AC", "AFR_AC", "AFR_MALE_AC", "AFR_FEMALE_AC", 
                                                        "AMR_AC", "AMR_MALE_AC", "AMR_FEMALE_AC", "EAS_AC", "EAS_MALE_AC", 
                                                        "EAS_FEMALE_AC", "EUR_AC", "EUR_MALE_AC", "EUR_FEMALE_AC" , "OTH_AC", 
                                                        "OTH_MALE_AC", "OTH_FEMALE_AC")] 
          x.ACs <- info(gnomadSV.LOF[names(gr.per.gene[x])])[c("MALE_AC", "FEMALE_AC", "AFR_AC", "AFR_MALE_AC", "AFR_FEMALE_AC", 
                                                      "AMR_AC", "AMR_MALE_AC", "AMR_FEMALE_AC", "EAS_AC", "EAS_MALE_AC", 
                                                      "EAS_FEMALE_AC", "EUR_AC", "EUR_MALE_AC", "EUR_FEMALE_AC" , "OTH_AC", 
                                                      "OTH_MALE_AC", "OTH_FEMALE_AC")]
          percentACmatching <- length(unlist(which(z.ACs==x.ACs)))/17
          if(percentACmatching == cutoff.percentACmatching & percentOverlap > cutoff.percentoverlap & (gr.per.gene$AF[x] < gr.per.gene$AF[z]))
          {TRUE} else if (percentACmatching == cutoff.percentACmatching & percentOverlap > cutoff.percentoverlap & (gr.per.gene$AF[x] == gr.per.gene$AF[z]) & (x<z))
          {TRUE} else {FALSE}
        })) 
        #print(paste0("   for subject names: ", paste(names(gr.per.gene[pair.idx]), collapse = ",") , "; logic vector is, ", 
        #             paste(remove.x.logic.vector, collapse = ",")))
        remove.x.logic.vector %>% any
        
      }
      #print(paste0(remove.x.or.not, " is the decision for index ", x, " removal. "))
      if(!remove.x.or.not) {gr.per.gene[x]} else {NULL}
    })))
  }))
}

i=1
gr.gene.recessive.grch37.gnomadSV.LOF3 <- gr.gene.recessive.grch37.gnomadSV.LOF2

repeat{
  print(paste0("iteration:", i))
  input.length <- length(gr.gene.recessive.grch37.gnomadSV.LOF3)
  gr.gene.recessive.grch37.gnomadSV.LOF4 <- remove.gnomadSV.overlapping.duplicates(gr.gene.recessive.grch37.gnomadSV.LOF3)
  output.length <- length(gr.gene.recessive.grch37.gnomadSV.LOF4)
  if(input.length!=output.length)
  {
    gr.gene.recessive.grch37.gnomadSV.LOF3 <<- gr.gene.recessive.grch37.gnomadSV.LOF4
  } else {
    gr.gene.recessive.grch37.gnomadSV.LOF3 <<- gr.gene.recessive.grch37.gnomadSV.LOF4
    break}
  i=i+1
  print(paste0("input length: ", input.length))
  print(paste0("output length: ", output.length))
}

## cleaned up gnomadSV with LOF in recessive genes
gnomadSV.af.by.region.exclude.NAHR <- as_tibble(cbind(elementMetadata(gr.gene.recessive.grch37.gnomadSV.LOF3)[,6:11], 
                                                      data.frame(CSRA=names(gr.gene.recessive.grch37.gnomadSV.LOF3), 
                                                                 stringsAsFactors = F) ))
names(gnomadSV.af.by.region.exclude.NAHR)[6] <- "GeneSymbol"

# validationCnv = c("gnomAD-SV_v2.1_DEL_16_152725", "gnomAD-SV_v2.1_DEL_22_181424", "gnomAD-SV_v2.1_DEL_2_21008", "gnomAD-SV_v2.1_DEL_2_21010", "gnomAD-SV_v2.1_DEL_15_146284", "gnomAD-SV_v2.1_DEL_15_149078", "gnomAD-SV_v2.1_DEL_15_146513")
# sapply(validationCnv, function(x){x %in% gnomadSV.af.by.region.exclude.NAHR$CSRA})
# validationCnv[1] %in% names(gr.gene.recessive.grch37.gnomadSV.LOF2)
# validationCnv[1] %in% names(gr.gene.recessive.grch37.gnomadSV.LOF3)
# validationCnv[1] %in% names.to.include
# gr.gene.recessive.grch37.gnomadSV.LOF2[names(gr.gene.recessive.grch37.gnomadSV.LOF2) == validationCnv[3]]

## convert variant based summary to gene based summary
gnomadSV.af.by.region.exclude.NAHR.multigene <- gnomadSV.af.by.region.exclude.NAHR[grepl(",",gnomadSV.af.by.region.exclude.NAHR$GeneSymbol),]
gnomadSV.af.by.region.exclude.NAHR.multigene.sep <- bind_rows(lapply(1:nrow(gnomadSV.af.by.region.exclude.NAHR.multigene), function(x)
  {
   genes = strsplit(gnomadSV.af.by.region.exclude.NAHR.multigene$GeneSymbol[x], ",")[[1]]
   tibble(AF=rep(gnomadSV.af.by.region.exclude.NAHR.multigene$AF[x], length(genes)), 
          AFR_AF=rep(gnomadSV.af.by.region.exclude.NAHR.multigene$AFR_AF[x], length(genes)), 
          AMR_AF=rep(gnomadSV.af.by.region.exclude.NAHR.multigene$AMR_AF[x], length(genes)), 
          EAS_AF=rep(gnomadSV.af.by.region.exclude.NAHR.multigene$EAS_AF[x], length(genes)), 
          EUR_AF=rep(gnomadSV.af.by.region.exclude.NAHR.multigene$EUR_AF[x], length(genes)), 
          GeneSymbol= genes, 
          CSRA= rep(gnomadSV.af.by.region.exclude.NAHR.multigene$CSRA[x], length(genes))
         )
  })) %>% mutate(key= paste0(CSRA, "_", GeneSymbol)) %>%
  filter(key %in% all.tx.coding$key) %>%   ## remove instances in which the gene is included but not LOF
  dplyr::select(-key)
gnomadSV.af.by.region.exclude.NAHR.singlegene <- gnomadSV.af.by.region.exclude.NAHR[!grepl(",",gnomadSV.af.by.region.exclude.NAHR$GeneSymbol),]
gnomadSV.af.by.region.exclude.NAHR.clean <- bind_rows(gnomadSV.af.by.region.exclude.NAHR.singlegene, gnomadSV.af.by.region.exclude.NAHR.multigene.sep)

## add NAHR large CNV back
gnomadSV.af.by.region <- rbind(gnomadSV.af.by.region.exclude.NAHR.clean, 
                               do.call(rbind, lapply(1:nrow(Recurrentdel.w.freq), function(x){
                                 genenames.inCNV <- strsplit(Recurrentdel.w.freq$biallelic[x], ",")[[1]]
                                 if(!is.na(genenames.inCNV))
                                 {
                                   data.frame(AF=Recurrentdel.w.freq$freq[x]/1e6, 
                                              AFR_AF=Recurrentdel.w.freq$freq[x]/1e6, 
                                              AMR_AF=Recurrentdel.w.freq$freq[x]/1e6, 
                                              EAS_AF=Recurrentdel.w.freq$freq[x]/1e6, 
                                              EUR_AF=Recurrentdel.w.freq$freq[x]/1e6, 
                                              GeneSymbol=genenames.inCNV, 
                                              CSRA= Recurrentdel.w.freq$name[x])
                                 }
                               })) 
                              )    

# convert.byregion.to.bygene <- function(gnomadSV.af.by.region)
# {
#   do.call(rbind,lapply(1:nrow(gnomadSV.af.by.region), function(x){
#     genenames <- strsplit(gnomadSV.af.by.region$GeneSymbol[x], ",")[[1]]
#     data.frame(AF=rep(gnomadSV.af.by.region$AF[x], length(genenames)), 
#                genes=genenames)
#   })) %>% arrange(AF) %>% group_by(genes) %>% 
#     summarise(AFsum=sum(AF)) %>% ungroup() %>% mutate(genes=as.character(genes)) %>% arrange(desc(AFsum))
# }
# 
# gnomadSV.af.by.gene <- convert.byregion.to.bygene(gnomadSV.af.by.region)
# 
# gnomadSV.af.by.gene.exclude.NAHR <- convert.byregion.to.bygene(gnomadSV.af.by.region.exclude.NAHR)
# 
# gnomadSV.af.by.gene.exclude.NAHR.rename <- gnomadSV.af.by.gene.exclude.NAHR %>% mutate(AFsum.exc.NAHR=AFsum) %>% dplyr::select(genes, AFsum.exc.NAHR)


## clinvar, everything needs to be done in GRCh38

ensembl <- useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", 
                      dataset="hsapiens_gene_ensembl",
                      mirror = "uswest")

gene.know.AR.bm <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position","external_gene_name",
                                     "hgnc_symbol", "hgnc_id"), 
                        filter="hgnc_id", values=Biallelic.gene.omim.ddd.clingen, mart=ensembl) %>% filter(chromosome_name %in% c(1:22, "X", "Y"))

gr.genes.biallelic <- with(gene.know.AR.bm, GRanges(seqnames = chromosome_name, IRanges(start= start_position, end = end_position)))

## there are three input files for clinvar
## 1. vcf is used to obtain AF
## 2. submission_summary is used to obtain curation text to facilitate filtering
## 3. variant_summary is used as the main file 

clinvar.raw.ar=readVcf("./input.bigfiles/clinvar_20210119.vcf.gz", genome = "GRCh38", param = reduce(gr.genes.biallelic))
writeVcf(clinvar.raw.ar, "./input.bigfiles/clinvar_20210119_AR.vcf")
# In order to get AF, annotate with VEP web portal. 
# ./vep --af_gnomad --appris --biotype --buffer_size 500 --check_existing --distance 5000 --mane --polyphen b --pubmed --refseq --regulatory --sift b --species homo_sapiens --symbol --transcript_version --tsl --cache --input_file [input_data] --output_file [output_file]
# check RefSeq transcripts. below is the result loaded
clinvar.rec <- readVcf("./input.bigfiles/clinvar_20210119_AR.annotated.vcf.gz", genome = "GRCh38")
## select pathogenic/likely pathogenic variants
idx.plp <- which(paste(info(clinvar.rec)[["CLNSIG"]], collapse = ";") %>% 
                   unlist %in% c("Pathogenic", "Pathogenic/Likely_pathogenic", "Likely_pathogenic", 
                                 "Pathogenic/Likely pathogenic, drug response", "Pathogenic, drug response", "Likely pathogenic, drug response") &
                   paste(info(clinvar.rec)[["CLNREVSTAT"]], collapse = "_") %>% unlist %in%
                   c("criteria_provided__single_submitter", "criteria_provided__multiple_submitters__no_conflicts",
                     "criteria_provided__conflicting_interpretations",  "reviewed_by_expert_panel", "practice_guideline"  ) 
                 # this selection is equivilant to one star or more, i.e. excluding "no_assertion_criteria_provided", "no_interpretation_for_the_single_variant", "no_assertion_provided"
)
clinvar.rec.plp <- clinvar.rec[idx.plp]
gr.clinvar.rec.plp <- rowRanges(clinvar.rec.plp)
## 32, 33, 34, 36, 38 are the number of column for AF, AFR_AF, AMR_AF, EAS_AF, NFE_AF (I changed this to EUR_AF) within the VEP annotation in this particular annotation file
afs.gr.clinvar.rec.plp <- bind_rows(lapply(1:length(clinvar.rec.plp), function(x)
  {
    if(x%%1000==0) {print(round(x/53304*100))}
    as.data.frame(t(as.numeric(strsplit(unlist(info(clinvar.rec.plp[x])[["CSQ"]]) , "\\|")[[1]][c(32, 33, 34, 36, 38)])))
  }))
afs.gr.clinvar.rec.plp[is.na(afs.gr.clinvar.rec.plp)] <- 0
gr.clinvar.rec.plp$AF <- afs.gr.clinvar.rec.plp$V1
gr.clinvar.rec.plp$AFR_AF <- afs.gr.clinvar.rec.plp$V2
gr.clinvar.rec.plp$AMR_AF <- afs.gr.clinvar.rec.plp$V3
gr.clinvar.rec.plp$EAS_AF <- afs.gr.clinvar.rec.plp$V4
gr.clinvar.rec.plp$EUR_AF <- afs.gr.clinvar.rec.plp$V5


csra.gnomad <- data.frame(CSRA= paste(seqnames(gr.clinvar.rec.plp), start(gr.clinvar.rec.plp), as.character(gr.clinvar.rec.plp$REF), 
                                      as.character(unlist(gr.clinvar.rec.plp$ALT)), sep="_"), 
                          AF=gr.clinvar.rec.plp$AF, 
                          AFR_AF=gr.clinvar.rec.plp$AFR_AF, 
                          AMR_AF=gr.clinvar.rec.plp$AMR_AF, 
                          EAS_AF=gr.clinvar.rec.plp$EAS_AF, 
                          EUR_AF=gr.clinvar.rec.plp$EUR_AF,
                          stringsAsFactors = F)

## obtain variant curation description from the submission summary file downloaded from the clinvar FTP
clinvar.subsum <- read_tsv("./input.bigfiles/submission_summary_2021-01.txt.gz", skip=15)
clinvar.subsum1 <- clinvar.subsum %>% mutate(VariationID=`#VariationID`) %>% 
  group_by(VariationID) %>% summarise(Description=paste(Description, collapse="||"))

## merge the descriptions to the variant summary file downloaded from the clinvar FTP website. 
clinvar.varsum <- read_tsv("./input.bigfiles/variant_summary_2021-01.txt.gz") %>% 
  filter(Assembly == "GRCh38" & ClinicalSignificance %in% 
           c("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic", 
             "Pathogenic/Likely pathogenic, drug response", "Pathogenic, drug response", "Likely pathogenic, drug response") &
           ReviewStatus %in% c("criteria provided, multiple submitters, no conflicts", "criteria provided, single submitter",
                               "practice guideline", "reviewed by expert panel")) %>%
  left_join(clinvar.subsum1)

## select only recessive genes
clinvar.varsum.recessive <- clinvar.varsum %>% filter(HGNC_ID %in% Biallelic.gene.omim.ddd.clingen) %>%
  mutate(CSRA=paste0(Chromosome, "_", PositionVCF, "_", ReferenceAlleleVCF, "_", AlternateAlleleVCF))

# remove some variants that are not recessive alleles 
idx.sel1 <- do.call(c,sapply(1:nrow(clinvar.varsum.recessive), function(x){
  if(grepl("dominant negative", clinvar.varsum.recessive$Description[x])){
    x
  }
}))

idx.sel2 <- do.call(c,sapply(1:nrow(clinvar.varsum.recessive), function(x){
  if(grepl("his variant .* an autosomal dominant", clinvar.varsum.recessive$Description[x])){
    x
  }
}))

clinvar.varsum.recessive <- clinvar.varsum.recessive[-c(idx.sel1, idx.sel2),]

## merge with AF annotation baesd on CSRA unique match
clinvar.varsum.recessive.gnomad <- left_join(clinvar.varsum.recessive,csra.gnomad)

## exclude variants with relatively high AF but low confidence curation
low.confidence.exclude <- with(clinvar.varsum.recessive.gnomad, 
                               ((AF>0.005 & NumberSubmitters<=3) | 
                                  (AF>0.001 &
                                     ReviewStatus %in% c("no assertion criteria provided", 
                                                         "criteria provided, single submitter"))))
clinvar.varsum.recessive.gnomad <- clinvar.varsum.recessive.gnomad[!low.confidence.exclude,]

clinvar.varsum.recessive.gnomad.reformat <- clinvar.varsum.recessive.gnomad %>% dplyr::select(GeneSymbol, CSRA, AF, AFR_AF, AMR_AF, EAS_AF, EUR_AF) %>% unique
clinvar.varsum.recessive.gnomad.reformat %>% nrow

## this is eveything combined showing frequency by allele
clinvar.gnomadsv.recessive <- rbind(clinvar.varsum.recessive.gnomad.reformat, 
                                              gnomadSV.af.by.region)



## compare ar gene allele freq with CNV allele freq
Recurrentdel.highAF <- Recurrentdel %>% filter(seqnames !="X" & freq>0) %>% arrange(desc(freq), desc(freqCMA))
NAHR.ARgenes <- do.call(c, lapply(1:nrow(Recurrentdel.highAF), function(x){
  strsplit(Recurrentdel.highAF$biallelic[x], ",")[[1]]})) %>% unique %>% setdiff(., NA)

## all alleles for recessive genes within NAHR regions
clinvar.gnomadsv.recessive.NAHR <- clinvar.gnomadsv.recessive %>%  filter(GeneSymbol %in% NAHR.ARgenes)

metrics.by.population <- function(AFname)
{
  clinvar.gnomadsv.recessive.pop <- clinvar.gnomadsv.recessive %>% dplyr::select(GeneSymbol, CSRA, AFname)
  names(clinvar.gnomadsv.recessive.pop)[3] <- "AFname1"
  clinvar.gnomadsv.recessive.by.gene <- clinvar.gnomadsv.recessive.pop %>% group_by(GeneSymbol) %>% summarise(AFsum=sum(AFname1)) %>% arrange(desc(AFsum)) %>% ungroup
  
  # all NAHR region enclosed genes- total carrier frequency (with NAHR CNV freq included)
  NAHR.ARgenes.freq <- clinvar.gnomadsv.recessive.by.gene %>% filter(GeneSymbol %in% NAHR.ARgenes)
  #  all NAHR region enclosed genes- carrier freq from NAHR CNVs only
  cnv.recessive.by.gene <- do.call(rbind, lapply(1:nrow(Recurrentdel.highAF), function(x){
    genenames.inCNV <- strsplit(Recurrentdel.highAF$biallelic[x], ",")[[1]]
    if(!is.na(genenames.inCNV))
    {
      data.frame(GeneSymbol=genenames.inCNV, AFnahr=Recurrentdel.highAF$freq[x]/1e6, region= Recurrentdel.highAF$name[x])
    }
  })) 
  
  pct.change.NAHR.ARgenes <- left_join(NAHR.ARgenes.freq, cnv.recessive.by.gene) %>% mutate(percentNAHR= AFnahr/AFsum)
  
  ## this is percent of ar genes significantly affected by NAHR for carrier frequencies. 
  percentGenesnahr10 <- length(pct.change.NAHR.ARgenes$percentNAHR[pct.change.NAHR.ARgenes$percentNAHR>0.1])/length(Biallelic.gene.omim.ddd.clingen)
  print(percentGenesnahr10)
  genesNAHR10ranked <- pct.change.NAHR.ARgenes[pct.change.NAHR.ARgenes$percentNAHR>0.1,] %>% arrange(desc(percentNAHR)) %>% as.data.frame()
  write_tsv(genesNAHR10ranked, path = paste0("./output/genesNAHR10ranked_", AFname, "_percGeneis_", round(percentGenesnahr10*100, 2), ".tsv"))
}

metrics.by.population("AF")
metrics.by.population("AFR_AF")
metrics.by.population("AMR_AF")
metrics.by.population("EAS_AF")
metrics.by.population("EUR_AF")

