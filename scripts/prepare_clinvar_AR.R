library(dplyr)
library(readr)
library(GenomicRanges)
library(biomaRt)
library(VariantAnnotation)
library(parallel)

## Extract P and LP variants from ClinVar. There are three input files for ClinVar
## 1. vcf is used to obtain AF
## 2. submission_summary is used to obtain curation text to facilitate filtering
## 3. variant_summary is used as the main file 
Biallelic.gene.omim.ddd.clingen <- read_lines("./output/biallelicgenes.txt")

ensembl <- useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", 
                      dataset="hsapiens_gene_ensembl")

gene.know.AR.bm <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position","external_gene_name",
                                      "hgnc_symbol", "hgnc_id"), 
                         filter="hgnc_id", values=Biallelic.gene.omim.ddd.clingen, mart=ensembl) %>% filter(chromosome_name %in% c(1:22, "X", "Y"))

gr.genes.biallelic <- with(gene.know.AR.bm, GRanges(seqnames = chromosome_name, IRanges(start= start_position, end = end_position), hgnc_symbol=hgnc_symbol))
clinvar.raw.ar=readVcf("./input.bigfiles/clinvar_20210119.vcf.gz", genome = "GRCh38", param = reduce(gr.genes.biallelic))
writeVcf(clinvar.raw.ar, "./input.bigfiles/clinvar_20210119_AR.vcf")
system("bgzip ./input.bigfiles/clinvar_20210119_AR.vcf")
# In order to get AF, annotate with VEP web portal. 
# ./vep --af_gnomad --appris --biotype --buffer_size 500 --check_existing --distance 5000 --mane --polyphen b --pubmed --refseq --regulatory --sift b --species homo_sapiens --symbol --transcript_version --tsl --cache --input_file [input_data] --output_file [output_file]
# check RefSeq transcripts. below is the result loaded
clinvar.rec <- readVcf("./input.bigfiles/clinvar_20210119_AR_annotated.vcf.gz", genome = "GRCh38")
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
## get the number of columns for AF, AFR_AF, AMR_AF, EAS_AF, NFE_AF (I changed this to EUR_AF) within the VEP annotation in this particular annotation file
gnomad.col.idx <- sapply(c("gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_EAS_AF", "gnomAD_NFE_AF"), function(x){which(strsplit(info(header(clinvar.rec.plp))["CSQ","Description"], "\\|")[[1]]==x)}) %>% unname

afs.gr.clinvar.rec.plp <- bind_rows(mclapply(1:length(clinvar.rec.plp), function(x)
{
  as.data.frame(t(as.numeric(strsplit(unlist(info(clinvar.rec.plp[x])[["CSQ"]]) , "\\|")[[1]][gnomad.col.idx])))
}, mc.cores = 16))

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

clinvar.varsum.recessive.gnomad.reformat <- clinvar.varsum.recessive.gnomad %>% dplyr::select(hgnc_id=HGNC_ID, CSRA, AF, AFR_AF, AMR_AF, EAS_AF, EUR_AF) %>% unique
## the two PRRT2 variants in the current annotation do not have accurate frequencies. Removing them and adding back the entires with edited frequencies
## excluding the RBM8A hypomorphic variant. This is not purely a recessive allele. Homozygotes of this allele do not cause disease. 
clinvar.varsum.recessive.gnomad.reformat <- clinvar.varsum.recessive.gnomad.reformat %>% filter(!CSRA %in% c("16_29813694_GC_G", "16_29813694_G_GC", "1_145927328_C_G"))
## change gene name for HBA2 into HBA1
clinvar.varsum.recessive.gnomad.reformat <- clinvar.varsum.recessive.gnomad.reformat %>% 
  mutate(hgnc_id = replace(hgnc_id, hgnc_id == "HGNC:4824", "HGNC:4823"))
## adding more variants curated to be P but not currently in gnomadSV. Variant frequency source is from PMID:27533158. SMN1, HBA1/HBA2, CYP21A2
## adding hypomorphic allels from RBM8A and TBX6
## adding a SLC25A1 variant that is curated to be pathogenic by us but is located in a low complexity region
clinvar.varsum.recessive.gnomad.reformat <- read_tsv("./input/clinvar.gnomadsv.recessive.supplement.tsv") %>% 
  bind_rows(clinvar.varsum.recessive.gnomad.reformat, .) %>% unique

hgnc_table <- read_tsv("./input/hgnc_table.tsv")
clinvar.varsum.recessive.gnomad.reformat <- clinvar.varsum.recessive.gnomad.reformat %>% 
  left_join(hgnc_table) %>% 
  arrange(CSRA) %>% 
  dplyr::rename(., GeneSymbol = hgnc_symbol)

write_tsv(clinvar.varsum.recessive.gnomad.reformat, "./output/clinvar.recessive.table.tsv")
