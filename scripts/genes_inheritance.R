library(tidyr)
library(stringr)
library(dplyr)
library(readr)
library(biomaRt)
library(GenomicRanges)
library(parallel)


ensembl <- useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", 
                      dataset="hsapiens_gene_ensembl")

attributes <- c("ensembl_gene_id","chromosome_name","start_position","end_position","hgnc_symbol", "hgnc_id", "ensembl_transcript_id", 
                "entrezgene_id", "gene_biotype") 

## try to use hgnc_id as unique keys for genes. Make a convert table between hgnc_id and hgnc_symbol
hgnc_table <- getBM(attributes=c("hgnc_symbol", "hgnc_id"), values="*", mart=ensembl) 
write_tsv(hgnc_table, "./input/hgnc_table.tsv")
coding_bm <- getBM(attributes=c("gene_biotype", "hgnc_id"), values="*", mart=ensembl)
coding_HGNC_ID <- coding_bm %>% filter(gene_biotype=="protein_coding") %>% .$hgnc_id %>% unique %>% setdiff("")
write_lines(coding_HGNC_ID, "./input/coding_HGNC_ID.txt")



# load input data
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
                    " Autosomal recessive, Autosomal dominant", 
                    " Autosomal recessive, Multifactorial, Autosomal dominant",   
                    " Autosomal recessive, Digenic dominant, Autosomal dominant", 
                    " Autosomal recessive, Autosomal dominant, Isolated cases" , 
                    " Autosomal recessive, Autosomal dominant, Somatic mutation", 
                    " Multifactorial, Autosomal recessive, Autosomal dominant", 
                    " Autosomal recessive, Autosomal dominant, Digenic dominant", 
                    " Autosomal recessive, Somatic mutation, Autosomal dominant")

## convert all gene names to hgnc gene id

Biallelic.OMIM.gene.value <- Genemap2_2 %>% filter(Inheritance %in% Biallelic.OMIM) %>% .$`Entrez Gene ID` %>% unique
Biallelic.OMIM.gene.bm <- getBM(attributes=attributes, filter="entrezgene_id", values=Biallelic.OMIM.gene.value, mart=ensembl)
Biallelic.OMIM.gene <- Biallelic.OMIM.gene.bm$hgnc_id %>% unique %>% c(., "HGNC:4823", "HGNC:4824") ## HBA1 and HBA2 are not labeled as AR in omim. Adding them back. 

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

Biallelic.gene.omim.ddd.clingen <- setdiff(unique(c(Biallelic.OMIM.gene, Biallelic.DDD.gene, Biallelic.clingen.gene)), "")
MonoAndBi.gene.omim.ddd.clingen <- setdiff(unique(c(MonoAndBi.OMIM.gene, MonoAndBi.DDD.gene, MonoandBiallelic.clingen.gene)), "")
Biallelic.gene.omim.ddd.clingen <- c(Biallelic.gene.omim.ddd.clingen, "HGNC:20718") ## add OTUD7A, new disese gene from this study
MonoAndBi.gene.omim.ddd.clingen <- c(MonoAndBi.gene.omim.ddd.clingen, "HGNC:20718") ## add OTUD7A, new disese gene from this study

write_lines(Biallelic.gene.omim.ddd.clingen, "./output/biallelicgenes.txt")
write_lines(MonoAndBi.gene.omim.ddd.clingen, "./output/biallelicmonoallelicgenes.txt")


## NAHR genes
## this is done after the NAHR coordinates are calculated. 
Recurrentdel <- read_tsv("./output/NAHR_GRCh38.withgenes.tsv")  %>% mutate(seqnames=gsub("chr", "",seqnames))
## genes in the 17_GL000258v2_alt haplotype do not map properly. Therefore the region needs to be converted to the chromosome 17 equivilant. 
idx.17q21.3 <- which(Recurrentdel$seqnames=="17_GL000258v2_alt" & Recurrentdel$end == 1355937)
Recurrentdel$start[idx.17q21.3] <- 45607798	
Recurrentdel$end[idx.17q21.3] <- 46210001
Recurrentdel$seqnames[idx.17q21.3] <-"17"
Recurrentdel.w.freq <- Recurrentdel %>% filter(freq>0)
Recurrentdel.w.freqCMA <- Recurrentdel %>% filter(freqCMA>0)
gr.Recurrentdel.w.freq <-  with(Recurrentdel.w.freq, GRanges(seqnames=seqnames, IRanges(start = start, end= end ), name= name))
gr.Recurrentdel.w.freqCMA <-  with(Recurrentdel.w.freqCMA, GRanges(seqnames=seqnames, IRanges(start = start, end= end ), name= name))

gene.bm <- getBM(attributes=attributes, values="*", mart=ensembl)

gr.gene.all <- with(gene.bm, GRanges( seqnames = chromosome_name,  IRanges(start=start_position, end=end_position), 
                                      hgnc_symbol = hgnc_symbol, hgnc_id = hgnc_id, 
                                      ensembl_transcript_id=ensembl_transcript_id, 
                                      contig_hgnc=paste0(chromosome_name, "_", hgnc_id))) # created this key to separate groups of genes on differnt contigs

gr.gene.grouped <- gr.gene.all %>% splitAsList(gr.gene.all$contig_hgnc) 
gr.gene.grouped.maxidx <- gr.gene.grouped %>% width %>% which.max
#save only longest transcript for each gene
gr.gene.all <- do.call(c, mclapply(1:length(gr.gene.grouped.maxidx), function(x)
{
  gr.gene.grouped[[names(gr.gene.grouped.maxidx[x])]][unname(gr.gene.grouped.maxidx[x])]
}, mc.cores= 16))

hits1 <- findOverlaps(gr.Recurrentdel.w.freq, gr.gene.all)
NAHR_HGNC_ID <- data.frame(name = gr.Recurrentdel.w.freq$name[queryHits(hits1)], hgnc_id = gr.gene.all$hgnc_id[subjectHits(hits1)], stringsAsFactors = F)
hits2 <- findOverlaps(gr.Recurrentdel.w.freqCMA, gr.gene.all)
NAHR_HGNC_ID_CMA <- data.frame(name = gr.Recurrentdel.w.freqCMA$name[queryHits(hits2)], hgnc_id = gr.gene.all$hgnc_id[subjectHits(hits2)], stringsAsFactors = F)

## top 30 most prevalent recurrent CNVs
NAHR_coding_HGNC_ID <- NAHR_HGNC_ID[NAHR_HGNC_ID$hgnc_id %in% coding_HGNC_ID, ] 
write_tsv(NAHR_coding_HGNC_ID, "./input/NAHR_coding_HGNC_ID.tsv")
write_tsv(NAHR_HGNC_ID, "./input/NAHR_HGNC_ID.tsv")

## top 51
NAHR_coding_HGNC_ID_CMA <- NAHR_HGNC_ID_CMA[NAHR_HGNC_ID_CMA$hgnc_id %in% coding_HGNC_ID, ] 
write_tsv(NAHR_coding_HGNC_ID_CMA, "./input/NAHR_coding_HGNC_ID_CMA.tsv")
write_tsv(NAHR_HGNC_ID_CMA, "./input/NAHR_HGNC_ID_CMA.tsv")
