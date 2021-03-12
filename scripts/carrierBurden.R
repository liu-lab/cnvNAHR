library(dplyr)
library(readr)
library(reshape2)
library(scales)
library(parallel)
library(ggplot2)


## load gene IDs
Biallelic.gene.omim.ddd.clingen <- read_lines("./output/biallelicgenes.txt")
hgnc_table <- read_tsv("./input/hgnc_table.tsv")
NAHR_coding_HGNC_ID <- read_tsv("./input/NAHR_coding_HGNC_ID.tsv")

## load NAHR coordinates
Recurrentdel <- read_tsv("./output/NAHR_GRCh38.withgenes.tsv")  %>% mutate(seqnames=gsub("chr", "",seqnames))
## genes in the 17_GL000258v2_alt haplotype do not map properly. Therefore the region needs to be converted to the chromosome 17 equivilant. 
idx.17q21.3 <- which(Recurrentdel$seqnames=="17_GL000258v2_alt" & Recurrentdel$end == 1355937)
Recurrentdel$start[idx.17q21.3] <- 45607798	
Recurrentdel$end[idx.17q21.3] <- 46210001
Recurrentdel$seqnames[idx.17q21.3] <-"17"
Recurrentdel.w.freq <- Recurrentdel %>% filter(freq>0)

## load gene-leve ClinVar recessive variant frequencies
clinvar.biallelic <- read_tsv("./output/clinvar.recessive.table.tsv")  %>% filter(hgnc_id %in% Biallelic.gene.omim.ddd.clingen)

## load gene-level gnomad small variant frequencies
gnomadLoF.genome.table <- read_tsv("./output/gnomadLoF.genome.table.tsv")

## load gene-level gnomadSV frequencies
gnomadSV.genome.table <- read_tsv("./output/gnomadSV.genome.table.tsv")

## compute gene-level NAHR deletion frequencies
NAHR.table <- NAHR_coding_HGNC_ID %>% left_join(., hgnc_table) %>% left_join(., Recurrentdel.w.freq) %>% 
  mutate(CSRA = paste0("cnv_", name), 
         AF = freq/1e6, 
         AFR_AF = freq/1e6, 
         AMR_AF = freq/1e6, 
         EAS_AF = freq/1e6, 
         EUR_AF = freq/1e6, 
         GeneSymbol = hgnc_symbol) %>% dplyr::select(CSRA, hgnc_id, AF, AFR_AF, AMR_AF, EAS_AF, EUR_AF, GeneSymbol)


## combining all deleterious alleles except Clinvar
dz.alleles <- bind_rows(gnomadLoF.genome.table, gnomadSV.genome.table, NAHR.table)

## recessive disease alleles. Intentionally put NAHR alleles in the last for ease of following calculations
biallelic.dz.alleles <- dz.alleles %>% filter(hgnc_id %in% Biallelic.gene.omim.ddd.clingen) %>%
  bind_rows(clinvar.biallelic[!clinvar.biallelic$CSRA %in% .$CSRA, ], .) %>% filter(!is.na(AF))

biallelic.dz.alleles %>% nrow # 85063

## all biallelic disease alleles
write_tsv(biallelic.dz.alleles, "./output/all.recessive.disease.alleles.tsv")
## allele burden by gene
write_tsv(biallelic.dz.alleles %>% filter(!startsWith(CSRA, "hypomorph")) %>% group_by(GeneSymbol) %>% summarise(AFsum=sum(AF)) %>% arrange(desc(AFsum)) %>% ungroup, 
          "./output/all.recessive.gene.carrier.burden.tsv")

## known biallelic genes only in NAHR regions
biallelic.NAHR.dz.alleles <- biallelic.dz.alleles[ biallelic.dz.alleles$hgnc_id %in% NAHR_coding_HGNC_ID$hgnc_id, ]

#  all NAHR region enclosed genes- carrier freq from NAHR CNVs only
cnv.recessive.by.gene <- NAHR.table %>% 
  filter(hgnc_id %in% biallelic.NAHR.dz.alleles$hgnc_id) %>% 
  dplyr::select(GeneSymbol, AF_NAHR = AF, region = CSRA) %>% 
  mutate(region = gsub("cnv_", "", region))


## functions to calculate allele contributions to disease burdens
## to calculate fraction of NAHR contributed patients among all patients; 
## for NAHR deletions whose homozygous loss are imcompatible with live birth
calculate_Fd <- function(variants, csra_names, region_name = csra_names[which(startsWith(csra_names, "cnv_"))])
{
  ## supplement with 10 hypothetical variants each accounting for 1% of total allele frequency. 
  ## if there are hypomorphic alleles, they don't count into this total number
  variants <- c(variants, rep(sum(variants[!startsWith(csra_names, "hypomorph")])/90, 10)) 
  csra_names <- c(csra_names,paste0("hypotheticalAllele_", 1:10))
  m <- variants%*%t(variants)
  idx.allele <- which(csra_names %in% region_name)

  ## need to subtract hmz NAHR alleles because they are incompatible with live birth. 2q13 and 15q13.3 are exceptions. 
  idx.nohmz.nahr <- which(startsWith(csra_names, "cnv") & !csra_names %in% c("cnv_15q13.3_BP4-BP5", "cnv_2q13_NPHP1"))
  total.allele.burden.to.subtract <- ifelse(identical(idx.nohmz.nahr, integer(0)), 
                                            0, 
                                            sum(diag(m)[idx.nohmz.nahr]))
  ## need to subtract biallelic hypomorphic alleles since they don't cause a disease phenotype
  idx.hypomorph <- which(startsWith(csra_names, "hypomorph"))
  if(identical(idx.hypomorph, integer(0)))
  {
    hypomorph.to.subtract <- 0
  } else 
  {
    m_h <- variants[idx.hypomorph]%*%t(variants[idx.hypomorph])
    hypomorph.to.subtract <- sum(m_h)
  }
  
  total.burden <- sum(m) - total.allele.burden.to.subtract - hypomorph.to.subtract
  nahr.hmz.lethal.to.subtract <- ifelse(startsWith(region_name, "cnv_") & !region_name %in% c("cnv_15q13.3_BP4-BP5", "cnv_2q13_NPHP1"), 
                                           diag(m)[idx.allele], 
                                           0)
  if(!startsWith(region_name[1], "hypomorph"))
  {
    nahr.allele.burden <- 2 * rowSums(m)[idx.allele] - diag(m)[idx.allele] - nahr.hmz.lethal.to.subtract ## subtract the hmz NAHR lethal burden 
  } else {
    nahr.allele.burden <- 2 * sum(m[!startsWith(csra_names, "hypomorph"), idx.allele])
  }
  
  frac.nahr <- signif(nahr.allele.burden/total.burden*100, 2)
  data.frame(region = region_name, dz_NAHR_contribution_perc = frac.nahr, stringsAsFactors = F)
}


dz_burden <-function(variants, csra_names)
{
  variants <- c(variants, rep(sum(variants[!startsWith(csra_names, "hypomorph")])/90, 10)) # supplement with 10 hypothetical variants each accounting for 1% of total allele frequency. 
  m <- variants%*%t(variants)
  ## need to subtract hmz NAHR alleles because they are incompatible with live birth. 2q13 and 15q13.3 are exceptions. 
  idx.nohmz.nahr <- which(startsWith(csra_names, "cnv") & !csra_names %in% c("cnv_15q13.3_BP4-BP5", "cnv_2q13_NPHP1"))
  region_name <- csra_names[startsWith(csra_names, "cnv")]
  total.burden.to.subtract <- ifelse(identical(idx.nohmz.nahr, integer(0)), 
                                     0, 
                                     sum(diag(m)[idx.nohmz.nahr]))
  ## need to subtract biallelic hypomorphic alleles since they don't cause a disease phenotype
  idx.hypomorph <- which(startsWith(csra_names, "hypomorph"))
  if(identical(idx.hypomorph, integer(0)))
  {
    hypomorph.to.subtract <- 0
  } else 
  {
    m_h <- variants[idx.hypomorph]%*%t(variants[idx.hypomorph])
    hypomorph.to.subtract <- sum(m_h)
  }
  total.burden <- sum(m) - total.burden.to.subtract - hypomorph.to.subtract
  data.frame(region = region_name, dz_burden = total.burden, stringsAsFactors = F)
}

calculate_Fd_top3 <- function(variants, csra_names)
{
  variants <- c(variants, rep(sum(variants[!startsWith(csra_names, "hypomorph")])/90, 10)) 
  csra_names <- c(csra_names,paste0("hypotheticalAllele_", 1:10))
  top3.idx <- intersect(order(-variants), which(!startsWith(csra_names, "hypomorph")))[1:3]
  m <- variants%*%t(variants)
  ## need to subtract hmz NAHR alleles because they are incompatible with live birth. 2q13 and 15q13.3 are exceptions. 
  idx.nohmz.nahr <- which(startsWith(csra_names, "cnv") & !csra_names %in% c("cnv_15q13.3_BP4-BP5", "cnv_2q13_NPHP1"))
  region_name <- csra_names[startsWith(csra_names, "cnv")]
  idx.nahr <- which(csra_names %in% region_name)
  nahr.hmz.lethal.to.subtract <- ifelse(top3.idx %in% idx.nohmz.nahr, diag(m)[top3.idx] , 0)
  ## need to subtract biallelic hypomorphic alleles since they don't cause a disease phenotype
  idx.hypomorph <- which(startsWith(csra_names, "hypomorph"))
  if(identical(idx.hypomorph, integer(0)))
  {
    hypomorph.to.subtract <- 0
  } else 
  {
    m_h <- variants[idx.hypomorph]%*%t(variants[idx.hypomorph])
    hypomorph.to.subtract <- sum(m_h)
  }
  total.burden <- sum(m) - sum(nahr.hmz.lethal.to.subtract) - hypomorph.to.subtract
  top1.allele.burden <- 2* rowSums(m)[top3.idx[1]] - diag(m)[top3.idx[1]] - nahr.hmz.lethal.to.subtract[1]
  top2.allele.burden <- 2* rowSums(m)[top3.idx[2]] - diag(m)[top3.idx[2]] - nahr.hmz.lethal.to.subtract[2]
  top3.allele.burden <- 2* rowSums(m)[top3.idx[3]] - diag(m)[top3.idx[3]] - nahr.hmz.lethal.to.subtract[3]
  string_top_1_2_3_contribution <- paste(round(top1.allele.burden/total.burden, 2)*100, 
                                           round(top2.allele.burden/total.burden, 2)*100, 
                                           round(top3.allele.burden/total.burden, 2)*100, sep=", ")
  data.frame(region = region_name, dz_contribution_by_allele_top3 = string_top_1_2_3_contribution, stringsAsFactors = F)
}

## NAHR-deletion Impact to Recessive Disease traits (NIRD)
calculate_NIRD <- function(variants1, csra_names1, region_name1 = csra_names1[startsWith(csra_names1, "cnv")])
{
  sumsq.frac.alleles <- bind_rows(mclapply(c(csra_names1, paste0("hypotheticalAllele_", 1:10)), 
                                            function(x)
                                              {calculate_Fd(variants1, csra_names1, region_name = x)}
                                            , mc.cores = 5)) %>% 
                            arrange(desc(dz_NAHR_contribution_perc))
  ## hypomorphic alleles do not participate in the ranking
  sumsq.frac.alleles <- sumsq.frac.alleles[!startsWith(sumsq.frac.alleles$region, "hypomorph"), ]
  sumsq.all.alleles <- sum(sumsq.frac.alleles$dz_NAHR_contribution_perc)
  odds.alleles <- sumsq.frac.alleles %>% 
                    mutate(odds = dz_NAHR_contribution_perc / (sumsq.all.alleles - dz_NAHR_contribution_perc) )
  
  idx.cumsum <- which(cumsum(sumsq.frac.alleles$dz_NAHR_contribution_perc) > 0.90 * sumsq.all.alleles)[1]
  idx.allele <- which(odds.alleles$region %in% region_name1)
  odd.median <- median(odds.alleles$odds[setdiff(1:idx.cumsum, idx.allele)])
  odd.allele <- odds.alleles$odds[idx.allele]
  if(!0 %in% odd.allele)
  {
    OR <- odd.allele / odd.median
    NIRD <- signif(log2(OR), 2)
  } else
  {NIRD <- -100}
  data.frame(region = region_name1, NIRD = NIRD, stringsAsFactors = F)
}

rank_NAHR_allele <-function(variants, csra_names)
{
  variants <- c(variants, rep(sum(variants[!startsWith(csra_names, "hypomorph")])/90, 10))  
  csra_names <- c(csra_names,paste0("hypotheticalAllele_", 1:10))
  variants1 <- variants[!startsWith(csra_names, "hypomorph")]
  csra_names1 <- csra_names[!startsWith(csra_names, "hypomorph")]
  idx.nahr <- which(startsWith(csra_names1, "cnv"))
  region_name <- csra_names1[idx.nahr]
  rank.nahr <- rank(-variants1)[idx.nahr]
  data.frame(region = region_name, NAHR_allele_contribution_rank = rank.nahr, stringsAsFactors = F)
}

all_allele_burden <-function(variants, csra_names)
{
  allele.burden <- sum(variants[!startsWith(csra_names, "hypomorph")])/9*10 ## add the 10% extra supplemented hypothetical alleles
  idx.nahr <- which(startsWith(csra_names, "cnv"))
  region_name <- csra_names[idx.nahr]
  data.frame(region = region_name, all_carrier_allele_burden = allele.burden, stringsAsFactors = F)
}
  
## main function to calcualte different metrics using data from a specific population
allele_architecture <- function(AFname, frac.cutoff=0.20)
{
  biallelic.dz.alleles.pop <- biallelic.NAHR.dz.alleles %>% dplyr::select(GeneSymbol, CSRA, AFname1 = AFname)
  
  NAHR.ARgenes.freq <- biallelic.dz.alleles.pop %>% 
    group_by(GeneSymbol) %>% 
    do(calculate_Fd(.$AFname1, .$CSRA))
  
  NAHR.ARgenes.top3.freq <- biallelic.dz.alleles.pop %>% 
    group_by(GeneSymbol) %>% 
    do(calculate_Fd_top3(.$AFname1, .$CSRA))
  
  NAHR.ARgenes.OR <- biallelic.dz.alleles.pop %>% 
    group_by(GeneSymbol) %>% 
    do(calculate_NIRD(.$AFname1, .$CSRA))
  
  dz.burden <- biallelic.dz.alleles.pop %>% 
    group_by(GeneSymbol) %>% 
    do(dz_burden(.$AFname1, .$CSRA))
  
  biallelic.dz.alleles.by.gene <- biallelic.dz.alleles.pop %>% 
    group_by(GeneSymbol) %>%
    do(all_allele_burden(.$AFname1, .$CSRA))
  
  
  NAHR.ARgenes.allele.rank <- biallelic.dz.alleles.pop %>% 
                              group_by(GeneSymbol) %>% 
                              do(rank_NAHR_allele(.$AFname1, .$CSRA))
  
  
  metrics.NAHRgenes <- left_join(NAHR.ARgenes.freq, NAHR.ARgenes.top3.freq) %>%
                              left_join(., NAHR.ARgenes.OR) %>% 
                              left_join(., biallelic.dz.alleles.by.gene) %>%
                              left_join(., dz.burden) %>% 
                              left_join(., NAHR.ARgenes.allele.rank) %>%
                              mutate(region = gsub("cnv_", "", region)) %>%
                              left_join(., cnv.recessive.by.gene) %>%
                              mutate(allele_NAHR_contribution_perc = signif(AF_NAHR / all_carrier_allele_burden *100, 2)) %>%
                              arrange(desc(dz_NAHR_contribution_perc, all_carrier_allele_burden)) %>% 
                              dplyr::select(GeneSymbol, region, AF_NAHR,  
                                            all_carrier_allele_burden, 
                                            dz_burden, 
                                            dz_NAHR_contribution_perc, 
                                            allele_NAHR_contribution_perc, 
                                            NAHR_allele_contribution_rank, 
                                            dz_contribution_by_allele_top3, 
                                            NIRD)
  
  
  ## this is percent of ar genes significantly affected by NAHR for carrier frequencies. 
  percentGenesnahr <- round(length(unique(metrics.NAHRgenes$GeneSymbol[metrics.NAHRgenes$dz_NAHR_contribution_perc > frac.cutoff*100]))/length(Biallelic.gene.omim.ddd.clingen)*100, 2)
  print(percentGenesnahr)
  write_tsv(metrics.NAHRgenes, path = paste0("./output/NAHRgenesRanked_", AFname, "_", percentGenesnahr, "_percGenesOver_", frac.cutoff, "_dzContribution.tsv"))
}

allele_architecture("AF")
allele_architecture("AFR_AF")
allele_architecture("AMR_AF")
allele_architecture("EAS_AF")
allele_architecture("EUR_AF")

## make population comparison plots
metrics_AF <- read_tsv("./output/NAHRgenesRanked_AF_2.11_percGenesOver_0.2_dzContribution.tsv") %>% 
                mutate(key = paste0(GeneSymbol, " (", region, ")")) %>% 
                dplyr::select(key, pop = NIRD)
metrics_AF_AFR <- read_tsv("./output/NAHRgenesRanked_AFR_AF_2.18_percGenesOver_0.2_dzContribution.tsv") %>% 
                mutate(key = paste0(GeneSymbol, " (", region, ")")) %>% 
                dplyr::select(key, AFR = NIRD)
metrics_AF_AMR <- read_tsv("./output/NAHRgenesRanked_AMR_AF_2.3_percGenesOver_0.2_dzContribution.tsv") %>% 
                mutate(key = paste0(GeneSymbol, " (", region, ")")) %>% 
                dplyr::select(key, AMR = NIRD)
metrics_AF_EAS <- read_tsv("./output/NAHRgenesRanked_EAS_AF_2.26_percGenesOver_0.2_dzContribution.tsv") %>% 
                mutate(key = paste0(GeneSymbol, " (", region, ")")) %>% 
                dplyr::select(key, EAS = NIRD)
metrics_AF_EUR <- read_tsv("./output/NAHRgenesRanked_EUR_AF_2.11_percGenesOver_0.2_dzContribution.tsv") %>% 
                mutate(key = paste0(GeneSymbol, " (", region, ")")) %>% 
                dplyr::select(key, EUR = NIRD)

pop_comp <- left_join(metrics_AF, metrics_AF_AFR) %>% 
  left_join(., metrics_AF_AMR) %>% 
  left_join(., metrics_AF_EAS) %>% 
  left_join(., metrics_AF_EUR)

pop_comp <- pop_comp %>% 
              mutate(relative_AFR = signif(AFR - pop, 2), 
                     relative_AMR = signif(AMR - pop, 2), 
                     relative_EAS = signif(EAS - pop, 2), 
                     relative_EUR = signif(EUR - pop, 2))  %>% 
              filter(!(relative_AFR==0 & relative_AMR==0 & relative_EAS==0 & relative_EUR==0))

reshape_pop_comp <- pop_comp %>% dplyr::select(key, AFR = relative_AFR, AMR = relative_AMR, EAS = relative_EAS, EUR = relative_EUR) %>% melt
reshape_pop_comp$key <- factor(reshape_pop_comp$key, levels = rev(unique(reshape_pop_comp$key)))

plot_pop_comp <- ggplot(reshape_pop_comp, 
       aes(x=variable, y=key)) + 
  geom_tile(aes(fill = value), color = "white")+
  geom_text(aes(label = value), colour = "grey60", size=2.8) +
  scale_fill_gradientn(colors =  c("seagreen4", "white", "orangered2"), 
                       values = rescale(c(min(reshape_pop_comp$value), 0, max(reshape_pop_comp$value))), 
                       name="orange: CNV+SNV\ncyan: SNV+SNV") +
  xlab("") + ylab("") 

ggsave(plot_pop_comp, filename = "./output/pop_diff_biallelic_var_type.pdf", width = 8, height = 11)

## adopt the NIRD algorithm to calculate contribution of a given Allele's Impact Recessive Disease (AIRD). 
calculate_AIRD <- function(gene_name, allele_name, AFname = "AF")
{
  biallelic.dz.alleles.sel <- biallelic.dz.alleles %>% dplyr::select(GeneSymbol, CSRA, AFname0 = AFname) %>% filter(GeneSymbol == gene_name)
  calculate_NIRD(biallelic.dz.alleles.sel$AFname0, biallelic.dz.alleles.sel$CSRA, region_name = allele_name)
}
  
## delta F508
calculate_AIRD("CFTR", "7_117559590_ATCT_A")
calculate_AIRD("CFTR", "7_117559590_ATCT_A", AFname = "AFR_AF")
## Tay Sachs AJ 2nd highest mutation c.1421+1G>C
calculate_AIRD("HEXA", "15_72346234_C_G")
## AJ familial dysautonomia 
calculate_AIRD("ELP1", "9_108899816_A_G")
## AJ Canavan top allele
calculate_AIRD("ASPA", "17_3499000_A_C")
## AJ Fanconi anemia group C top allele
calculate_AIRD("FANCC", "9_95172033_T_A")
## AJ Tay-Sachs top allele
calculate_AIRD("HEXA", "15_72346579_G_GGATA")
## AJ NPC A top allele
calculate_AIRD("NPC1", "18_23536736_A_G")
## AJ Bloom top allele; this is absent in gnomAD
calculate_AIRD("BLM", "15_90766923_ATCTGA_TAGATTC")
## AJ Mucolipidosis IV top allele
calculate_AIRD("MCOLN1", "19_7526759_A_G")

calculate_AIRD("SMPD1", "11_6392055_TC_T")

calculate_AIRD("NPHP1", "cnv_2q13_NPHP1")
calculate_AIRD("SMN1", "SMN1_del")


## analysis for expected ratios of hmz, comhet small variants, comhet CNV + small variant
## compare them with meta-analysis observation for patient counts
calculate_exp_comhet_comhetdel_ratio <- function(gene_name, obs_hmz, obs_comhet, obs_comhetdel, AFname="AF")
{
  biallelic.dz.alleles.pop.sel <- biallelic.NAHR.dz.alleles %>% 
    filter(GeneSymbol == gene_name) %>%
    dplyr::select(GeneSymbol, CSRA, AFname1 = AFname)
  variants0 <- biallelic.dz.alleles.pop.sel$AFname1
  csra_names0 <- biallelic.dz.alleles.pop.sel$CSRA
  
  variants <- c(variants0, rep(sum(variants0[!startsWith(csra_names0, "hypomorph")])/90, 10)) 
  csra_names <- c(csra_names0,paste0("hypotheticalAllele_", 1:10))
  region_name <- csra_names[startsWith(csra_names, "cnv")]
  m <- variants%*%t(variants)
  idx.allele <- which(csra_names %in% region_name)
  
  ## need to subtract biallelic hypomorphic alleles since they don't cause a disease phenotype
  idx.hypomorph <- which(startsWith(csra_names, "hypomorph"))
  if(identical(idx.hypomorph, integer(0)))
  {
    hypomorph.to.subtract <- 0
    hypomorph.diag.to.subtract <- 0
  } else 
  {
    m_h <- variants[idx.hypomorph]%*%t(variants[idx.hypomorph])
    hypomorph.to.subtract <- 2 * sum(m_h[lower.tri(m_h)])
    hypomorph.diag.to.subtract <- sum(diag(m_h))
  }
  nahr.hmz.lethal.to.subtract <- ifelse(startsWith(region_name, "cnv_") & !region_name %in% c("cnv_15q13.3_BP4-BP5", "cnv_2q13_NPHP1"), 
                                           diag(m)[idx.allele], 
                                           0) %>% sum
  
  hmz <- sum(diag(m)) - hypomorph.diag.to.subtract - nahr.hmz.lethal.to.subtract
  comhet_plus_comhetdel <- 2 * sum(m[lower.tri(m)]) - hypomorph.to.subtract
  
  comhetdel <-sum( 2 * rowSums(m)[idx.allele] - 2* diag(m)[idx.allele] ) 
 
  comhet <- comhet_plus_comhetdel - comhetdel
  ratio_hmz_in_all = hmz/(hmz+comhet_plus_comhetdel)

  data.frame(gene_name = gene_name, 
             exp_frac_hmz_in_all = signif(ratio_hmz_in_all, 2), 
             act_frac_hmz_in_all = obs_hmz/(obs_comhet + obs_comhetdel + obs_hmz), 
             exp_ratio_comhet_comhetdel = signif(comhet/comhetdel, 2), 
             act_ratio_comhet_comhetdel = signif(obs_comhet/obs_comhetdel, 2), 
             stringsAsFactors = F)
}

metaanal <- read_tsv("./input/metaanalysis.tsv")
metaanal_exp <- bind_rows(lapply(1:nrow(metaanal), function(x){calculate_exp_comhet_comhetdel_ratio(metaanal$GeneSymbol[x], metaanal$hmz[x], metaanal$comhet[x], metaanal$comhetdel[x])}))
idx.RBM8A <- which(metaanal_exp$gene_name=="RBM8A")
## calcualte the statistics of the actual fraction of hmz relative to the expected fraction of hmz per gene
summary(metaanal_exp$act_frac_hmz_in_all[-idx.RBM8A]/metaanal_exp$exp_frac_hmz_in_all[-idx.RBM8A])  # RBM8A is excluded because it fellows a different pattern compared to the rest of the genes
## calcualte the statistics of the expected ratio of comhetdel / comhet
summary(1/metaanal_exp$exp_ratio_comhet_comhetdel[-idx.RBM8A])
## calcualte the statistics of the actual ratio of comhetdel / comhet
summary(1/metaanal_exp$act_ratio_comhet_comhetdel[-idx.RBM8A])
