# library(devtools)
# install_github("lindenb/rbcf")
library(dplyr)
library(readr)
library(rbcf)
library(parallel)


#### Download the gnomAD v3.1
chrs = c(-(-22:(-1)), "X", "Y")
# dl_gnomad = function(chr_id){
#   # vcf_file = paste0("https://storage.googleapis.com/gnomad-public/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chr", chr_id,".vcf.bgz")
#   tbi_file = paste0("https://storage.googleapis.com/gnomad-public/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chr", chr_id,".vcf.bgz.tbi")
# 
#   # system(paste0('wget ' , vcf_file -P ./input.bigfiles/gnomad))
#   system(paste0('wget ', tbi_file -P ./input.bigfiles/gnomad))
# }
# 
# parallel::mclapply(chrs,
#                    dl_gnomad,
#                    mc.cores = 12)

system("grep '^#' ./input.bigfiles/gnomad/gnomad.genomes.v3.1.sites.chr22.vcf > ./input.bigfiles/gnomad/gnomad.genomes.v3.1.sites.header.vcf")
system("sed 's/=vep/=CSQ/g' ./input.bigfiles/gnomad/gnomad.genomes.v3.1.sites.header.vcf > ./input.bigfiles/gnomad/gnomad.genomes.v3.1.sites.header.CSQ.vcf")

convert.CSQ.HC = function(chr_id)
{
  system(paste0("cp ./input.bigfiles/gnomad/gnomad.genomes.v3.1.sites.header.CSQ.vcf ./input.bigfiles/gnomad/gnomad.genomes.v3.1.sites.chr", chr_id,".CSQ.HC.vcf"))
  system(paste0("grep '|HC|||' ./input.bigfiles/gnomad/gnomad.genomes.v3.1.sites.chr", chr_id,".vcf.bgz | awk '{if($7==\"PASS\") print}' | sed 's/vep=/CSQ=/g' >> ./input.bigfiles/gnomad/gnomad.genomes.v3.1.sites.chr", chr_id,".PASS.CSQ.HC.vcf"))
}
parallel::mclapply(chrs,
                   convert.CSQ.HC,
                   mc.cores = 12)


clean.HC = function(chr_id)
{
  filename <- paste0("./input.bigfiles/gnomad/gnomad.genomes.v3.1.sites.chr", chr_id,".PASS.CSQ.HC.vcf")
  filenameout <- paste0("./input.bigfiles/gnomad/gnomad.genomes.v3.1.sites.chr", chr_id,".PASS.CSQ.HCclean.vcf")
  fp <- bcf.open(filename,FALSE)
  # error on opening (exit 0 for tests)
  if(is.null(fp)) quit(save="no",status=0,runLast=FALSE)
  # current variant
  vc <- NULL
  out <- bcf.new.writer(fp,filenameout)
  while(!is.null(vc<-bcf.next(fp))) 
  {
    df.vepparse <- variant.vep(vc)
    ## fix parsing problem in the `LoF_info"` column
    if(nrow(df.vepparse)>1)
    {
      df.vepparse.edit <- df.vepparse[1,]
      idx.vep <- 1
      while((idx.vep <- idx.vep +1) != nrow(df.vepparse)+1)
      {
        if(df.vepparse$Allele[idx.vep] == df.vepparse$Allele[1])
        {df.vepparse.edit <- rbind(df.vepparse.edit, df.vepparse[idx.vep,])} else
        {
          df.vepparse.edit$`LoF_info"`[nrow(df.vepparse.edit)] <- paste0(df.vepparse.edit$`LoF_info"`[nrow(df.vepparse.edit)], ",", df.vepparse$Allele[idx.vep])
        }
      }
    } else 
    {
      df.vepparse.edit <- df.vepparse
    }
    
    hgnc_id = df.vepparse.edit$HGNC_ID[which(df.vepparse.edit$LoF=="HC" & df.vepparse.edit$BIOTYPE=="protein_coding")] %>% unique %>% head(1)
    if(!(identical(hgnc_id, character(0)) | identical(hgnc_id, "")) )
    {
      df.vepparse.selTx <- df.vepparse.edit[df.vepparse.edit$HGNC_ID == hgnc_id &
                                              startsWith(df.vepparse.edit$Gene, "ENS") &
                                              df.vepparse.edit$BIOTYPE == "protein_coding", ]
      if(nrow(df.vepparse.selTx)>0 & variant.int.attribute(vc = vc, att = "AN")!=0)
      {
        veplogic <- !any(c(any(df.vepparse.selTx$LoF!="HC"), 
                           any(df.vepparse.selTx$LoF_filter!=""), 
                           any(df.vepparse.selTx$LoF_flags!=""), 
                           any(grepl("50_BP_RULE:FAIL", df.vepparse.selTx$`LoF_info"`))))
        if(veplogic &
           !variant.flag.attribute(vc = vc, att = "lcr") &
           variant.int.attribute(vc = vc, att = "AN") > 7.5e4 &
           variant.int.attribute(vc = vc, att = "QUALapprox") < 1e5 &
           variant.float.attribute(vc = vc, att = "AF") < 0.01 &
           variant.int.attribute(vc = vc, att = "nhomalt") == 0)
        {
          bcf.write.variant(out,vc)
        }
      }
    }
  }
  
  
  # dispose the vcf reader
  bcf.close(fp)
  # show
  bcf.close(out)
  
}

parallel::mclapply(chrs,
                   clean.HC,
                   mc.cores = 12)


extract_from_vcf = function(chr_id)
{
  filename <- paste0("./input.bigfiles/gnomad/gnomad.genomes.v3.1.sites.chr", chr_id,".PASS.CSQ.HCclean.vcf")
  fp <- bcf.open(filename,FALSE)
  # error on opening (exit 0 for tests)
  if(is.null(fp)) quit(save="no",status=0,runLast=FALSE)
  # current variant
  vc <- NULL
  df.HClof <- NULL
  while(!is.null(vc<-bcf.next(fp))) 
  {
    df.vepparse <- variant.vep(vc)
    hgnc_id = df.vepparse$HGNC_ID[which(df.vepparse$LoF=="HC" & df.vepparse$BIOTYPE=="protein_coding")] %>% unique %>% head(1)
    CSRA = paste(gsub("chr", "", variant.chrom(vc)), variant.start(vc), variant.reference(vc), variant.alt.alleles(vc), sep="_")
    AF = variant.float.attribute(vc = vc, att = "AF") 
    AFR_AF = variant.float.attribute(vc = vc, att = "AF-afr")
    AMR_AF = variant.float.attribute(vc = vc, att = "AF-amr") 
    EAS_AF = variant.float.attribute(vc = vc, att = "AF-eas")
    EUR_AF = variant.float.attribute(vc = vc, att = "AF-nfe")
    df.HClof <- rbind(df.HClof, data.frame(CSRA, AF, AFR_AF, AMR_AF, EAS_AF, EUR_AF, hgnc_id, stringsAsFactors = F))
  }
  # dispose the vcf reader
  bcf.close(fp)
  df.HClof
}


gnomadLoF.genome.table <- do.call(rbind, mclapply(chrs,
                        extract_from_vcf,
                        mc.cores = 12)) 
hgnc_table <- read_tsv("./input/hgnc_table.tsv")
gnomadLoF.genome.table <- gnomadLoF.genome.table %>% 
                            left_join(hgnc_table) %>% 
                            arrange(CSRA) %>% 
                            dplyr::rename(., GeneSymbol = hgnc_symbol)

write_tsv(gnomadLoF.genome.table, "./output/gnomadLoF.genome.table.tsv")
