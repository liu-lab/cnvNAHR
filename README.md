# cnvNAHR

Different sections of the scripts should be run in the following order. 

I. To preprocess the data source needed for the analyses:
  1. ./scripts/genes_inheritance.R
        This is used to define recessive and dominant gene IDs. 
  2. ./scripts/prepare_gnomadSNV.R
     ./scripts/prepare_gnomadSV.R
        These are used to prepare high-quality loss-of-function variants from gnomAD. 
  3. ./scripts/prepare_clinvar_AR.R
        This is used to prepare high-confidence pathogenic or likely pathogenic variants for recessive traits from ClinVar. 
      
II. To generate a genome-wide predicted NAHR map and to annotate the map, the following scripts can be run:
      ./scripts/NAHRmap.GRCh38.R
      ./scripts/NAHRmap.GRCh37.R
    Then, the following script is used to plot the genome-wide NAHR-deletion map:
      ./scripts/deletionmap.R
  
III. The following is the main code that produces the major analytical results of the manuscript. 
      ./scripts/carrierBurden.R