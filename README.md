# Sequencing individual genomes with recurrent deletions reveals allelic architecture and disease loci for autosomal recessive traits

The preprint of the manucript is available at https://doi.org/10.1101/2021.02.16.21251842. 

The work was selected as a platform presentation at the American Socienty of Human Genetics 2021 Annual Meeting. Below is the presentation. 

https://user-images.githubusercontent.com/43894048/138635349-0eeaa189-47bd-4f24-a8bd-1d48d8882238.mp4

A flowchart of the computational modeling processes described in this GitHub page is illustrated below.  
![flowchart](https://user-images.githubusercontent.com/43894048/138635674-02728f38-348f-4a7c-8e07-419d894e27b4.jpg)

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
      
      
The folder ./input.bigfiles is not included in the current GitHub page due to large file sizes. Files in this folder can be either downloaded from publicly available websites or generated as intermediate files during the analysis as described in the code above. 
  1. gnomAD vcf files: available for download at https://gnomad.broadinstitute.org/downloads. 
  2. ClinVar vcf files, variant summary files, and submission summary files at available at https://ftp.ncbi.nlm.nih.gov/pub/clinvar/. The data release of 01/19/2021 is used for the manuscript. 
  3. annotated vcf files are generated as intermediate files as described in the code above. 
