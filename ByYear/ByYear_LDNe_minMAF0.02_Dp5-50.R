# Wild delta smelt Ne 1995-2020 
# Filter SNPs and generate input file for NeEstimator to calc LD Ne - WITH minMAF0.02, depth 5-100
# Updated May 2025 by MEFL

library(snpR)
library(dplyr)


#### Filter SNPs and generate input files for NeEstimator ####

setwd("ByYear")
byYear <- data.frame(file=list.files(pattern=".geno"))
byYear$pop <- substr(sub("\\..*", "", byYear$file),1,6)
byYear$yr <- as.numeric(substr(byYear$pop,3,6))
byYear$n.loci <- NA


for (i in list.files(pattern=".geno")){
  geno <- read.table(i,header=F)
  basename <- sub("\\.geno.*", "", i)
  pop <- substr(basename,1,6)
  yr <- as.numeric(substr(basename,3,6))
  
  ## Exclude loci on scaffolds (only want loci on 26 assembled chromosomes)
  geno <- geno[grepl("lg", geno$V1),]
  
  ## Set up data in snpR format
  sample_meta <- data.frame(pop=rep(pop,ncol(geno)-2), stringsAsFactors=F)
  snp_meta <- geno[,1:2]
  names(snp_meta) <- c("chr","position")
  gen_data <- geno[,3:ncol(geno)]
  x.snpR <- import.snpR.data(genotypes = gen_data, 
                             snp.meta = snp_meta, 
                             sample.meta = sample_meta)
  
  ## Filter loci
  x.snpR <- filter_snps(x.snpR,
                        maf=0.02,
                        hwe=0.99,
                        min_ind=0.75)
  
  # Add number of filtered loci to results table
  byYear$n.loci[byYear$file==i] <- nrow(x.snpR)
  
  ## Save output
  #genepop format for PCA:
  format_snps(x.snpR,output="genepop",outfile=paste0(basename,"_filtered.gen"),interpolate=F,facets="pop")
  #save as data frame for easy import back into snpR object:
  #format_snps(x.snpR,output="NN",outfile=paste0(basename,"_filtered.geno"),interpolate=F)
  #export chromosome/locus list
  write.table(x.snpR@snp.meta,file=paste0(basename,"_loci.txt"),quote=F,row.names=F)
}

#save results
byYear
write.csv(byYear,"nloci_ByYearANGSD_minMAF0.02_Dp5-50.csv",row.names=F,quote=F)
