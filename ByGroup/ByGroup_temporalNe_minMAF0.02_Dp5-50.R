# Wild delta smelt Ne 1995-2020
# Filter SNPs and generate input file for NeEstimator to calc temporal Ne - WITH minMAF0.02, depth 5-50
# Updated Apr 2025 by MEFL

library(snpR)
library(dplyr)


#### Filter SNPs and generate input files for NeEstimator ####
setwd("ByGroup")

## Pick group of years you're targeting
#grp <- "1995-1999"
#grp <- "2002-2008"
#grp <- "2008-2020"


## Import data (genepop output from ANGSD, bamlist for sample IDs)
#geno<-read.table("BY1995-1999_minMAF0.02_minDepth5.geno",header=F)
geno <- read.table(list.files(pattern = paste0("BY", grp, ".*\\.geno$")),header=F)
ids <- read.table(paste0("ForNe_",grp,"_NoOutliers.bamlist"),header=F,stringsAsFactors=F)
ids$short <- sub(".*/(Ht[^.]+)\\..*", "\\1", ids$V1)
ids$pop <- paste0("BY",sub(".*_([0-9]{4})_.*", "\\1", ids$short))
basename <- sub("\\.geno.*", "", list.files(pattern = paste0("BY", grp, ".*\\.geno$")))

## Exclude loci on scaffolds (so only keep loci on 26 assembled chromosomes to match LD method,
#  even though this probably doesn't matter much for the temporal method)
#geno <- geno[grepl("lg", geno$V1),] #would rather keep more loci with new depth filter

#set up data in snpR format
sample_meta <- data.frame(pop=ids$pop, stringsAsFactors=F)
snp_meta <- geno[,1:2]
names(snp_meta) <- c("chr","position")
gen_data <- geno[,3:ncol(geno)]
ds <- import.snpR.data(genotypes = gen_data, 
                       snp.meta = snp_meta, 
                       sample.meta = sample_meta)

#export locus list
write.table(ds@snp.meta,file=paste0(basename,"_ANGSDloci.txt"),quote=F,row.names=F)

ds <- filter_snps(ds,
                  maf=0.02, #lowered from 0.05 at Robin's suggestion to retain more loci
                  hwe=0.99,
                  min_ind=0.75) #want this to be by facets, but it's not an option...
#NOTE: these 3 filters are the same as filters in ANGSD, but ANGSD is working with genotype
# likelihoods, so more data retained there, whereas here we have called hard genotypes and
# the filtering is redone

nrow(ds)
#Number of loci:
#1995-1999: 2568
#2002-2008: 8679
#2008-2020: 1791

## Write genepop file
#genepop format for PCA and NeEstimator:
format_snps(ds,output="genepop",outfile=paste0(basename,"_filtered.gen"),interpolate=F,facets="pop")
#save as data frame for easy import back into snpR object:
#format_snps(ds,output="NN",outfile=paste0(basename,"_filtered.geno"),interpolate=F)
#export chromosome/locus list
write.table(ds@snp.meta,file=paste0(basename,"_FilteredLoci.txt"),quote=F,row.names=F)
# #format requested by Robin:
# format_snps(ds,output="sn",outfile=paste0("BY",grp,"_filtered_012.txt"),interpolate=F,facets="pop")






#### Prep NeEstimator output for MaxTemp ####

setwd("ByGroup")

##### 1995-1999 ####
x <- read.table("tNe_BY1995-1999_in.txt",header=F)
x <- x[,c(2,4:13)]
names(x) <- c("gen1","gen2","H.mean.size","L","F","Fprime","Ne","lCIp","uCIp","lCIj","uCIj")
#add 1 to each gen
x$gen1 <- x$gen1 + 1
x$gen2 <- x$gen2 + 1
#add yr1, yr2, t, Sprime
x$yr1 <- x$gen1 + 1994
x$yr2 <- x$gen2 + 1994
x$t <- x$gen2 - x$gen1
x

#add Sprime
x$Sprime <- 1/(x$F-x$Fprime)

#change Infinite to Inf
x[] <- lapply(x, function(col) {
  if (is.character(col)) {
    col[col == "Infinite"] <- NA  # temporarily replace to allow numeric coercion
    col <- as.numeric(col)
    col[is.na(col)] <- Inf  # assign Inf to the former "Infinite" entries
  }
  return(col)
})

#rearrange
x <- x[,c("yr1","yr2","gen1","gen2","t","H.mean.size","L","F","Fprime","Ne","lCIp","uCIp","lCIj","uCIj","Sprime")]

## Export results
write.csv(x,"PreMaxTemp_BY1995-1999.csv",row.names=F,quote=F)



##### 2002-2008 ####
x <- read.table("tNe_BY2002-2008_in.txt",header=F)
x <- x[,c(2,4:13)]
names(x) <- c("gen1","gen2","H.mean.size","L","F","Fprime","Ne","lCIp","uCIp","lCIj","uCIj")
#add 1 to each gen
x$gen1 <- x$gen1 + 1
x$gen2 <- x$gen2 + 1
#add yr1, yr2, t, Sprime
x$yr1 <- x$gen1 + 2001
x$yr2 <- x$gen2 + 2001
x$t <- x$gen2 - x$gen1
x

#add Sprime
x$Sprime <- 1/(x$F-x$Fprime)

#change Infinite to Inf
x[] <- lapply(x, function(col) {
  if (is.character(col)) {
    col[col == "Infinite"] <- NA  # temporarily replace to allow numeric coercion
    col <- as.numeric(col)
    col[is.na(col)] <- Inf  # assign Inf to the former "Infinite" entries
  }
  return(col)
})

#rearrange
x <- x[,c("yr1","yr2","gen1","gen2","t","H.mean.size","L","F","Fprime","Ne","lCIp","uCIp","lCIj","uCIj","Sprime")]
x

## Export results
write.csv(x,"PreMaxTemp_BY2002-2008.csv",row.names=F,quote=F)



##### 2008-2020 ####

x <- read.table("tNe_BY2008-2017_minInd0.9_in.txt",header=F)
y <- read.table("tNe_BY2011-2020_minInd0.9_in.txt",header=F)

x <- x[,c(2,4:13)]
names(x) <- c("gen1","gen2","H.mean.size","L","F","Fprime","Ne","lCIp","uCIp","lCIj","uCIj")
#add 1 to each gen
x$gen1 <- x$gen1 + 1
x$gen2 <- x$gen2 + 1
#add yr1, yr2, t, Sprime
x$yr1 <- x$gen1 + 2007
x$yr2 <- x$gen2 + 2007
x$t <- x$gen2 - x$gen1
x

#add Sprime
x$Sprime <- 1/(x$F-x$Fprime)

#change Infinite to Inf
x[] <- lapply(x, function(col) {
  if (is.character(col)) {
    col[col == "Infinite"] <- NA  # temporarily replace to allow numeric coercion
    col <- as.numeric(col)
    col[is.na(col)] <- Inf  # assign Inf to the former "Infinite" entries
  }
  return(col)
})

#rearrange
x <- x[,c("yr1","yr2","gen1","gen2","t","H.mean.size","L","F","Fprime","Ne","lCIp","uCIp","lCIj","uCIj","Sprime")]
x


## Repeat for 2nd set of years
y <- y[,c(2,4:13)]
names(y) <- c("gen1","gen2","H.mean.size","L","F","Fprime","Ne","lCIp","uCIp","lCIj","uCIj")
#add 1 to each gen
y$gen1 <- y$gen1 + 4
y$gen2 <- y$gen2 + 4
#add yr1, yr2, t, Sprime
y$yr1 <- y$gen1 + 2007
y$yr2 <- y$gen2 + 2007
y$t <- y$gen2 - y$gen1
y

#add Sprime
y$Sprime <- 1/(y$F-y$Fprime)

#change Infinite to Inf
y[] <- lapply(y, function(col) {
  if (is.character(col)) {
    col[col == "Infinite"] <- NA  # temporarily replace to allow numeric coercion
    col <- as.numeric(col)
    col[is.na(col)] <- Inf  # assign Inf to the former "Infinite" entries
  }
  return(col)
})

#rearrange
y <- y[,c("yr1","yr2","gen1","gen2","t","H.mean.size","L","F","Fprime","Ne","lCIp","uCIp","lCIj","uCIj","Sprime")]
y



## Combine the two datasets
x <- rbind(x,y)
x <- x[!duplicated(x), ]
x <- x[order(x$yr1,x$yr2),]

## Export results
write.csv(x,"PreMaxTemp_BY2008-2020.csv",row.names=F,quote=F)




#### Prep NeEstimator output for MaxTemp - MULTIPLE methods & pcrit ####

setwd("ByGroup")

##### 1995-1999 ####
x <- read.table("tNe_BY1995-1999_in.txt",header=F)
x <- x[,c(2,4:ncol(x))]
names(x) <- c("gen1","gen2","method","pcrit","H.mean.size","L","F","Fprime","Ne","lCIp","uCIp","lCIj","uCIj")
#add 1 to each gen
x$gen1 <- x$gen1 + 1
x$gen2 <- x$gen2 + 1
#add yr1, yr2, t, Sprime
x$yr1 <- x$gen1 + 1994
x$yr2 <- x$gen2 + 1994
x$t <- x$gen2 - x$gen1
x

#add Sprime
x$Sprime <- 1/(x$F-x$Fprime)

#change Infinite to Inf
x[names(x) != "method"] <- lapply(x[names(x) != "method"], function(col) {
  if (is.character(col)) {
    col[col == "Infinite"] <- NA  # temporarily replace to allow numeric coercion
    col <- as.numeric(col)
    col[is.na(col)] <- Inf  # assign Inf to the former "Infinite" entries
  }
  return(col)
})

#rearrange
x <- x[,c("method","pcrit","yr1","yr2","gen1","gen2","t","H.mean.size","L","F","Fprime","Ne","lCIp","uCIp","lCIj","uCIj","Sprime")]

## Export results
write.csv(x,"PreMaxTemp_BY1995-1999.csv",row.names=F,quote=F)



##### 2002-2008 ####
x <- read.table("tNe_BY2002-2008_in.txt",header=F)
x <- x[,c(2,4:ncol(x))]
names(x) <- c("gen1","gen2","method","pcrit","H.mean.size","L","F","Fprime","Ne","lCIp","uCIp","lCIj","uCIj")
#add 1 to each gen
x$gen1 <- x$gen1 + 1
x$gen2 <- x$gen2 + 1
#add yr1, yr2, t, Sprime
x$yr1 <- x$gen1 + 2001
x$yr2 <- x$gen2 + 2001
x$t <- x$gen2 - x$gen1
x

#add Sprime
x$Sprime <- 1/(x$F-x$Fprime)

#change Infinite to Inf
x[names(x) != "method"] <- lapply(x[names(x) != "method"], function(col) {
  if (is.character(col)) {
    col[col == "Infinite"] <- NA  # temporarily replace to allow numeric coercion
    col <- as.numeric(col)
    col[is.na(col)] <- Inf  # assign Inf to the former "Infinite" entries
  }
  return(col)
})

#rearrange
x <- x[,c("method","pcrit","yr1","yr2","gen1","gen2","t","H.mean.size","L","F","Fprime","Ne","lCIp","uCIp","lCIj","uCIj","Sprime")]
x

## Export results
write.csv(x,"PreMaxTemp_BY2002-2008.csv",row.names=F,quote=F)



##### 2008-2020 ####
x <- read.table("tNe_BY2008-2017_in.txt",header=F)
x <- x[,c(2,4:ncol(x))]
names(x) <- c("gen1","gen2","method","pcrit","H.mean.size","L","F","Fprime","Ne","lCIp","uCIp","lCIj","uCIj")
#add 1 to each gen
x$gen1 <- x$gen1 + 1
x$gen2 <- x$gen2 + 1
#add yr1, yr2, t, Sprime
x$yr1 <- x$gen1 + 2007
x$yr2 <- x$gen2 + 2007
x$t <- x$gen2 - x$gen1
x

#add Sprime
x$Sprime <- 1/(x$F-x$Fprime)

#change Infinite to Inf
x[names(x) != "method"] <- lapply(x[names(x) != "method"], function(col) {
  if (is.character(col)) {
    col[col == "Infinite"] <- NA  # temporarily replace to allow numeric coercion
    col <- as.numeric(col)
    col[is.na(col)] <- Inf  # assign Inf to the former "Infinite" entries
  }
  return(col)
})

#rearrange
x <- x[,c("method","pcrit","yr1","yr2","gen1","gen2","t","H.mean.size","L","F","Fprime","Ne","lCIp","uCIp","lCIj","uCIj","Sprime")]
x


## Repeat for 2nd set of years
y <- read.table("tNe_BY2011-2020_in.txt",header=F)
y <- y[,c(2,4:ncol(y))]
names(y) <- c("gen1","gen2","method","pcrit","H.mean.size","L","F","Fprime","Ne","lCIp","uCIp","lCIj","uCIj")

#add 1 to each gen
y$gen1 <- y$gen1 + 4
y$gen2 <- y$gen2 + 4
#add yr1, yr2, t, Sprime
y$yr1 <- y$gen1 + 2007
y$yr2 <- y$gen2 + 2007
y$t <- y$gen2 - y$gen1
y

#add Sprime
y$Sprime <- 1/(y$F-y$Fprime)

#change Infinite to Inf
y[names(y) != "method"] <- lapply(y[names(y) != "method"], function(col) {
  if (is.character(col)) {
    col[col == "Infinite"] <- NA  # temporarily replace to allow numeric coercion
    col <- as.numeric(col)
    col[is.na(col)] <- Inf  # assign Inf to the former "Infinite" entries
  }
  return(col)
})

#rearrange
y <- y[,c("method","pcrit","yr1","yr2","gen1","gen2","t","H.mean.size","L","F","Fprime","Ne","lCIp","uCIp","lCIj","uCIj","Sprime")]
y



## Combine the two datasets
x <- rbind(x,y)
x <- x[!duplicated(x), ]
x <- x[order(x$yr1,x$yr2),]

## Export results
write.csv(x,"PreMaxTemp_BY2008-2020.csv",row.names=F,quote=F)
