#### Script for calculating marker density based on
#### in silico or empirical number of loci and genome size
## 09/06/2021
## H. Christiansen
## v1.0

#### load packages
library(here) # to shorten file paths



#library(bioanalyzeR) # to read bioanalyzer files
#library(tidyverse) # to arrange data
#library(gridExtra) # to combine plots
#library(SimRAD) # for in silico digestion
#source(here("scripts/recto_REs_and_functions.R")) # custom functions

#### prepare for marker density calculation
#### and re-map to the actual sample names
#####
## define constant variables
reads_hiseq4000 <- 300000000
reads_hiseq2500 <- 200000000
reads_novaseq <- 10000000000 # need to check this
  
ind_per_lib <- 96
read_length_hiseq4000 <- 150
read_length_hiseq2500 <- 125
read_length_novaseq <- 100 # need to check this
#####

#### function to calculate marker density
#####
marker_density <- function(fragments, genome_size, sequencer = "HiSeq2500", paired_end, SNP_density){
  ## check whether paired end or single end sequences should be estimated
  if (paired_end == T) {
    ratio <- 2
  } else if (paired_end == F) {
    ratio <- 1
  } else {
    stop("paired_end must be either true or false")
  }
  ## define read length
  if (sequencer == "HiSeq4000") {
    read_length <- read_length_hiseq4000
  } else if (sequencer == "HiSeq2500") {
    read_length <- read_length_hiseq2500
  } else if (sequencer == "NovaSeq") {
    read_length <- read_length_novaseq
  } else {
    stop("sequencer must be either HiSeq2500, HiSeq4000 or NovaSeq")
  }
  results <- data.frame()
  sequenced_bases <- ratio*read_length*fragments
  density <- sequenced_bases*SNP_density
  portion <- round((sequenced_bases/genome_size*100), digits = 2)
  everyother <- round((genome_size/density), digits = 0)
  results <- data.frame(fragments, sequenced_bases, density, portion, everyother)
  ## make the results a bit more appealing/easier to read
  names(results)[1] <- "number_of_fragments"
  names(results)[2] <- "number_of_bases_sequenced"
  names(results)[3] <- "marker_density"
  names(results)[4] <- "portion_of_genome_sequenced"
  names(results)[5] <- "one_snp_every_other"
  return(results)  
}

test <- marker_density(fragments = 88550, genome_size = 250000000, sequencer = "HiSeq4000", paired_end = T, 0.01)  
View(test)
