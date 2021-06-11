#### Script for calculating marker density based on
#### in silico or empirical number of loci and genome size
## 09/06/2021
## H. Christiansen
## v1.0

#### load packages
library(here) # to shorten file paths
library(scales) # to display percentages
source(here("scripts/recto_REs_and_functions.R")) # custom functions

#### prepare for marker density calculation
#####
## define constant variables
## (these ones are probably unnecessary?)
reads_hiseq4000 <- 300000000
reads_hiseq2500 <- 200000000
reads_novaseq <- 10000000000 # need to check this
ind_per_lib <- 96

## these ones are important
read_length_hiseq4000 <- 150
read_length_hiseq2500 <- 125
read_length_novaseq <- 100 # need to check this
#####

#### function to calculate marker density
#### this function calculates how many bp are sequenced in one individual,
#### provided the number of fragments and assuming that all of these fragments
#### are sequenced in full and retained through the bioinformatic processing
#### (that means in theory each fragment must receive at least 1x sequencing coverage,
####  but more realistically minimum 10x coverage or so to pass bioinformatics)
#### this also assumes that each fragment is at least as large as two times the read-length

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
  portion <- round((sequenced_bases/genome_size), digits = 2)
  everyother <- genome_size/density
  results <- data.frame(format(fragments, big.mark = ","), format(sequenced_bases, big.mark = ","),
                        format(density, big.mark = ","), percent(portion), paste(round(everyother, digits = 0), "bp"))
  ## make the results a bit more appealing/easier to read
  names(results)[1] <- "number_of_fragments"
  names(results)[2] <- "number_of_bases_sequenced"
  names(results)[3] <- "number_of_SNPs"
  names(results)[4] <- "portion_of_genome_sequenced"
  names(results)[5] <- "one_SNP_every"
  return(results)  
}
#####

test <- marker_density(fragments = 88550, genome_size = 250000000, sequencer = "HiSeq4000", paired_end = T, 0.01)  
View(test)
