#### Script for calculating marker density based on
#### in silico or empirical number of loci and genome size
## 09/06/2021
## H. Christiansen
## v1.0

#### load packages
library(here) # to shorten file paths
library(scales) # to caclulate percentages
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

#### calculate marker density for different cases
#####

## ostracoda
## based on: in silico with apek1, 200-350 bp, hiseq 4000, C. torosa genome
ostracoda1 <- marker_density(fragments = 88550, genome_size = 250000000, sequencer = "HiSeq4000", paired_end = T, 0.01)  
## based on: in silico with apek1, 200-350 bp, hiseq 2500, C. torosa genome
ostracoda2 <- marker_density(fragments = 65244, genome_size = 250000000, sequencer = "HiSeq2500", paired_end = T, 0.01)  
## based on: empirical with apek1, 200-350 bp, hiseq 2500, Macrocyprididae genome size guess
ostracoda3 <- marker_density(fragments = 69817, genome_size = 250000000, sequencer = "HiSeq2500", paired_end = T, 0.01)  
ostracoda <- rbind(ostracoda1, ostracoda2, ostracoda3)


