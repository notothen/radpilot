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
## (perhaps add later a function to estimate how many inds. to pool)
reads_hiseq4000 <- 300000000
reads_hiseq2500 <- 200000000
reads_novaseq <- 10000000000 # need to check this
ind_per_lib <- 96

## these ones are important
read_length_hiseq4000 <- 150
read_length_hiseq2500 <- 125
read_length_novaseq <- 100 # need to check this
#####

#### read in some statistics from test libraries
#### and create a function to select the data per species
#####
coverage <- read.csv(here("data/test_libraries/coverage_stats.csv"), header = T)
coverage$opt_M <- as.factor(coverage$opt_M)

## function to calculate mean number of loci across individuals
mean_n_loci <- function(coverage, species_sel) {
  mean_n_loci <- coverage %>%
    filter(species == species_sel) %>%
    summarise_at(vars(n_loci), list(mean))
  return(round(deframe(mean_n_loci)))
}
#####

#### calculate marker density for different cases
#####

## ostracoda
## based on: in silico with apek1, 200-350 bp, hiseq 4000, C. torosa genome
ostracoda1 <- marker_density(fragments = 88550, genome_size = 250000000, sequencer = "HiSeq4000", paired_end = T, 0.01)  
## based on: in silico with apek1, 200-350 bp, hiseq 2500, C. torosa genome
ostracoda2 <- marker_density(fragments = 65244, genome_size = 250000000, sequencer = "HiSeq2500", paired_end = T, 0.01)  
## based on: empirical with apek1, 200-350 bp, hiseq 2500, Macrocyprididae genome size guess
ostracoda3 <- marker_density(fragments = mean_n_loci(coverage, "Macrocyprididae"), genome_size = 250000000, sequencer = "HiSeq2500", paired_end = T, 0.01)  
ostracoda <- rbind(ostracoda1, ostracoda2, ostracoda3)
ostracoda <- add_column(ostracoda, input = c("in-silico_apek1_200-350", "in-silico_apek1_200-350", "empirical_apek1_200-350"), .before = 1)
ostracoda




#####


#### read in loci statistics from empirical test libraries
#####
loci_stats <- read.csv(here("data/test_libraries/loci_stats.csv"), header = T)
loci_stats <- as_tibble(loci_stats)

loci_stats %>%
  filter(species == "Macrocyprididae" & filtered == "no") %>%
  select(number_of_loci)
#####


