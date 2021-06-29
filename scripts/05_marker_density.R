#### Script for calculating marker density based on
#### in silico or empirical number of loci and genome size
## 29/06/2021
## H. Christiansen
## v2.0

#### load packages
library(here) # to shorten file paths
library(scales) # to caclulate percentages
library(tidyverse) # for data manipulation
source(here("scripts/recto_REs_and_functions.R")) # custom functions

#### prepare for marker density calculation
#####
## define constant variables
reads_hiseq4000 <- 300000000 # conservative estimate
reads_hiseq2500 <- 200000000 # conservative estimate
reads_novaseq <- 1800000000 # adjust this depending on what exactly you plan to sequence on!
ind_per_lib <- 96 # adjust this depending on how many individuals you plan to pool per library!
read_length_hiseq4000 <- 150
read_length_hiseq2500 <- 125
read_length_novaseq <- 100 # adjust this depending on what exactly you plan to sequence on!
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
## based on: in silico with apek1, 200-500 bp, hiseq 4000, C. torosa genome
ostracoda1 <- marker_density(fragments = 88550, genome_size = 250000000, sequencer = "HiSeq4000", paired_end = T, 0.01)  
## based on: in silico with apek1, 200-350 bp, hiseq 2500, C. torosa genome
ostracoda2 <- marker_density(fragments = 65244, genome_size = 250000000, sequencer = "HiSeq2500", paired_end = T, 0.01)  
## based on: empirical with apek1, 200-350 bp, hiseq 2500, Macrocyprididae genome size guess
ostracoda3 <- marker_density(fragments = mean_n_loci(coverage, "Macrocyprididae"), genome_size = 250000000, sequencer = "HiSeq2500", paired_end = T, 0.01)  
ostracoda <- rbind(ostracoda1, ostracoda2, ostracoda3)
## add metadata
ostracoda <- add_column(ostracoda, input = c("in-silico_apek1_200-350", "in-silico_apek1_200-350", "empirical_apek1_200-350"), .before = 1)
ostracoda <- add_column(ostracoda, species = rep("Macrocyprididae", 3), .before = 1)
ostracoda

## malacostraca
## eusiridae
## based on: in silico with ecor1_sph1, 250-350 bp, hiseq 4000, H. azteca genome
malacostraca1 <- marker_density(fragments = 101900, genome_size = 7000000000, sequencer = "HiSeq4000", paired_end = T, 0.01)  
## based on: in silico with ecor1_sph1, 200-260 bp, hiseq 2500, H. azteca genome
malacostraca2 <- marker_density(fragments = 63572, genome_size = 7000000000, sequencer = "HiSeq2500", paired_end = T, 0.01) 
## lysianassidae
## based on: in silico with sbf1_msp1, 250-450 bp, hiseq 4000, H. azteca genome
malacostraca3 <- marker_density(fragments = 91927, genome_size = 27000000000, sequencer = "HiSeq4000", paired_end = T, 0.01)  
## based on: in silico with sbf1_msp1, 200-320 bp, hiseq 2500, H. azteca genome
malacostraca4 <- marker_density(fragments = 64094, genome_size = 27000000000, sequencer = "HiSeq2500", paired_end = T, 0.01) 
malacostraca <- rbind(malacostraca1, malacostraca2, malacostraca3, malacostraca4)
## add metadata
malacostraca <- add_column(malacostraca, input = c("in-silico_ecor1-sph1_250-350", "in-silico_ecor1-sph1_200-260", "in-silico_sbf1-msp1_250-450", "in-silico_sbf1-msp1_200-320"), .before = 1)
malacostraca <- add_column(malacostraca, species = c(rep("Eusiridae", 2), rep("Lyssianassidae", 2)), .before = 1)
malacostraca
marker_density <- rbind(ostracoda, malacostraca)
marker_density

## bivalvia
## based on: in silico with apek1, 250-350 bp, hiseq 4000, P. martensii genome
bivalvia1 <- marker_density(fragments = 102333, genome_size = 3000000000, sequencer = "HiSeq4000", paired_end = T, 0.01)  
## based on: in silico with apek1, 250-350 bp, hiseq 4000, C. gigas genome
bivalvia2 <- marker_density(fragments = 105349, genome_size = 3000000000, sequencer = "HiSeq4000", paired_end = T, 0.01)
## based on: in silico with apek1, 250-350 bp, hiseq 4000, B. platifrons genome
bivalvia3 <- marker_density(fragments = 83580, genome_size = 3000000000, sequencer = "HiSeq4000", paired_end = T, 0.01)
## based on: in silico with apek1, 250-350 bp, hiseq 4000, P. martensii genome
bivalvia4 <- marker_density(fragments = 102333, genome_size = 3000000000, sequencer = "HiSeq4000", paired_end = T, 0.01)  
## based on: in silico with apek1, 250-350 bp, hiseq 4000, C. gigas genome
bivalvia5 <- marker_density(fragments = 105349, genome_size = 3000000000, sequencer = "HiSeq4000", paired_end = T, 0.01)
## based on: in silico with apek1, 250-350 bp, hiseq 4000, B. platifrons genome
bivalvia6 <- marker_density(fragments = 83580, genome_size = 3000000000, sequencer = "HiSeq4000", paired_end = T, 0.01)

## laternula

## aequiyoldia


bivalvia <- rbind(bivalvia1, bivalvia2, bivalvia3, bivalvia4)



## add metadata
bivalvia <- add_column(bivalvia, input = c("in-silico_apek1_250-350", "in-silico_apek1_200-260", "in-silico_apek1_250-450", "in-silico_apek1_200-320"), .before = 1)
bivalvia <- add_column(bivalvia, species = c(rep("Laternula_elliptica", 2), rep("Aequiyoldia_eightsii", 2)), .before = 1)
bivalvia
marker_density <- rbind(marker_density, bivalvia)
marker_density


#####


#### read in loci statistics from empirical test libraries
#####
density_stats <- read.csv(here("data/test_libraries/density_stats.csv"), header = T)
density_stats <- as_tibble(density_stats)

densitu_stats %>%
  filter(species == "Macrocyprididae" & filtered == "no") %>%
  select(total_number_of_loci)
#####


