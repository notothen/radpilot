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
## based on: empirical with apek1, 200-350 bp, hiseq 2500, Macrocyprididae
ostracoda3 <- marker_density(fragments = mean_n_loci(coverage, "Macrocyprididae"), genome_size = 250000000, sequencer = "HiSeq2500", paired_end = T, 0.01)  
ostracoda <- rbind(ostracoda1, ostracoda2, ostracoda3)
## add metadata
ostracoda <- add_column(ostracoda, input = c("in-silico_apek1_200-350", "in-silico_apek1_200-350", "empirical_apek1_200-350"), .before = 1)
ostracoda <- add_column(ostracoda, estimated_species = c(rep("C. torosa", 2), "Macrocyprididae"), .before = 1)
ostracoda <- add_column(ostracoda, target_class = c(rep("Ostracoda", 3)), .before = 1)
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
malacostraca <- add_column(malacostraca, estimated_species = c(rep("H. azteca (for Eusiridae)", 2), rep("H. azteca (for Lyssianassidae)", 2)), .before = 1)
malacostraca <- add_column(malacostraca, target_class = c(rep("Malacostraca", 4)), .before = 1)
malacostraca
marker_density_results <- rbind(ostracoda, malacostraca)
marker_density_results

## bivalvia
## based on: in silico with apek1, 250-350 bp, hiseq 4000, P. imbricata genome
bivalvia1 <- marker_density(fragments = 102333, genome_size = 3000000000, sequencer = "HiSeq4000", paired_end = T, 0.01)  
## based on: in silico with apek1, 250-350 bp, hiseq 4000, C. gigas genome
bivalvia2 <- marker_density(fragments = 105349, genome_size = 3000000000, sequencer = "HiSeq4000", paired_end = T, 0.01)
## based on: in silico with apek1, 250-350 bp, hiseq 4000, B. platifrons genome
bivalvia3 <- marker_density(fragments = 83580, genome_size = 3000000000, sequencer = "HiSeq4000", paired_end = T, 0.01)
## based on: in silico with apek1, 200-260 bp, hiseq 2500, P. imbricata genome
bivalvia4 <- marker_density(fragments = 64486, genome_size = 3000000000, sequencer = "HiSeq2500", paired_end = T, 0.01)  
## based on: in silico with apek1, 200-260 bp, hiseq 2500, C. gigas genome
bivalvia5 <- marker_density(fragments = 69027, genome_size = 3000000000, sequencer = "HiSeq2500", paired_end = T, 0.01)
## based on: in silico with apek1, 200-260 bp, hiseq 2500, B. platifrons genome
bivalvia6 <- marker_density(fragments = 53399, genome_size = 3000000000, sequencer = "HiSeq2500", paired_end = T, 0.01)
## based on: empirical with apek1, 200-260 bp, hiseq 2500, L. elliptica
bivalvia7 <- marker_density(fragments = mean_n_loci(coverage, "L. elliptica"), genome_size = 3000000000, sequencer = "HiSeq2500", paired_end = T, 0.01)
## based on: empirical with apek1, 200-260 bp, hiseq 2500, A. eightsii
bivalvia8 <- marker_density(fragments = mean_n_loci(coverage, "A. eightsii"), genome_size = 3000000000, sequencer = "HiSeq2500", paired_end = T, 0.01)
bivalvia <- rbind(bivalvia1, bivalvia2, bivalvia3, bivalvia4, bivalvia5, bivalvia6, bivalvia7, bivalvia8)
## add metadata
bivalvia <- add_column(bivalvia, input = c(rep("in-silico_apek1_250-350", 3), rep("in-silico_apek1_200-260", 3), rep("empirical_apek1_200-260", 2)), .before = 1)
bivalvia <- add_column(bivalvia, estimated_species = c(rep(c("P. imbricata", "C. gigas", "B. platifrons"), 2), "L. elliptica", "A. eightsii"), .before = 1)
bivalvia <- add_column(bivalvia, target_class = c(rep("Bivalvia", 8)), .before = 1)
bivalvia
marker_density_results <- rbind(marker_density_results, bivalvia)
marker_density_results

## asteroidea
## based on: in silico with apek1, 250-400 bp, hiseq 4000, A. planci genome
asteroidea1 <- marker_density(fragments = 98911, genome_size = 500000000, sequencer = "HiSeq4000", paired_end = T, 0.01)  
## based on: in silico with apek1, 250-400 bp, hiseq 4000, P. miniata genome
asteroidea2 <- marker_density(fragments = 79380, genome_size = 500000000, sequencer = "HiSeq4000", paired_end = T, 0.01)
## based on: in silico with apek1, 250-400 bp, hiseq 4000, P. regularis genome
asteroidea3 <- marker_density(fragments = 83222, genome_size = 500000000, sequencer = "HiSeq4000", paired_end = T, 0.01)
## based on: in silico with apek1, 200-300 bp, hiseq 2500, A. planci genome
asteroidea4 <- marker_density(fragments = 76988, genome_size = 500000000, sequencer = "HiSeq2500", paired_end = T, 0.01)  
## based on: in silico with apek1, 200-300 bp, hiseq 2500, P. miniata genome
asteroidea5 <- marker_density(fragments = 64466, genome_size = 500000000, sequencer = "HiSeq2500", paired_end = T, 0.01)
## based on: in silico with apek1, 200-300 bp, hiseq 2500, P. regularis genome
asteroidea6 <- marker_density(fragments = 62272, genome_size = 500000000, sequencer = "HiSeq2500", paired_end = T, 0.01)
## based on: empirical with apek1, 200-300 bp, hiseq 2500, B. loripes
asteroidea7 <- marker_density(fragments = mean_n_loci(coverage, "B. loripes"), genome_size = 500000000, sequencer = "HiSeq2500", paired_end = T, 0.01)
## based on: empirical with apek1, 200-300 bp, hiseq 2500, P. charcoti
asteroidea8 <- marker_density(fragments = mean_n_loci(coverage, "P. charcoti"), genome_size = 500000000, sequencer = "HiSeq2500", paired_end = T, 0.01)
asteroidea <- rbind(asteroidea1, asteroidea2, asteroidea3, asteroidea4, asteroidea5, asteroidea6, asteroidea7, asteroidea8)
## add metadata
asteroidea <- add_column(asteroidea, input = c(rep("in-silico_apek1_250-400", 3), rep("in-silico_apek1_200-300", 3), rep("empirical_apek1_200-300", 2)), .before = 1)
asteroidea <- add_column(asteroidea, estimated_species = c(rep(c("A. planci", "P. miniata", "P. regularis"), 2), "B. loripes", "P. charcoti"), .before = 1)
asteroidea <- add_column(asteroidea, target_class = c(rep("Asteroidea", 8)), .before = 1)
asteroidea
marker_density_results <- rbind(marker_density_results, asteroidea)
marker_density_results

## actinopterygii
## based on: in silico with ecor1-msp1, 200-600 bp, hiseq 4000, N. coriiceps genome
actinopterygii1 <- marker_density(fragments = 101138, genome_size = 1500000000, sequencer = "HiSeq4000", paired_end = T, 0.01)  
## based on: in silico with ecor1-msp1, 200-450 bp, hiseq 2500, N. coriiceps genome
actinopterygii2 <- marker_density(fragments = 81605, genome_size = 1500000000, sequencer = "HiSeq2500", paired_end = T, 0.01)  
## based on: empirical with ecor1-msp1, 200-450 bp, hiseq 2500, T. bernacchii
actinopterygii3 <- marker_density(fragments = mean_n_loci(coverage, "T. bernacchii"), genome_size = 1500000000, sequencer = "HiSeq2500", paired_end = T, 0.01)  
## based on: empirical with ecor1-msp1, 200-450 bp, hiseq 2500, T. loennbergii
actinopterygii4 <- marker_density(fragments = mean_n_loci(coverage, "T. loennbergii"), genome_size = 1500000000, sequencer = "HiSeq2500", paired_end = T, 0.01)  
actinopterygii <- rbind(actinopterygii1, actinopterygii2, actinopterygii3, actinopterygii4)
## add metadata
actinopterygii <- add_column(actinopterygii, input = c("in-silico_ecor1-msp1_200-600", "in-silico_ecor1-msp1_200-450", "empirical_ecor1-msp1_200-450", "empirical_ecor1-msp1_200-450"), .before = 1)
actinopterygii <- add_column(actinopterygii, estimated_species = c(rep("N. coriiceps", 2), "T. bernacchii", "T. loennbergii"), .before = 1)
actinopterygii <- add_column(actinopterygii, target_class = c(rep("Actinopterygii", 4)), .before = 1)
actinopterygii
marker_density_results <- rbind(marker_density_results, actinopterygii)
marker_density_results

## aves
## based on: in silico with pst1, 250-400 bp, hiseq 4000, F. glacialis genome
aves1 <- marker_density(fragments = 92422, genome_size = 1500000000, sequencer = "HiSeq4000", paired_end = T, 0.01)  
## based on: in silico with pst1, 200-300 bp, hiseq 2500, F. glacialis genome
aves2 <- marker_density(fragments = 66258, genome_size = 1500000000, sequencer = "HiSeq2500", paired_end = T, 0.01)  
## based on: empirical with pst1, 200-300 bp, hiseq 2500, P. nivea
aves3 <- marker_density(fragments = mean_n_loci(coverage, "P. nivea"), genome_size = 1500000000, sequencer = "HiSeq2500", paired_end = T, 0.01)  
aves <- rbind(aves1, aves2, aves3)
## add metadata
aves <- add_column(aves, input = c("in-silico_pst1_250-400", "in-silico_pst1_200-300", "empirical_pst1_200-300"), .before = 1)
aves <- add_column(aves, estimated_species = c(rep("F. glacialis", 2), "P. nivea"), .before = 1)
aves <- add_column(aves, target_class = c(rep("Aves", 3)), .before = 1)
aves
marker_density_results <- rbind(marker_density_results, aves)
marker_density_results

## export results
write.csv(marker_density_results, file = here("data/marker_density_results.csv"))

## clean up
rm(actinopterygii, actinopterygii1, actinopterygii2, actinopterygii3, actinopterygii4, asteroidea, asteroidea1,
   asteroidea2, asteroidea3, asteroidea4, asteroidea5, asteroidea6, asteroidea7, asteroidea8, aves, aves1, aves2, aves3,
   bivalvia, bivalvia1, bivalvia2, bivalvia3, bivalvia4, bivalvia5, bivalvia6, bivalvia7, bivalvia8, malacostraca,
   malacostraca1, malacostraca2, malacostraca3, malacostraca4, coverage, marker_density_results, ostracoda, ostracoda1, ostracoda2,
   ostracoda3, ind_per_lib, read_length_hiseq2500, read_length_hiseq4000, read_length_novaseq, reads_hiseq2500, reads_hiseq4000, 
   reads_novaseq, mean_n_loci)
#####


