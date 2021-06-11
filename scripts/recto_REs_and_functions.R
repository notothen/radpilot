#### Source script containing restriction enzymes,
#### size windows and functions 
#### for the RRS pilot experiment within RECTO
## 11/06/2021
## H. Christiansen
## v4.0

## This is a source script containing R objects and functions that do most of the work.
## The scripts 01_digests.R and 02_empirical_digests.R depend on this source script.

#### load packages required for these functions
#### install them first if you don't have them
library(here) # to shorten file paths
library(SimRAD) # for in silico digestion
library(bioanalyzeR) # to read bioanalyzer files
library(tidyverse) # to arrange data and plot
library(scales) # to calculate percentages

#### load restriction enzymes
#####
## SbfI : 5'--CCTGCA  GG--3
sbf1_5 <- "CCTGCA"
sbf1_3 <- "GG"
## EcoRI : 5'--G  AATTC--3
ecor1_5 <- "G"
ecor1_3 <- "AATTC"
## SphI : 5'--GCATG  C--3
sph1_5 <- "GCATG"
sph1_3 <- "C"
## PstI : 5'--CTGCA  G--3
pst1_5 <- "CTGCA"
pst1_3 <- "G"
## ApeKI : 5'--G  CWGC--3; W = A OR T
apek1_5 <- "G"
apek1_3a <- "CAGC"
apek1_3t <- "CTGC"
## MspI : 5'--C CGG--3
msp1_5 <- "C"
msp1_3 <- "CGG"
## MseI : 5'--T  TAA--3
mse1_5 <- "T"
mse1_3 <- "TAA"

## prepare data frame
run <- c(1:12)
names <- c("sbf1", "ecor1", "sph1", "pst1", "apek1", "msp1", "mse1",
           "sbf1-sph1", "sbf1-msp1", "pst1-msp1", "ecor1-sph1", "ecor1-msp1")
forward <- c(sbf1_5, ecor1_5, sph1_5, pst1_5, apek1_5, msp1_5, mse1_5,
             sbf1_5, sbf1_5, pst1_5, ecor1_5, ecor1_5)
reverse <- c(sbf1_3, ecor1_3, sph1_3, pst1_3, apek1_3a, msp1_3, mse1_3,
             sbf1_3, sbf1_3, pst1_3, ecor1_3, ecor1_3)
forward2 <- c((rep('', 4)), apek1_5, '', '',
              sph1_5, msp1_5, msp1_5, sph1_5, msp1_5)
reverse2 <- c((rep('', 4)), apek1_3t, '', '',
              sph1_3, msp1_3, msp1_3, sph1_3, msp1_3)

## combine
recto_REs <- data.frame(run, names, forward, reverse, forward2, reverse2)

## add some enzymes also as simple data frame
## needed as input for plotting functions
run <- 1
forward <- ecor1_5
reverse <- ecor1_3
ecor1 <- data.frame(run, forward, reverse)

forward <- msp1_5
reverse <- msp1_3
msp1 <- data.frame(run, forward, reverse)

forward <- pst1_5
reverse <- pst1_3
pst1 <- data.frame(run, forward, reverse)

forward <- apek1_5
reverse <- apek1_3a
apek1a <- data.frame(run, forward, reverse)

forward <- apek1_5
reverse <- apek1_3t
apek1t <- data.frame(run, forward, reverse)


## add size windows
lower_size <- c(210, 240, 0, 100, 200, 300, 400, 500, 600, 700, 800)
upper_size <- c(260, 340, 100, 200, 300, 400, 500, 600, 700, 800, 900)

## clean up
rm(apek1_3a, apek1_3t, apek1_5, ecor1_3, ecor1_5, forward, forward2, mse1_3, mse1_5, msp1_3, msp1_5,
   names, pst1_3, pst1_5, reverse, reverse2, run, sbf1_3, sbf1_5, sph1_3, sph1_5)
#####

#### create digest function
#####
recto_digest <- function(genome, enzyme, minsize, maxsize, ratio){
    dig1 <- insilico.digest(genome, enzyme$forward[1], enzyme$reverse[1], verbose = F)
    dig1_1 <- (length(dig1)-1)*ratio
    dig1_2 <- (length(size.select(dig1, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig1_3 <- (length(size.select(dig1, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig1_4 <- (length(size.select(dig1, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig1_5 <- (length(size.select(dig1, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig1_6 <- (length(size.select(dig1, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig1_7 <- (length(size.select(dig1, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig1_8 <- (length(size.select(dig1, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig1_9 <- (length(size.select(dig1, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig1_10 <- (length(size.select(dig1, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig1_11 <- (length(size.select(dig1, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig1_12 <- (length(size.select(dig1, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio
  
    
    dig2 <- insilico.digest(genome, enzyme$forward[2], enzyme$reverse[2], verbose = F)
    dig2_1 <- (length(dig2)-1)*ratio
    dig2_2 <- (length(size.select(dig2, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig2_3 <- (length(size.select(dig2, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig2_4 <- (length(size.select(dig2, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig2_5 <- (length(size.select(dig2, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig2_6 <- (length(size.select(dig2, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig2_7 <- (length(size.select(dig2, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig2_8 <- (length(size.select(dig2, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig2_9 <- (length(size.select(dig2, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig2_10 <- (length(size.select(dig2, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig2_11 <- (length(size.select(dig2, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig2_12 <- (length(size.select(dig2, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio
    
    dig3 <- insilico.digest(genome, enzyme$forward[3], enzyme$reverse[3], verbose = F)
    dig3_1 <- (length(dig3)-1)*ratio
    dig3_2 <- (length(size.select(dig3, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig3_3 <- (length(size.select(dig3, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig3_4 <- (length(size.select(dig3, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig3_5 <- (length(size.select(dig3, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig3_6 <- (length(size.select(dig3, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig3_7 <- (length(size.select(dig3, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig3_8 <- (length(size.select(dig3, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig3_9 <- (length(size.select(dig3, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig3_10 <- (length(size.select(dig3, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig3_11 <- (length(size.select(dig3, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig3_12 <- (length(size.select(dig3, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio
    
    dig4 <- insilico.digest(genome, enzyme$forward[4], enzyme$reverse[4], verbose = F)
    dig4_1 <- (length(dig4)-1)*ratio
    dig4_2 <- (length(size.select(dig4, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig4_3 <- (length(size.select(dig4, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig4_4 <- (length(size.select(dig4, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig4_5 <- (length(size.select(dig4, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig4_6 <- (length(size.select(dig4, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig4_7 <- (length(size.select(dig4, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig4_8 <- (length(size.select(dig4, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig4_9 <- (length(size.select(dig4, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig4_10 <- (length(size.select(dig4, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig4_11 <- (length(size.select(dig4, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig4_12 <- (length(size.select(dig4, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio
    
    dig5 <- insilico.digest(genome, enzyme$forward[5], enzyme$reverse[5], enzyme$forward2[5], enzyme$reverse2[5], verbose = F)
    dig5 <- adapt.select(dig5, type = "AB+BA", enzyme$forward[5], as.character(enzyme$reverse[5]), enzyme$forward2[5], as.character(enzyme$reverse2[5]))
    dig5_1 <- (length(dig5)-1)*ratio
    dig5_2 <- (length(size.select(dig5, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig5_3 <- (length(size.select(dig5, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig5_4 <- (length(size.select(dig5, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig5_5 <- (length(size.select(dig5, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig5_6 <- (length(size.select(dig5, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig5_7 <- (length(size.select(dig5, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig5_8 <- (length(size.select(dig5, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig5_9 <- (length(size.select(dig5, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig5_10 <- (length(size.select(dig5, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig5_11 <- (length(size.select(dig5, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig5_12 <- (length(size.select(dig5, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio
    
    dig6 <- insilico.digest(genome, enzyme$forward[6], enzyme$reverse[6], verbose = F)
    dig6_1 <- (length(dig6)-1)*ratio
    dig6_2 <- (length(size.select(dig6, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig6_3 <- (length(size.select(dig6, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig6_4 <- (length(size.select(dig6, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig6_5 <- (length(size.select(dig6, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig6_6 <- (length(size.select(dig6, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig6_7 <- (length(size.select(dig6, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig6_8 <- (length(size.select(dig6, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig6_9 <- (length(size.select(dig6, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig6_10 <- (length(size.select(dig6, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig6_11 <- (length(size.select(dig6, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig6_12 <- (length(size.select(dig6, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio
    
    dig7 <- insilico.digest(genome, enzyme$forward[7], enzyme$reverse[7], verbose = F)
    dig7_1 <- (length(dig7)-1)*ratio
    dig7_2 <- (length(size.select(dig7, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig7_3 <- (length(size.select(dig7, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig7_4 <- (length(size.select(dig7, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig7_5 <- (length(size.select(dig7, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig7_6 <- (length(size.select(dig7, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig7_7 <- (length(size.select(dig7, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig7_8 <- (length(size.select(dig7, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig7_9 <- (length(size.select(dig7, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig7_10 <- (length(size.select(dig7, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig7_11 <- (length(size.select(dig7, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig7_12 <- (length(size.select(dig7, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio
    
    dig8 <- insilico.digest(genome, enzyme$forward[8], enzyme$reverse[8], enzyme$forward2[8], enzyme$reverse2[8], verbose = F)
    dig8 <- adapt.select(dig8, type = "AB+BA", enzyme$forward[8], as.character(enzyme$reverse[8]), enzyme$forward2[8], as.character(enzyme$reverse2[8]))
    dig8_1 <- (length(dig8)-1)*ratio
    dig8_2 <- (length(size.select(dig8, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig8_3 <- (length(size.select(dig8, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig8_4 <- (length(size.select(dig8, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig8_5 <- (length(size.select(dig8, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig8_6 <- (length(size.select(dig8, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig8_7 <- (length(size.select(dig8, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig8_8 <- (length(size.select(dig8, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig8_9 <- (length(size.select(dig8, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig8_10 <- (length(size.select(dig8, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig8_11 <- (length(size.select(dig8, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig8_12 <- (length(size.select(dig8, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio
    
    dig9 <- insilico.digest(genome, enzyme$forward[9], enzyme$reverse[9], enzyme$forward2[9], enzyme$reverse2[9], verbose = F)
    dig9 <- adapt.select(dig9, type = "AB+BA", enzyme$forward[9], as.character(enzyme$reverse[9]), enzyme$forward2[9], as.character(enzyme$reverse2[9]))
    dig9_1 <- (length(dig9)-1)*ratio
    dig9_2 <- (length(size.select(dig9, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig9_3 <- (length(size.select(dig9, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig9_4 <- (length(size.select(dig9, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig9_5 <- (length(size.select(dig9, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig9_6 <- (length(size.select(dig9, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig9_7 <- (length(size.select(dig9, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig9_8 <- (length(size.select(dig9, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig9_9 <- (length(size.select(dig9, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig9_10 <- (length(size.select(dig9, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig9_11 <- (length(size.select(dig9, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig9_12 <- (length(size.select(dig9, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio
    
    dig10 <- insilico.digest(genome, enzyme$forward[10], enzyme$reverse[10], enzyme$forward2[10], enzyme$reverse2[10], verbose = F)
    dig10 <- adapt.select(dig10, type = "AB+BA", enzyme$forward[10], as.character(enzyme$reverse[10]), enzyme$forward2[10], as.character(enzyme$reverse2[10]))
    dig10_1 <- (length(dig10)-1)*ratio
    dig10_2 <- (length(size.select(dig10, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig10_3 <- (length(size.select(dig10, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig10_4 <- (length(size.select(dig10, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig10_5 <- (length(size.select(dig10, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig10_6 <- (length(size.select(dig10, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig10_7 <- (length(size.select(dig10, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig10_8 <- (length(size.select(dig10, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig10_9 <- (length(size.select(dig10, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig10_10 <- (length(size.select(dig10, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig10_11 <- (length(size.select(dig10, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig10_12 <- (length(size.select(dig10, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio
    
    dig11 <- insilico.digest(genome, enzyme$forward[11], enzyme$reverse[11], enzyme$forward2[11], enzyme$reverse2[11], verbose = F)
    dig11 <- adapt.select(dig11, type = "AB+BA", enzyme$forward[11], as.character(enzyme$reverse[11]), enzyme$forward2[11], as.character(enzyme$reverse2[11]))
    dig11_1 <- (length(dig11)-1)*ratio
    dig11_2 <- (length(size.select(dig11, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig11_3 <- (length(size.select(dig11, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig11_4 <- (length(size.select(dig11, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig11_5 <- (length(size.select(dig11, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig11_6 <- (length(size.select(dig11, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig11_7 <- (length(size.select(dig11, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig11_8 <- (length(size.select(dig11, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig11_9 <- (length(size.select(dig11, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig11_10 <- (length(size.select(dig11, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig11_11 <- (length(size.select(dig11, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig11_12 <- (length(size.select(dig11, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio
    
    dig12 <- insilico.digest(genome, enzyme$forward[12], enzyme$reverse[12], enzyme$forward2[12], enzyme$reverse2[12], verbose = F)
    dig12 <- adapt.select(dig12, type = "AB+BA", enzyme$forward[12], as.character(enzyme$reverse[12]), enzyme$forward2[12], as.character(enzyme$reverse2[12]))
    dig12_1 <- (length(dig12)-1)*ratio
    dig12_2 <- (length(size.select(dig12, min.size = minsize[1], max.size = maxsize[1], graph = F, verbose = F))-1)*ratio
    dig12_3 <- (length(size.select(dig12, min.size = minsize[2], max.size = maxsize[2], graph = F, verbose = F))-1)*ratio
    dig12_4 <- (length(size.select(dig12, min.size = minsize[3], max.size = maxsize[3], graph = F, verbose = F))-1)*ratio
    dig12_5 <- (length(size.select(dig12, min.size = minsize[4], max.size = maxsize[4], graph = F, verbose = F))-1)*ratio
    dig12_6 <- (length(size.select(dig12, min.size = minsize[5], max.size = maxsize[5], graph = F, verbose = F))-1)*ratio
    dig12_7 <- (length(size.select(dig12, min.size = minsize[6], max.size = maxsize[6], graph = F, verbose = F))-1)*ratio
    dig12_8 <- (length(size.select(dig12, min.size = minsize[7], max.size = maxsize[7], graph = F, verbose = F))-1)*ratio
    dig12_9 <- (length(size.select(dig12, min.size = minsize[8], max.size = maxsize[8], graph = F, verbose = F))-1)*ratio
    dig12_10 <- (length(size.select(dig12, min.size = minsize[9], max.size = maxsize[9], graph = F, verbose = F))-1)*ratio
    dig12_11 <- (length(size.select(dig12, min.size = minsize[10], max.size = maxsize[10], graph = F, verbose = F))-1)*ratio
    dig12_12 <- (length(size.select(dig12, min.size = minsize[11], max.size = maxsize[11], graph = F, verbose = F))-1)*ratio
    
    
    fragments_all <- c(dig1_1, dig2_1, dig3_1, dig4_1, dig5_1, dig6_1, dig7_1, dig8_1, dig9_1, dig10_1, dig11_1, dig12_1)
    fragments_2 <- c(dig1_2, dig2_2, dig3_2, dig4_2, dig5_2, dig6_2, dig7_2, dig8_2, dig9_2, dig10_2, dig11_2, dig12_2)
    fragments_3 <- c(dig1_3, dig2_3, dig3_3, dig4_3, dig5_3, dig6_3, dig7_3, dig8_3, dig9_3, dig10_3, dig11_3, dig12_3)
    fragments_4 <- c(dig1_4, dig2_4, dig3_4, dig4_4, dig5_4, dig6_4, dig7_4, dig8_4, dig9_4, dig10_4, dig11_4, dig12_4)
    fragments_5 <- c(dig1_5, dig2_5, dig3_5, dig4_5, dig5_5, dig6_5, dig7_5, dig8_5, dig9_5, dig10_5, dig11_5, dig12_5)
    fragments_6 <- c(dig1_6, dig2_6, dig3_6, dig4_6, dig5_6, dig6_6, dig7_6, dig8_6, dig9_6, dig10_6, dig11_6, dig12_6)
    fragments_7 <- c(dig1_7, dig2_7, dig3_7, dig4_7, dig5_7, dig6_7, dig7_7, dig8_7, dig9_7, dig10_7, dig11_7, dig12_7)
    fragments_8 <- c(dig1_8, dig2_8, dig3_8, dig4_8, dig5_8, dig6_8, dig7_8, dig8_8, dig9_8, dig10_8, dig11_8, dig12_8)
    fragments_9 <- c(dig1_9, dig2_9, dig3_9, dig4_9, dig5_9, dig6_9, dig7_9, dig8_9, dig9_9, dig10_9, dig11_9, dig12_9)
    fragments_10 <- c(dig1_10, dig2_10, dig3_10, dig4_10, dig5_10, dig6_10, dig7_10, dig8_10, dig9_10, dig10_10, dig11_10, dig12_10)
    fragments_11 <- c(dig1_11, dig2_11, dig3_11, dig4_11, dig5_11, dig6_11, dig7_11, dig8_11, dig9_11, dig10_11, dig11_11, dig12_11)
    fragments_12 <- c(dig1_12, dig2_12, dig3_12, dig4_12, dig5_12, dig6_12, dig7_12, dig8_12, dig9_12, dig10_12, dig11_12, dig12_12)
    
    x <- data.frame(enzyme$run, enzyme$names, enzyme$forward, enzyme$reverse, enzyme$forward2, enzyme$reverse2,
                    fragments_all, fragments_2, fragments_3, fragments_4, fragments_5, fragments_6, fragments_7,
                    fragments_8, fragments_9, fragments_10, fragments_11, fragments_12)
    names(x)[names(x) == "fragments_2"] <- c(paste("fragments", minsize[1], maxsize[1], sep = "_"))
    names(x)[names(x) == "fragments_3"] <- c(paste("fragments", minsize[2], maxsize[2], sep = "_"))
    names(x)[names(x) == "fragments_4"] <- c(paste("fragments", minsize[3], maxsize[3], sep = "_"))
    names(x)[names(x) == "fragments_5"] <- c(paste("fragments", minsize[4], maxsize[4], sep = "_"))
    names(x)[names(x) == "fragments_6"] <- c(paste("fragments", minsize[5], maxsize[5], sep = "_"))
    names(x)[names(x) == "fragments_7"] <- c(paste("fragments", minsize[6], maxsize[6], sep = "_"))
    names(x)[names(x) == "fragments_8"] <- c(paste("fragments", minsize[7], maxsize[7], sep = "_"))
    names(x)[names(x) == "fragments_9"] <- c(paste("fragments", minsize[8], maxsize[8], sep = "_"))
    names(x)[names(x) == "fragments_10"] <- c(paste("fragments", minsize[9], maxsize[9], sep = "_"))
    names(x)[names(x) == "fragments_11"] <- c(paste("fragments", minsize[10], maxsize[10], sep = "_"))
    names(x)[names(x) == "fragments_12"] <- c(paste("fragments", minsize[11], maxsize[11], sep = "_"))
    names(x)[names(x) == "enzyme.run"] <- "run"
    return(x)
}
#####

#### create digest plotting function
#### modified from SimRAD source code
#####
recto_digest_plot <- function(genome, enzyme, minsize, maxsize){
  sequences <- insilico.digest(genome, enzyme$forward[1], enzyme$reverse[1], verbose = F)
  ssel <- sequences[width(sequences) < maxsize & width(sequences) > minsize]
  bk <- hist(width(sequences), breaks = length(sequences)/20, plot = F)$breaks
  hist(width(sequences), border = "grey75", col = "grey75", breaks = bk, main = "", xlab = "Locus size (bp)", ylab = "Number of loci", xlim = c(0, 2200))
  hist(width(ssel), border = "red", col = "red", add = T, breaks = bk, xlim = c(0, 2200))
  text(mean(c(minsize, maxsize)), max(hist(width(ssel), breaks = bk, plot = F)$counts), pos = 4, labels = paste(length(ssel), " loci between ", minsize, " and ", maxsize, " bp", sep = ""), col = "red", cex = 0.9, font = 2)
  return(length(ssel))
}
#####

#### alternative: create a ggplot digest plotting function
#### similar to qplot.electrophoresis output
#####
## function for in silico digestion of several genomes per enzyme
## to compare with the empirical tests
recto_digest_2 <- function(genomes, enzyme){
  sequences1 <- insilico.digest(genomes[, 1], enzyme$forward[1], enzyme$reverse[1], verbose = F)
  sequences2 <- insilico.digest(genomes[, 2], enzyme$forward[1], enzyme$reverse[1], verbose = F)
  sequences3 <- insilico.digest(genomes[, 3], enzyme$forward[1], enzyme$reverse[1], verbose = F)
  seq1 <- data.frame(genome = rep(paste(names(genomes[1])), each = length(sequences1)), size = width(sequences1))
  seq2 <- data.frame(genome = rep(paste(names(genomes[2])), each = length(sequences2)), size = width(sequences2))
  seq3 <- data.frame(genome = rep(paste(names(genomes[3])), each = length(sequences3)), size = width(sequences3))
  seqwidth <- rbind(seq1, seq2, seq3)
  return(seqwidth)
}

## same function, but for digest with apek1
recto_digest_3 <- function(genomes, enzyme1, enzyme2){
  sequences1 <- insilico.digest(genomes[, 1], enzyme1$forward[1], enzyme1$reverse[1], enzyme2$forward[1], enzyme2$reverse[1], verbose = F)
  sequences2 <- insilico.digest(genomes[, 2], enzyme1$forward[1], enzyme1$reverse[1], enzyme2$forward[1], enzyme2$reverse[1], verbose = F)
  sequences3 <- insilico.digest(genomes[, 3], enzyme1$forward[1], enzyme1$reverse[1], enzyme2$forward[1], enzyme2$reverse[1], verbose = F)
  seq1 <- data.frame(genome = rep(paste(names(genomes[1])), each = length(sequences1)), size = width(sequences1))
  seq2 <- data.frame(genome = rep(paste(names(genomes[2])), each = length(sequences2)), size = width(sequences2))
  seq3 <- data.frame(genome = rep(paste(names(genomes[3])), each = length(sequences3)), size = width(sequences3))
  seqwidth <- rbind(seq1, seq2, seq3)
  return(seqwidth)
}

## same function, but for double digest
recto_doubledigest <- function(genomes, enzyme1, enzyme2){
  sequences1 <- insilico.digest(genomes[, 1], enzyme1$forward[1], enzyme1$reverse[1], enzyme2$forward[1], enzyme2$reverse[1], verbose = F)
  sequences2 <- insilico.digest(genomes[, 2], enzyme1$forward[1], enzyme1$reverse[1], enzyme2$forward[1], enzyme2$reverse[1], verbose = F)
  sequences3 <- insilico.digest(genomes[, 3], enzyme1$forward[1], enzyme1$reverse[1], enzyme2$forward[1], enzyme2$reverse[1], verbose = F)
  seq1sel <- adapt.select(sequences1, type = "AB+BA", enzyme1$forward[1], enzyme1$reverse[1], enzyme2$forward[1], enzyme2$reverse[1])
  seq2sel <- adapt.select(sequences2, type = "AB+BA", enzyme1$forward[1], enzyme1$reverse[1], enzyme2$forward[1], enzyme2$reverse[1])
  seq3sel <- adapt.select(sequences3, type = "AB+BA", enzyme1$forward[1], enzyme1$reverse[1], enzyme2$forward[1], enzyme2$reverse[1])
  seq1 <- data.frame(genome = rep(paste(names(genomes[1])), each = length(seq1sel)), size = width(seq1sel))
  seq2 <- data.frame(genome = rep(paste(names(genomes[2])), each = length(seq2sel)), size = width(seq2sel))
  seq3 <- data.frame(genome = rep(paste(names(genomes[3])), each = length(seq3sel)), size = width(seq3sel))
  seqwidth <- rbind(seq1, seq2, seq3)
  return(seqwidth)
}

## ggplot function
recto_digest_ggplot <- function(seqwidth, bins) {
  ggplot2::ggplot(data = seqwidth, aes(x = size, fill = genome)) + 
    geom_histogram(binwidth = length(seqwidth$size)/bins, alpha = 0.4, position = "identity") +
    xlim(0, 2200) +
    xlab("locus size (bp)") +
    ylab("number of loci") +
    theme_classic() +
    theme(legend.position = c(0.8, 0.8)) +
    theme(legend.title = element_blank())
}
#####

#### add also a customized qplot function
#### and a function to save outputs
#####
## version 1, all three replicates in one plot
recto_qplot_1 <- function(runsubset, ylim1, ylim2, legendpos1, legendpos2) {
  qplot.electrophoresis(runsubset, geom = "area", region.alpha = NA, facets = NULL, area.alpha = 0.4, xlim = c(0, 2200), ylim = c(ylim1, ylim2)) + 
  theme_classic() +
  theme(legend.position = c(legendpos1, legendpos2)) +
  theme(legend.title = element_blank())
}
## version 2, next to each other
recto_qplot_2 <- function(runsubset, ylim1, ylim2) {
  qplot.electrophoresis(runsubset, geom = "area", peak.fill = NA, region.alpha = NA, xlim = c(0, 2200), ylim = c(ylim1, ylim2)) + 
  theme_classic() + 
  theme(legend.position = "none") +
  theme(strip.background =  element_blank())
}
## saving output in a standardized way and multiple formats
recto_ggsave <- function(name, plot) {
  ggsave(paste0(name, ".tiff"), plot, path = here('figures'),
         width = 10, height = 5, units = 'in', dpi = 300)
  ggsave(paste0(name, ".pdf"), plot, path = here('figures'),
         width = 10, height = 5, units = 'in', dpi = 1200)
  ggsave(paste0(name, ".jpg"), plot, path = here('figures'),
         width = 10, height = 5, units = 'in', dpi = 1200)
  ggsave(paste0(name, ".eps"), plot, path = here('figures'),
         width = 10, height = 5, units = 'in', dpi = 1200)
  dev.off()
}
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
  results <- data.frame(format(fragments, big.mark = ","), genome_size, sequencer, paired_end, SNP_density,
                        format(sequenced_bases, big.mark = ","),
                        format(density, big.mark = ","), percent(portion), paste(round(everyother, digits = 0), "bp"))
  ## make the results a bit more appealing/easier to read
  names(results)[1] <- "number_of_fragments"
  names(results)[6] <- "number_of_bases_sequenced"
  names(results)[7] <- "number_of_SNPs"
  names(results)[8] <- "portion_of_genome_sequenced"
  names(results)[9] <- "one_SNP_every"
  return(results)  
}
#####
