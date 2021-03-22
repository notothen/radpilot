#### Script for in silico digestion of
#### various genomes with different enzymes
#### for RADseq pilot experiment RECTO
## 22/03/2021
## H. Christiansen
## v2.3

#### load packages
library(here) # to shorten file paths
library(SimRAD) # for in silico digestion
library(seqinr) # for GC calculation
library(tidyverse) # to arrange data
source(here("scripts/recto_REs_and_functions.R")) # custom functions

#### Ostracoda
#####
## reference genome C. torosa
#### the reference genomes are too big to be hosted on github
#### you need to get a local copy and store it somewhere and change the file path accordingly
Ctorosa <- ref.DNAseq(here("../refgenomes/Cyprideis_torosa.raw_filtered.fasta"),
                               subselect.contigs = F)
width(Ctorosa)
# 286.1 Mb
GC(s2c(Ctorosa))
# 0.439 GC

# digest
ostracoda_ctorosa <- recto_digest(Ctorosa, recto_REs, lower_size, upper_size, 1)
ostracoda_ctorosa$class <- "ostracoda"
ostracoda_ctorosa$ref <- "cyprideis_torosa"
write.csv(ostracoda_ctorosa, file = here("data/in_silico_results/ostracoda_ctorosa.csv"))

## simulated genome, I
simI <- sim.DNAseq(size=100000000, GCfreq=0.439)
GC(s2c(simI))

genome_size <- 100000000 # genome size: 100 Mb
ratio_simI <- genome_size/width(simI)
ratio_simI

ostracoda_simI <- recto_digest(simI, recto_REs, lower_size, upper_size, ratio_simI)
ostracoda_simI
ostracoda_simI$class <- "ostracoda"
ostracoda_simI$ref <- "simI"
write.csv(ostracoda_simI, file = here("data/in_silico_results/ostracoda_simI.csv"))

## simulated genome, II
simII <- sim.DNAseq(size=200000000, GCfreq=0.439)
GC(s2c(simI))

genome_size <- 500000000 # genome size: 500 Mb
ratio_simII <- genome_size/width(simII)
ratio_simII

ostracoda_simII <- recto_digest(simII, recto_REs, lower_size, upper_size, ratio_simII)
ostracoda_simII
ostracoda_simII$class <- "ostracoda"
ostracoda_simII$ref <- "simII"
write.csv(ostracoda_simII, file = here("data/in_silico_results/ostracoda_simII.csv"))

ostracoda <- rbind(ostracoda_ctorosa, ostracoda_simI, ostracoda_simII)
write.csv(ostracoda, file = here("data/in_silico_results/ostracoda.csv"))

###########BEGIN
#######UPDATE CODE HERE
### calculate some additional size windows
lower_size <- c(210, 240, 0, 100, 200, 300, 400, 500, 600, 700, 800)
upper_size <- c(260, 340, 100, 200, 300, 400, 500, 600, 700, 800, 900)

ostracoda_ctorosa <- recto_digest(Ctorosa, recto_REs, c(200, 250), c(350, 500), 1)
#########END


rm(Ctorosa, simI, simII, genome_size, ratio_simI, ratio_simII, ostracoda,
   ostracoda_ctorosa, ostracoda_simI, ostracoda_simII)
#####

#### Amphipoda
#####
## reference genome H. azteca
Hazteca <- ref.DNAseq(here("data/refgenomes/GCF_000764305.1_Hazt_2.0_genomic.fna"),
                      subselect.contigs = F)
width(Hazteca)
# 550.9 Mb
GC(s2c(Hazteca))
# 0.385 GC
# digest
amphipoda_hazteca <- recto_digest(Hazteca, recto_REs, lower_size, upper_size, 1)
amphipoda_hazteca$class <- "amphipoda"
amphipoda_hazteca$ref <- "Hyalella_azteca"
write.csv(amphipoda_hazteca, file = here("data/in_silico_results/amphipoda_hazteca.csv"))
rm(Hazteca)

## reference genome P. hawaiensis
Phawaiensis <- ref.DNAseq(here("refgenomes/GCA_001587735.1_Phaw3.0_genomic.fna"),
                          subselect.contigs = T, prop.contigs = 0.25)
width(Phawaiensis)
# 1000.7 Mb * 4
#GC(s2c(Phawaiensis))
# 0.408 GC
# digest
amphipoda_phaw <- recto_digest(Phawaiensis, recto_REs, lower_size, upper_size, 4)
amphipoda_phaw$class <- "amphipoda"
amphipoda_phaw$ref <- "Parhyale hawaiensis"
write.csv(amphipoda_phaw, file = here("data/in_silico_results/amphipoda_phawaiensis.csv"))
rm(Phawaiensis)

## shotgun sequences E. perdentatus
Eperdentatus <- ref.DNAseq(here("refgenomes/Eperdentatus.fasta"),
                           subselect.contigs = F)
width(Eperdentatus)
# 966.2 Mb
GC(s2c(Eperdentatus))
# 0.4408 GC 
# digest
amphipoda_eper <- recto_digest(Eperdentatus, recto_REs, lower_size, upper_size, 1)
amphipoda_eper$class <- "amphipoda"
amphipoda_eper$ref <- "Eusirus perdentatus"
write.csv(amphipoda_eper, file = here("data/in_silico_results/amphipoda_eperdentatus.csv"))
rm(Eperdentatus)

## simulated genome, I
simI <- sim.DNAseq(size=100000000, GCfreq=0.385)
GC(s2c(simI))

genome_size <- 10000000000 # genome size: 10 000 Mb
ratio_simI <- genome_size/width(simI)
amphipoda_simI <- recto_digest(simI, recto_REs, lower_size, upper_size, ratio_simI)
amphipoda_simI$class <- "amphipoda"
amphipoda_simI$ref <- "simI"
write.csv(amphipoda_simI, file = here("data/in_silico_results/amphipoda_simI.csv"))
rm(simI)

## simulated genome, II
simII <- sim.DNAseq(size=200000000, GCfreq=0.408)
GC(s2c(simII))
genome_size <- 30000000000 # genome size: 30 000 Mb
ratio_simII <- genome_size/width(simII)
amphipoda_simII <- recto_digest(simII, recto_REs, lower_size, upper_size, ratio_simII)
amphipoda_simII$class <- "amphipoda"
amphipoda_simII$ref <- "simII"
write.csv(amphipoda_simII, file = here("data/in_silico_results/amphipoda_simII.csv"))
rm(simII)

## combine 
amphipoda <- rbind(amphipoda_hazteca, amphipoda_phaw, amphipoda_eper,
                   amphipoda_simI, amphipoda_simII)
write.csv(amphipoda, file = here("data/in_silico_results/amphipoda.csv"))
rm(genome_size, ratio_simI, ratio_simII, amphipoda,
   amphipoda_hazteca, amphipoda_phaw, amphipoda_eper,
   amphipoda_simI, amphipoda_simII)
#####

#### Bivalvia
#####
## reference genome P. martensii
Pmartensii <- ref.DNAseq(here("data/refgenomes/GCA_002216045.1_PinMar1.0_genomic.fna"),
                         subselect.contigs = F)
width(Pmartensii)
# 991.0 Mb
#GC(s2c(Pmartensii))
# 0.353 GC
# digest
bivalvia_pmar <- recto_digest(Pmartensii, recto_REs, lower_size, upper_size, 1)
bivalvia_pmar$class <- "bivalvia"
bivalvia_pmar$ref <- "Pinctada martensii"
write.csv(bivalvia_pmar, file = here("data/in_silico_results/bivalvia_pmartensii.csv"))
rm(Pmartensii)

## reference genome B. platifrons
Bplatifrons <- ref.DNAseq(here("data/refgenomes/GCA_002080005.1_Bpl_v1.0_genomic.fna"),
                          subselect.contigs = F)
width(Bplatifrons)
# 1658.2 Mb * 4
#GC(s2c(Bplatifrons))
# 0.342 GC
# digest
bivalvia_bpla <- recto_digest(Bplatifrons, recto_REs, lower_size, upper_size, 4)
bivalvia_bpla$class <- "bivalvia"
bivalvia_bpla$ref <- "Bathymodiolus platifrons"
write.csv(bivalvia_bpla, file = here("data/in_silico_results/bivalvia_bplatifrons.csv"))
rm(Bplatifrons)

## reference genome C. gigas
Cgigas <- ref.DNAseq(here("data/refgenomes/GCF_000297895.1_oyster_v9_genomic.fna"),
                     subselect.contigs = F)
width(Cgigas)
# 557.7 Mb
#GC(s2c(Cgigas))
# 0.334 GC 
# digest
bivalvia_cgig <- recto_digest(Cgigas, recto_REs, lower_size, upper_size, 1)
bivalvia_cgig$class <- "bivalvia"
bivalvia_cgig$ref <- "Crassostrea gigas"
write.csv(bivalvia_cgig, file = here("data/in_silico_results/bivalvia_cgigas.csv"))
rm(Cgigas)

## simulated genome, I
simI <- sim.DNAseq(size=100000000, GCfreq=0.353)
genome_size <- 1000000000 # genome size: 1000Mb
ratio_simI <- genome_size/width(simI)
bivalvia_simI <- recto_digest(simI, recto_REs, lower_size, upper_size, ratio_simI)
bivalvia_simI$class <- "bivalvia"
bivalvia_simI$ref <- "simI"
write.csv(bivalvia_simI, file = here("data/in_silico_results/bivalvia_simI.csv"))
rm(simI)

## simulated genome, II
simII <- sim.DNAseq(size=500000000, GCfreq=0.342)
GC(s2c(simII))
genome_size <- 5000000000 # genome size: 5000Mb
ratio_simII <- genome_size/width(simII)
bivalvia_simII <- recto_digest(simII, recto_REs, lower_size, upper_size, ratio_simII)
bivalvia_simII$class <- "bivalvia"
bivalvia_simII$ref <- "simII"
write.csv(bivalvia_simII, file = here("data/in_silico_results/bivalvia_simII.csv"))
rm(simII)

## combine 
bivalvia <- rbind(bivalvia_pmar, bivalvia_bpla, bivalvia_cgig,
                  bivalvia_simI, bivalvia_simII)
write.csv(bivalvia, file = here("data/in_silico_results/bivalvia.csv"))
rm(genome_size, ratio_simI, ratio_simII, bivalvia,
   bivalvia_pmar, bivalvia_bpla, bivalvia_cgig,
   bivalvia_simI, bivalvia_simII)
#####

#### Asteroidea
#####
## reference genome A. planci
Aplanci <- ref.DNAseq(here("data/refgenomes/GCF_001949145.1_OKI-Apl_1.0_genomic.fna"),
                      subselect.contigs = F)
width(Aplanci)
# 383.9 Mb
#GC(s2c(Aplanci))
# 0.413 GC
# digest
asteroidea_apla <- recto_digest(Aplanci, recto_REs, lower_size, upper_size, 1)
asteroidea_apla$class <- "asteroidea"
asteroidea_apla$ref <- "Acanthaster planci"
write.csv(asteroidea_apla, file = here("data/in_silico_results/asteroidea_aplanci.csv"))
rm(Aplanci)

## reference genome P. miniata
Pminiata <- ref.DNAseq(here("data/refgenomes/GCA_000285935.1_Pmin_1.0_genomic.fna"),
                       subselect.contigs = F)
width(Pminiata)
# 811.0
#GC(s2c(Pminiata))
# 0.402 GC
# digest
asteroidea_pmin <- recto_digest(Pminiata, recto_REs, lower_size, upper_size, 4)
asteroidea_pmin$class <- "asteroidea"
asteroidea_pmin$ref <- "Patiria miniata"
write.csv(asteroidea_pmin, file = here("data/in_silico_results/asteroidea_pminiata.csv"))
rm(Pminiata)

## reference genome P. regularis
Pregularis <- ref.DNAseq(here("data/refgenomes/GCA_900067625.1_Patiriella_regularis_genome_assembly_1.0_genomic.fna"),
                         subselect.contigs = F)
width(Pregularis)
# 949.3 Mb
#GC(s2c(Pregularis))
# 0.404 GC 
# digest
asteroidea_preg <- recto_digest(Pregularis, recto_REs, lower_size, upper_size, 1)
asteroidea_preg$class <- "asteroidea"
asteroidea_preg$ref <- "Patririella regularis"
write.csv(asteroidea_preg, file = here("data/in_silico_results/asteroidea_pregularis.csv"))
rm(Pregularis)

## simulated genome, I
simI <- sim.DNAseq(size=100000000, GCfreq=0.413)
genome_size <- 1000000000 # genome size: 1000Mb
ratio_simI <- genome_size/width(simI)
asteroidea_simI <- recto_digest(simI, recto_REs, lower_size, upper_size, ratio_simI)
asteroidea_simI$class <- "asteroidea"
asteroidea_simI$ref <- "simI"
write.csv(asteroidea_simI, file = here("data/in_silico_results/asteroidea_simI.csv"))
rm(simI)

## simulated genome, II
simII <- sim.DNAseq(size=200000000, GCfreq=0.404)
GC(s2c(simII))
genome_size <- 2000000000 # genome size: 2000Mb
ratio_simII <- genome_size/width(simII)
asteroidea_simII <- recto_digest(simII, recto_REs, lower_size, upper_size, ratio_simII)
asteroidea_simII$class <- "asteroidea"
asteroidea_simII$ref <- "simII"
write.csv(asteroidea_simII, file = here("data/in_silico_results/asteroidea_simII.csv"))
rm(simII)

## combine 
asteroidea <- rbind(asteroidea_apla, asteroidea_pmin, asteroidea_preg,
                    asteroidea_simI, asteroidea_simII)
write.csv(asteroidea, file = here("data/in_silico_results/asteroidea.csv"))
rm(genome_size, ratio_simI, ratio_simII, asteroidea,
   asteroidea_apla, asteroidea_preg, asteroidea_pmin,
   asteroidea_simI, asteroidea_simII)
#####

#### Actinopterygii
#####
## reference genome N. coriiceps
Ncoriiceps <- ref.DNAseq(here("data/refgenomes/GCF_000735185.1_NC01_genomic.fna"),
                         subselect.contigs = F)
width(Ncoriiceps)
#  Mb
#GC(s2c(Ncoriiceps))
# 0. GC
# digest
actinopterygii_ncor <- recto_digest(Ncoriiceps, recto_REs, lower_size, upper_size, 1)
actinopterygii_ncor$class <- "actinopterygii"
actinopterygii_ncor$ref <- "Notothenia coriiceps"
write.csv(actinopterygii_ncor, file = here("data/in_silico_results/actinopterygii_ncoriiceps.csv"))
rm(Ncoriiceps)

## simulated genome, I
simI <- sim.DNAseq(size=100000000, GCfreq=0.408)
genome_size <- 1000000000 # genome size: 1000Mb
ratio_simI <- genome_size/width(simI)
actinopterygii_simI <- recto_digest(simI, recto_REs, lower_size, upper_size, ratio_simI)
actinopterygii_simI$class <- "actinopterygii"
actinopterygii_simI$ref <- "simI"
write.csv(actinopterygii_simI, file = here("data/in_silico_results/actinopterygii_simI.csv"))
rm(simI)

## simulated genome, II
simII <- sim.DNAseq(size=200000000, GCfreq=0.408)
GC(s2c(simII))
genome_size <- 1800000000 # genome size: 1800Mb
ratio_simII <- genome_size/width(simII)
actinopterygii_simII <- recto_digest(simII, recto_REs, lower_size, upper_size, ratio_simII)
actinopterygii_simII$class <- "actinopterygii"
actinopterygii_simII$ref <- "simII"
write.csv(actinopterygii_simII, file = here("data/in_silico_results/actinopterygii_simII.csv"))
rm(simII)

## combine 
actinopterygii <- rbind(actinopterygii_ncor, actinopterygii_simI, actinopterygii_simII)
write.csv(actinopterygii, file = here("data/in_silico_results/actinopterygii.csv"))
rm(genome_size, ratio_simI, ratio_simII, actinopterygii,
   actinopterygii_ncor, actinopterygii_simI, actinopterygii_simII)
#####

#### Aves
#####
## reference genome F. glacialis
Fglacialis <- ref.DNAseq(here("data/refgenomes/GCF_000690835.1_ASM69083v1_genomic.fna"),
                         subselect.contigs = F)
width(Fglacialis)
#  Mb
#GC(s2c(Fglacialis))
# 0. GC
# digest
aves_fgla <- recto_digest(Fglacialis, recto_REs, lower_size, upper_size, 1)
aves_fgla$class <- "aves"
aves_fgla$ref <- "Fulmarus glacialis"
write.csv(aves_fgla, file = here("data/in_silico_results/aves_fglacialis.csv"))
rm(Fglacialis)

## simulated genome, I
simI <- sim.DNAseq(size=100000000, GCfreq=0.412)
genome_size <- 1500000000 # genome size: 1500Mb
ratio_simI <- genome_size/width(simI)
aves_simI <- recto_digest(simI, recto_REs, lower_size, upper_size, ratio_simI)
aves_simI$class <- "aves"
aves_simI$ref <- "simI"
write.csv(aves_simI, file = here("data/in_silico_results/aves_simI.csv"))
rm(simI)

## simulated genome, II
simII <- sim.DNAseq(size=200000000, GCfreq=0.412)
GC(s2c(simII))
genome_size <- 2000000000 # genome size: 2000Mb
ratio_simII <- genome_size/width(simII)
aves_simII <- recto_digest(simII, recto_REs, lower_size, upper_size, ratio_simII)
aves_simII$class <- "aves"
aves_simII$ref <- "simII"
write.csv(aves_simII, file = here("data/in_silico_results/aves_simII.csv"))
rm(simII)

## combine 
aves <- rbind(aves_fgla, aves_simI, aves_simII)
write.csv(aves, file = here("data/in_silico_results/aves.csv"))
rm(genome_size, ratio_simI, ratio_simII, aves,
   aves_fgla, aves_simI, aves_simII)
#####

