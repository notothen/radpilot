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

## calculate some additional size windows
lower_size <- c(250, 250, 300, 350, 250, 250, 250, 250, 200, 200, 200)
upper_size <- c(350, 400, 400, 450, 340, 330, 320, 300, 250, 260, 280)

ostracoda_ctorosa <- recto_digest(Ctorosa, recto_REs, lower_size, upper_size, 1)
ostracoda_ctorosa$class <- "ostracoda"
ostracoda_ctorosa$ref <- "cyprideis_torosa"
ostracoda_simI <- recto_digest(simI, recto_REs, lower_size, upper_size, ratio_simI)
ostracoda_simI$class <- "ostracoda"
ostracoda_simI$ref <- "simI"
ostracoda_simII <- recto_digest(simII, recto_REs, lower_size, upper_size, ratio_simII)
ostracoda_simII$class <- "ostracoda"
ostracoda_simII$ref <- "simII"
ostracoda2 <- rbind(ostracoda_ctorosa, ostracoda_simI, ostracoda_simII)
write.csv(ostracoda, file = here("data/in_silico_results/ostracoda2.csv"))

## and even more size windows
lower_size <- c(200, 200, 200, 200, 200, 200, 250, 250, 250, 200, 300)
upper_size <- c(300, 350, 400, 450, 500, 320, 420, 450, 500, 500, 500)

ostracoda_ctorosa <- recto_digest(Ctorosa, recto_REs, lower_size, upper_size, 1)
ostracoda_ctorosa$class <- "ostracoda"
ostracoda_ctorosa$ref <- "cyprideis_torosa"
ostracoda_simI <- recto_digest(simI, recto_REs, lower_size, upper_size, ratio_simI)
ostracoda_simI$class <- "ostracoda"
ostracoda_simI$ref <- "simI"
ostracoda_simII <- recto_digest(simII, recto_REs, lower_size, upper_size, ratio_simII)
ostracoda_simII$class <- "ostracoda"
ostracoda_simII$ref <- "simII"
ostracoda3 <- rbind(ostracoda_ctorosa, ostracoda_simI, ostracoda_simII)
write.csv(ostracoda, file = here("data/in_silico_results/ostracoda3.csv"))

## clean up
rm(Ctorosa, simI, simII, genome_size, ratio_simI, ratio_simII, ostracoda,
   ostracoda2, ostracoda3, ostracoda_ctorosa, ostracoda_simI, ostracoda_simII)
#####

#### Amphipoda
#####
## reference genome H. azteca
#### the reference genomes are too big to be hosted on github
#### you need to get a local copy and store it somewhere and change the file paths accordingly
Hazteca <- ref.DNAseq(here("../refgenomes/GCF_000764305.1_Hazt_2.0_genomic.fna"),
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
#### the reference genomes are too big to be hosted on github
#### you need to get a local copy and store it somewhere and change the file paths accordingly
Phawaiensis <- ref.DNAseq(here("../refgenomes/GCA_001587735.1_Phaw3.0_genomic.fna"),
                                  subselect.contigs = T, prop.contigs = 0.1)
width(Phawaiensis)
# 1000.7 Mb * 10
#GC(s2c(Phawaiensis))
# 0.408 GC
# digest
amphipoda_phaw <- recto_digest(Phawaiensis, recto_REs, lower_size, upper_size, 10)
amphipoda_phaw$class <- "amphipoda"
amphipoda_phaw$ref <- "Parhyale hawaiensis"
write.csv(amphipoda_phaw, file = here("data/in_silico_results/amphipoda_phawaiensis.csv"))
rm(Phawaiensis)

## shotgun sequences E. pontomedon
#Epontomedon <- ref.DNAseq(here("../refgenomes/Epontomedon.fasta"),
#                           subselect.contigs = F)
## this was only some shotgut sequencing data
## which is not available anymore
width(Epontomedon)
# 966.2 Mb
GC(s2c(Epontomedon))
# 0.4408 GC 
# digest
amphipoda_epon <- recto_digest(Epontomedon, recto_REs, lower_size, upper_size, 1)
amphipoda_epon$class <- "amphipoda"
amphipoda_epon$ref <- "Eusirus pontomedon"
write.csv(amphipoda_epon, file = here("data/in_silico_results/amphipoda_epontomedon.csv"))
rm(Epontomedon)

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
#amphipoda <- rbind(amphipoda_hazteca, amphipoda_phaw, amphipoda_epon,
#                  amphipoda_simI, amphipoda_simII)
amphipoda <- rbind(amphipoda_hazteca, amphipoda_phaw,
                   amphipoda_simI, amphipoda_simII)
write.csv(amphipoda, file = here("data/in_silico_results/amphipoda.csv"))

## calculate some additional size windows
lower_size <- c(250, 250, 300, 350, 250, 250, 250, 250, 200, 200, 200)
upper_size <- c(350, 400, 400, 450, 340, 330, 320, 300, 250, 260, 280)

Hazteca <- ref.DNAseq(here("../refgenomes/GCF_000764305.1_Hazt_2.0_genomic.fna"),
                      subselect.contigs = F)
amphipoda_hazteca <- recto_digest(Hazteca, recto_REs, lower_size, upper_size, 1)
amphipoda_hazteca$class <- "amphipoda"
amphipoda_hazteca$ref <- "Hyalella_azteca"
rm(Hazteca)
Phawaiensis <- ref.DNAseq(here("../refgenomes/GCA_001587735.1_Phaw3.0_genomic.fna"),
                          subselect.contigs = T, prop.contigs = 0.1)
amphipoda_phaw <- recto_digest(Phawaiensis, recto_REs, lower_size, upper_size, 4)
amphipoda_phaw$class <- "amphipoda"
amphipoda_phaw$ref <- "Parhyale hawaiensis"
rm(Phawaiensis)
simI <- sim.DNAseq(size=100000000, GCfreq=0.385)
genome_size <- 10000000000 # genome size: 10 000 Mb
ratio_simI <- genome_size/width(simI)
amphipoda_simI <- recto_digest(simI, recto_REs, lower_size, upper_size, ratio_simI)
amphipoda_simI$class <- "amphipoda"
amphipoda_simI$ref <- "simI"
rm(simI)
simII <- sim.DNAseq(size=200000000, GCfreq=0.408)
genome_size <- 30000000000 # genome size: 30 000 Mb
ratio_simII <- genome_size/width(simII)
amphipoda_simII <- recto_digest(simII, recto_REs, lower_size, upper_size, ratio_simII)
amphipoda_simII$class <- "amphipoda"
amphipoda_simII$ref <- "simII"
rm(simII)
amphipoda2 <- rbind(amphipoda_hazteca, amphipoda_phaw,
                   amphipoda_simI, amphipoda_simII)
write.csv(amphipoda2, file = here("data/in_silico_results/amphipoda2.csv"))

## and even more size windows
lower_size <- c(200, 200, 200, 200, 200, 200, 250, 250, 250, 200, 300)
upper_size <- c(300, 350, 400, 450, 500, 320, 420, 450, 500, 500, 500)

Hazteca <- ref.DNAseq(here("../refgenomes/GCF_000764305.1_Hazt_2.0_genomic.fna"),
                      subselect.contigs = F)
amphipoda_hazteca <- recto_digest(Hazteca, recto_REs, lower_size, upper_size, 1)
amphipoda_hazteca$class <- "amphipoda"
amphipoda_hazteca$ref <- "Hyalella_azteca"
rm(Hazteca)
Phawaiensis <- ref.DNAseq(here("../refgenomes/GCA_001587735.1_Phaw3.0_genomic.fna"),
                          subselect.contigs = T, prop.contigs = 0.1)
amphipoda_phaw <- recto_digest(Phawaiensis, recto_REs, lower_size, upper_size, 4)
amphipoda_phaw$class <- "amphipoda"
amphipoda_phaw$ref <- "Parhyale hawaiensis"
rm(Phawaiensis)
simI <- sim.DNAseq(size=100000000, GCfreq=0.385)
genome_size <- 10000000000 # genome size: 10 000 Mb
ratio_simI <- genome_size/width(simI)
amphipoda_simI <- recto_digest(simI, recto_REs, lower_size, upper_size, ratio_simI)
amphipoda_simI$class <- "amphipoda"
amphipoda_simI$ref <- "simI"
rm(simI)
simII <- sim.DNAseq(size=200000000, GCfreq=0.408)
genome_size <- 30000000000 # genome size: 30 000 Mb
ratio_simII <- genome_size/width(simII)
amphipoda_simII <- recto_digest(simII, recto_REs, lower_size, upper_size, ratio_simII)
amphipoda_simII$class <- "amphipoda"
amphipoda_simII$ref <- "simII"
rm(simII)
amphipoda3 <- rbind(amphipoda_hazteca, amphipoda_phaw,
                    amphipoda_simI, amphipoda_simII)
write.csv(amphipoda2, file = here("data/in_silico_results/amphipoda3.csv"))

## clean up
rm(genome_size, ratio_simI, ratio_simII, amphipoda, amphipoda2, amphipoda3,
   amphipoda_hazteca, amphipoda_phaw, amphipoda_eper,
   amphipoda_simI, amphipoda_simII)
#####

#### Bivalvia
#####
## reference genome P. martensii
#### the reference genomes are too big to be hosted on github
#### you need to get a local copy and store it somewhere and change the file paths accordingly
Pimbricata <- ref.DNAseq(here("../refgenomes/GCA_002216045.1_PinMar1.0_genomic.fna"),
                                 subselect.contigs = T, prop.contigs = 0.5)
width(Pimbricata)
# 446 Mb (223 * 2)
#GC(s2c(Pimbricata))
# 0.353 GC
# digest
bivalvia_pimb <- recto_digest(Pimbricata, recto_REs, lower_size, upper_size, 2)
bivalvia_pimb$class <- "bivalvia"
bivalvia_pimb$ref <- "Pinctada imbricata"
write.csv(bivalvia_pimb, file = here("data/in_silico_results/bivalvia_pimbricata.csv"))
rm(Pimbricata)

## reference genome B. platifrons
#### the reference genomes are too big to be hosted on github
#### you need to get a local copy and store it somewhere and change the file paths accordingly
Bplatifrons <- ref.DNAseq(here("../refgenomes/GCA_002080005.1_Bpl_v1.0_genomic.fna"),
                                       subselect.contigs = T, prop.contigs = 0.25)
width(Bplatifrons)
# 1622.4 Mb (405.6 * 4)
#GC(s2c(Bplatifrons))
# 0.342 GC
# digest
bivalvia_bpla <- recto_digest(Bplatifrons, recto_REs, lower_size, upper_size, 4)
bivalvia_bpla$class <- "bivalvia"
bivalvia_bpla$ref <- "Bathymodiolus platifrons"
write.csv(bivalvia_bpla, file = here("data/in_silico_results/bivalvia_bplatifrons.csv"))
rm(Bplatifrons)

## reference genome C. gigas
#### the reference genomes are too big to be hosted on github
#### you need to get a local copy and store it somewhere and change the file paths accordingly
Cgigas <- ref.DNAseq(here("../refgenomes/GCF_000297895.1_oyster_v9_genomic.fna"),
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
simI <- sim.DNAseq(size=10000000, GCfreq=0.353)
genome_size <- 1000000000 # genome size: 1000Mb
ratio_simI <- genome_size/width(simI)
bivalvia_simI <- recto_digest(simI, recto_REs, lower_size, upper_size, ratio_simI)
bivalvia_simI$class <- "bivalvia"
bivalvia_simI$ref <- "simI"
write.csv(bivalvia_simI, file = here("data/in_silico_results/bivalvia_simI.csv"))
rm(simI)

## simulated genome, II
simII <- sim.DNAseq(size=50000000, GCfreq=0.342)
GC(s2c(simII))
genome_size <- 5000000000 # genome size: 5000Mb
ratio_simII <- genome_size/width(simII)
bivalvia_simII <- recto_digest(simII, recto_REs, lower_size, upper_size, ratio_simII)
bivalvia_simII$class <- "bivalvia"
bivalvia_simII$ref <- "simII"
write.csv(bivalvia_simII, file = here("data/in_silico_results/bivalvia_simII.csv"))
rm(simII)

## combine 
bivalvia <- rbind(bivalvia_pimb, bivalvia_bpla, bivalvia_cgig,
                  bivalvia_simI, bivalvia_simII)
write.csv(bivalvia, file = here("data/in_silico_results/bivalvia.csv"))

## calculate some additional size windows
lower_size <- c(250, 250, 300, 350, 250, 250, 250, 250, 200, 200, 200)
upper_size <- c(350, 400, 400, 450, 340, 330, 320, 300, 250, 260, 280)

Pimbricata <- ref.DNAseq(here("../refgenomes/GCA_002216045.1_PinMar1.0_genomic.fna"),
                         subselect.contigs = T, prop.contigs = 0.5)
bivalvia_pimb <- recto_digest(Pimbricata, recto_REs, lower_size, upper_size, 2)
bivalvia_pimb$class <- "bivalvia"
bivalvia_pimb$ref <- "Pinctada imbricata"
rm(Pimbricata)
Bplatifrons <- ref.DNAseq(here("../refgenomes/GCA_002080005.1_Bpl_v1.0_genomic.fna"),
                          subselect.contigs = T, prop.contigs = 0.25)
bivalvia_bpla <- recto_digest(Bplatifrons, recto_REs, lower_size, upper_size, 4)
bivalvia_bpla$class <- "bivalvia"
bivalvia_bpla$ref <- "Bathymodiolus platifrons"
rm(Bplatifrons)
Cgigas <- ref.DNAseq(here("../refgenomes/GCF_000297895.1_oyster_v9_genomic.fna"),
                     subselect.contigs = F)
bivalvia_cgig <- recto_digest(Cgigas, recto_REs, lower_size, upper_size, 1)
bivalvia_cgig$class <- "bivalvia"
bivalvia_cgig$ref <- "Crassostrea gigas"
rm(Cgigas)
simI <- sim.DNAseq(size=10000000, GCfreq=0.353)
genome_size <- 1000000000 # genome size: 1000Mb
ratio_simI <- genome_size/width(simI)
bivalvia_simI <- recto_digest(simI, recto_REs, lower_size, upper_size, ratio_simI)
bivalvia_simI$class <- "bivalvia"
bivalvia_simI$ref <- "simI"
rm(simI)
simII <- sim.DNAseq(size=50000000, GCfreq=0.342)
genome_size <- 5000000000 # genome size: 5000Mb
ratio_simII <- genome_size/width(simII)
bivalvia_simII <- recto_digest(simII, recto_REs, lower_size, upper_size, ratio_simII)
bivalvia_simII$class <- "bivalvia"
bivalvia_simII$ref <- "simII"
rm(simII)
bivalvia2 <- rbind(bivalvia_pimb, bivalvia_bpla, bivalvia_cgig,
                  bivalvia_simI, bivalvia_simII)
write.csv(bivalvia2, file = here("data/in_silico_results/bivalvia.csv"))

## and even more size windows
lower_size <- c(200, 200, 200, 200, 200, 200, 250, 250, 250, 200, 300)
upper_size <- c(300, 350, 400, 450, 500, 320, 420, 450, 500, 500, 500)

Pimbricata <- ref.DNAseq(here("../refgenomes/GCA_002216045.1_PinMar1.0_genomic.fna"),
                         subselect.contigs = T, prop.contigs = 0.5)
bivalvia_pimb <- recto_digest(Pimbricata, recto_REs, lower_size, upper_size, 1)
bivalvia_pimb$class <- "bivalvia"
bivalvia_pimb$ref <- "Pinctada imbricata"
rm(Pimbricata)
Bplatifrons <- ref.DNAseq(here("../refgenomes/GCA_002080005.1_Bpl_v1.0_genomic.fna"),
                          subselect.contigs = T, prop.contigs = 0.25)
bivalvia_bpla <- recto_digest(Bplatifrons, recto_REs, lower_size, upper_size, 4)
bivalvia_bpla$class <- "bivalvia"
bivalvia_bpla$ref <- "Bathymodiolus platifrons"
rm(Bplatifrons)
Cgigas <- ref.DNAseq(here("../refgenomes/GCF_000297895.1_oyster_v9_genomic.fna"),
                     subselect.contigs = F)
bivalvia_cgig <- recto_digest(Cgigas, recto_REs, lower_size, upper_size, 1)
bivalvia_cgig$class <- "bivalvia"
bivalvia_cgig$ref <- "Crassostrea gigas"
rm(Cgigas)
simI <- sim.DNAseq(size=10000000, GCfreq=0.353)
genome_size <- 1000000000 # genome size: 1000Mb
ratio_simI <- genome_size/width(simI)
bivalvia_simI <- recto_digest(simI, recto_REs, lower_size, upper_size, ratio_simI)
bivalvia_simI$class <- "bivalvia"
bivalvia_simI$ref <- "simI"
rm(simI)
simII <- sim.DNAseq(size=50000000, GCfreq=0.342)
genome_size <- 5000000000 # genome size: 5000Mb
ratio_simII <- genome_size/width(simII)
bivalvia_simII <- recto_digest(simII, recto_REs, lower_size, upper_size, ratio_simII)
bivalvia_simII$class <- "bivalvia"
bivalvia_simII$ref <- "simII"
rm(simII)
bivalvia3 <- rbind(bivalvia_pimb, bivalvia_bpla, bivalvia_cgig,
                   bivalvia_simI, bivalvia_simII)
write.csv(bivalvia3, file = here("data/in_silico_results/bivalvia.csv"))

## clean up
rm(genome_size, ratio_simI, ratio_simII, bivalvia,
   bivalvia_pimb, bivalvia_bpla, bivalvia_cgig,
   bivalvia_simI, bivalvia_simII, bivalvia2, bivalvia3)
#####

#### Asteroidea
#####
## reference genome A. planci
#### the reference genomes are too big to be hosted on github
#### you need to get a local copy and store it somewhere and change the file paths accordingly
Aplanci <- ref.DNAseq(here("../refgenomes/GCF_001949145.1_OKI-Apl_1.0_genomic.fna"),
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
#### the reference genomes are too big to be hosted on github
#### you need to get a local copy and store it somewhere and change the file paths accordingly
Pminiata <- ref.DNAseq(here("../refgenomes/GCA_000285935.1_Pmin_1.0_genomic.fna"),
                       subselect.contigs = F)
width(Pminiata)
# 811.0
#GC(s2c(Pminiata))
# 0.402 GC
# digest
asteroidea_pmin <- recto_digest(Pminiata, recto_REs, lower_size, upper_size, 1)
asteroidea_pmin$class <- "asteroidea"
asteroidea_pmin$ref <- "Patiria miniata"
write.csv(asteroidea_pmin, file = here("data/in_silico_results/asteroidea_pminiata.csv"))
rm(Pminiata)

## reference genome P. regularis
#### the reference genomes are too big to be hosted on github
#### you need to get a local copy and store it somewhere and change the file paths accordingly
Pregularis <- ref.DNAseq(here("../refgenomes/GCA_900067625.1_Patiriella_regularis_genome_assembly_1.0_genomic.fna"),
                         subselect.contigs = T, prop.contigs = 0.5)
width(Pregularis)
# 949.3 Mb (474 724 778 * 2)
#GC(s2c(Pregularis))
# 0.404 GC 
# digest
asteroidea_preg <- recto_digest(Pregularis, recto_REs, lower_size, upper_size, 2)
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

## calculate some additional size windows
lower_size <- c(250, 250, 300, 350, 250, 250, 250, 250, 200, 200, 200)
upper_size <- c(350, 400, 400, 450, 340, 330, 320, 300, 250, 260, 280)

Aplanci <- ref.DNAseq(here("../refgenomes/GCF_001949145.1_OKI-Apl_1.0_genomic.fna"),
                      subselect.contigs = F)
asteroidea_apla <- recto_digest(Aplanci, recto_REs, lower_size, upper_size, 1)
asteroidea_apla$class <- "asteroidea"
asteroidea_apla$ref <- "Acanthaster planci"
rm(Aplanci)
Pminiata <- ref.DNAseq(here("../refgenomes/GCA_000285935.1_Pmin_1.0_genomic.fna"),
                       subselect.contigs = T, prop.contigs = 0.5)
asteroidea_pmin <- recto_digest(Pminiata, recto_REs, lower_size, upper_size, 2)
asteroidea_pmin$class <- "asteroidea"
asteroidea_pmin$ref <- "Patiria miniata"
rm(Pminiata)
Pregularis <- ref.DNAseq(here("../refgenomes/GCA_900067625.1_Patiriella_regularis_genome_assembly_1.0_genomic.fna"),
                         subselect.contigs = T, prop.contigs = 0.5)
asteroidea_preg <- recto_digest(Pregularis, recto_REs, lower_size, upper_size, 2)
asteroidea_preg$class <- "asteroidea"
asteroidea_preg$ref <- "Patririella regularis"
rm(Pregularis)
simI <- sim.DNAseq(size=100000000, GCfreq=0.413)
genome_size <- 1000000000 # genome size: 1000Mb
ratio_simI <- genome_size/width(simI)
asteroidea_simI <- recto_digest(simI, recto_REs, lower_size, upper_size, ratio_simI)
asteroidea_simI$class <- "asteroidea"
asteroidea_simI$ref <- "simI"
rm(simI)
simII <- sim.DNAseq(size=200000000, GCfreq=0.404)
genome_size <- 2000000000 # genome size: 2000Mb
ratio_simII <- genome_size/width(simII)
asteroidea_simII <- recto_digest(simII, recto_REs, lower_size, upper_size, ratio_simII)
asteroidea_simII$class <- "asteroidea"
asteroidea_simII$ref <- "simII"
rm(simII)
asteroidea2 <- rbind(asteroidea_apla, asteroidea_pmin, asteroidea_preg,
                    asteroidea_simI, asteroidea_simII)
write.csv(asteroidea2, file = here("data/in_silico_results/asteroidea2.csv"))

## and even more size windows
lower_size <- c(200, 200, 200, 200, 200, 200, 250, 250, 250, 200, 300)
upper_size <- c(300, 350, 400, 450, 500, 320, 420, 450, 500, 500, 500)

Aplanci <- ref.DNAseq(here("../refgenomes/GCF_001949145.1_OKI-Apl_1.0_genomic.fna"),
                      subselect.contigs = F)
asteroidea_apla <- recto_digest(Aplanci, recto_REs, lower_size, upper_size, 1)
asteroidea_apla$class <- "asteroidea"
asteroidea_apla$ref <- "Acanthaster planci"
rm(Aplanci)
Pminiata <- ref.DNAseq(here("../refgenomes/GCA_000285935.1_Pmin_1.0_genomic.fna"),
                       subselect.contigs = F)
asteroidea_pmin <- recto_digest(Pminiata, recto_REs, lower_size, upper_size, 1)
asteroidea_pmin$class <- "asteroidea"
asteroidea_pmin$ref <- "Patiria miniata"
rm(Pminiata)
Pregularis <- ref.DNAseq(here("../refgenomes/GCA_900067625.1_Patiriella_regularis_genome_assembly_1.0_genomic.fna"),
                         subselect.contigs = T, prop.contigs = 0.5)
asteroidea_preg <- recto_digest(Pregularis, recto_REs, lower_size, upper_size, 2)
asteroidea_preg$class <- "asteroidea"
asteroidea_preg$ref <- "Patririella regularis"
rm(Pregularis)
simI <- sim.DNAseq(size=100000000, GCfreq=0.413)
genome_size <- 1000000000 # genome size: 1000Mb
ratio_simI <- genome_size/width(simI)
asteroidea_simI <- recto_digest(simI, recto_REs, lower_size, upper_size, ratio_simI)
asteroidea_simI$class <- "asteroidea"
asteroidea_simI$ref <- "simI"
rm(simI)
simII <- sim.DNAseq(size=200000000, GCfreq=0.404)
genome_size <- 2000000000 # genome size: 2000Mb
ratio_simII <- genome_size/width(simII)
asteroidea_simII <- recto_digest(simII, recto_REs, lower_size, upper_size, ratio_simII)
asteroidea_simII$class <- "asteroidea"
asteroidea_simII$ref <- "simII"
rm(simII)
asteroidea3 <- rbind(asteroidea_apla, asteroidea_pmin, asteroidea_preg,
                     asteroidea_simI, asteroidea_simII)
write.csv(asteroidea3, file = here("data/in_silico_results/asteroidea3.csv"))

## clean up
rm(genome_size, ratio_simI, ratio_simII, asteroidea,
   asteroidea_apla, asteroidea_preg, asteroidea_pmin,
   asteroidea_simI, asteroidea_simII, asteroidea2, asteroidea3)
#####

#### Actinopterygii
#####
## reference genome N. coriiceps
#### the reference genomes are too big to be hosted on github
#### you need to get a local copy and store it somewhere and change the file paths accordingly
Ncoriiceps <- ref.DNAseq(here("../refgenomes/GCF_000735185.1_NC01_genomic.fna"),
                                   subselect.contigs = F)
width(Ncoriiceps)
# 636.6  Mb
GC(s2c(Ncoriiceps))
# 0.408 GC
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


## calculate some additional size windows
lower_size <- c(250, 250, 300, 350, 250, 250, 250, 250, 200, 200, 200)
upper_size <- c(350, 400, 400, 450, 340, 330, 320, 300, 250, 260, 280)

Ncoriiceps <- ref.DNAseq(here("../refgenomes/GCF_000735185.1_NC01_genomic.fna"),
                         subselect.contigs = F)
actinopterygii_ncor <- recto_digest(Ncoriiceps, recto_REs, lower_size, upper_size, 1)
actinopterygii_ncor$class <- "actinopterygii"
actinopterygii_ncor$ref <- "Notothenia coriiceps"
rm(Ncoriiceps)
simI <- sim.DNAseq(size=100000000, GCfreq=0.408)
genome_size <- 1000000000 # genome size: 1000Mb
ratio_simI <- genome_size/width(simI)
actinopterygii_simI <- recto_digest(simI, recto_REs, lower_size, upper_size, ratio_simI)
actinopterygii_simI$class <- "actinopterygii"
actinopterygii_simI$ref <- "simI"
rm(simI)
simII <- sim.DNAseq(size=200000000, GCfreq=0.408)
genome_size <- 1800000000 # genome size: 1800Mb
ratio_simII <- genome_size/width(simII)
actinopterygii_simII <- recto_digest(simII, recto_REs, lower_size, upper_size, ratio_simII)
actinopterygii_simII$class <- "actinopterygii"
actinopterygii_simII$ref <- "simII"
rm(simII)
actinopterygii2 <- rbind(actinopterygii_ncor, actinopterygii_simI, actinopterygii_simII)
write.csv(actinopterygii2, file = here("data/in_silico_results/actinopterygii2.csv"))

## and even more size windows
lower_size <- c(200, 200, 200, 200, 200, 200, 250, 250, 250, 200, 300)
upper_size <- c(300, 350, 400, 450, 500, 320, 420, 450, 500, 500, 500)

Ncoriiceps <- ref.DNAseq(here("../refgenomes/GCF_000735185.1_NC01_genomic.fna"),
                         subselect.contigs = F)
actinopterygii_ncor <- recto_digest(Ncoriiceps, recto_REs, lower_size, upper_size, 1)
actinopterygii_ncor$class <- "actinopterygii"
actinopterygii_ncor$ref <- "Notothenia coriiceps"
rm(Ncoriiceps)
simI <- sim.DNAseq(size=100000000, GCfreq=0.408)
genome_size <- 1000000000 # genome size: 1000Mb
ratio_simI <- genome_size/width(simI)
actinopterygii_simI <- recto_digest(simI, recto_REs, lower_size, upper_size, ratio_simI)
actinopterygii_simI$class <- "actinopterygii"
actinopterygii_simI$ref <- "simI"
rm(simI)
simII <- sim.DNAseq(size=200000000, GCfreq=0.408)
genome_size <- 1800000000 # genome size: 1800Mb
ratio_simII <- genome_size/width(simII)
actinopterygii_simII <- recto_digest(simII, recto_REs, lower_size, upper_size, ratio_simII)
actinopterygii_simII$class <- "actinopterygii"
actinopterygii_simII$ref <- "simII"
rm(simII)
actinopterygii3 <- rbind(actinopterygii_ncor, actinopterygii_simI, actinopterygii_simII)
write.csv(actinopterygii3, file = here("data/in_silico_results/actinopterygii3.csv"))

## clean up
rm(genome_size, ratio_simI, ratio_simII, actinopterygii, actinopterygii2, actinopterygii3,
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

