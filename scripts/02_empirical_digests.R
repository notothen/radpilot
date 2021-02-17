#### Script for reading in restriction enzyme
#### digestion data from bioanalyzer output files
## 17/02/2021
## H. Christiansen
## v2.3

#### NOT RUN
#### install package
#####
## latest version throws an error when trying to load the xml files
## USE THIS VERSION (0.5.1):
#install.packages("https://github.com/jwfoley/bioanalyzeR/releases/download/v0.5.1/bioanalyzeR_0.5.1-no_data.tar.gz", repos = NULL)
#####

#### load packages
library(here) # to shorten file paths
library(bioanalyzeR) # to read bioanalyzer files
library(tidyverse) # to arrange data
library(gridExtra) # to combine plots
library(SimRAD) # for in silico digestion
source(here("scripts/recto_REs_and_functions.R")) # custom functions

#### read in data from bioanalyzer
#### and re-map to the actual sample names
#####
## read in metadata from the runs
samples <- read.csv(here("data/bioanalyzer_results/run_overview.csv"))
samples <- as_tibble(samples)

## function to extract sample names per run
run_sample_names <-  function(samples, run_nr) {
  run_names <- samples %>%
    filter(run == run_nr) %>% 
    select(species, enzyme, replicate) %>%
    transmute(sample_names = str_c(species, enzyme, replicate, sep = "_"))
  return(deframe(run_names))
}

## read in bioanalyzer data for each run
run1 <- read.bioanalyzer(here("data/bioanalyzer_results/run01_2017-12-01_09-25-08.xml"))
run2 <- read.bioanalyzer(here("data/bioanalyzer_results/run02_2017-12-01_10-21-28.xml"))
run3 <- read.bioanalyzer(here("data/bioanalyzer_results/run03_2018-06-15_11-05-44.xml"))
run4 <- read.bioanalyzer(here("data/bioanalyzer_results/run04_2018-06-20_09-30-57.xml"))
run5 <- read.bioanalyzer(here("data/bioanalyzer_results/run05_2018-07-05_10-45-54.xml"))
run6 <- read.bioanalyzer(here("data/bioanalyzer_results/run06_2018-08-03_09-07-54.xml"))
run7 <- read.bioanalyzer(here("data/bioanalyzer_results/run07_2018-08-03_10-18-03_NEW.xml"))
run8 <- read.bioanalyzer(here("data/bioanalyzer_results/run08_2018-09-07_09-46-06.xml"))
run9 <- read.bioanalyzer(here("data/bioanalyzer_results/run09_2019-07-22_11-23-48.xml"))
run10 <- read.bioanalyzer(here("data/bioanalyzer_results/run10_2019-07-25_11-54-39.xml"))

## rename samples in bioanalyzer data
library(plyr)
default_names <- c("sample 1", "sample 2", "sample 3", "sample 4", "sample 5", "sample 6", "sample 7", "sample 8", "sample 9", "sample 10", "sample 11")
run1$samples$sample.name <- mapvalues(run1$samples$sample.name, from = default_names, to = run_sample_names(samples, 1))
run2$samples$sample.name <- mapvalues(run2$samples$sample.name, from = default_names, to = run_sample_names(samples, 2))
run3$samples$sample.name <- mapvalues(run3$samples$sample.name, from = default_names, to = run_sample_names(samples, 3))
run4$samples$sample.name <- mapvalues(run4$samples$sample.name, from = default_names, to = run_sample_names(samples, 4))
run5$samples$sample.name <- mapvalues(run5$samples$sample.name, from = default_names, to = run_sample_names(samples, 5))
run6$samples$sample.name <- mapvalues(run6$samples$sample.name, from = default_names, to = run_sample_names(samples, 6))
run7$samples$sample.name <- mapvalues(run7$samples$sample.name, from = default_names, to = run_sample_names(samples, 7))
run8$samples$sample.name <- mapvalues(run8$samples$sample.name, from = default_names, to = run_sample_names(samples, 8))
run9$samples$sample.name <- mapvalues(run9$samples$sample.name, from = default_names, to = run_sample_names(samples, 9))
run10$samples$sample.name <- mapvalues(run10$samples$sample.name, from = default_names, to = run_sample_names(samples, 10))
detach(package:plyr) # conflicts with package:here
rm(default_names)
#####

#### helper functions
#####
# function to identify which runs were used per species
runs_vs_species <-  function(samples, species_sel) {
  runs <- samples %>%
    filter(species %in% species_sel) %>% 
    select(run, species, enzyme, replicate)
  return(runs)
}

# function to grab samples names per species for subsetting
species_sample_names <-  function(samples, species_sel) {
  species_names <- samples %>%
    filter(species %in% species_sel) %>% 
    select(species, enzyme, replicate) %>%
    transmute(sample_names = str_c(species, enzyme, replicate, sep = "_"))
  return(deframe(species_names))
}
#####

#### ostracoda
#####
## subset
unique(samples$species)
runs_vs_species(samples, c("Macropyxis hornei", "Macrocyprina rocas", "Macroscapha falcis"))
ostra_msp1a <- subset(run2, sample.name %in% species_sample_names(samples, c("Macropyxis hornei", "Macrocyprina rocas", "Macroscapha falcis"))[1:3])
ostra_msp1b <- subset(run10, sample.name %in% species_sample_names(samples, c("Macropyxis hornei", "Macrocyprina rocas", "Macroscapha falcis"))[12:14])
ostra_ecor1a <- subset(run6, sample.name %in% species_sample_names(samples, c("Macropyxis hornei", "Macrocyprina rocas", "Macroscapha falcis"))[4:6])
ostra_ecor1b <- subset(run9, sample.name %in% species_sample_names(samples, c("Macropyxis hornei", "Macrocyprina rocas", "Macroscapha falcis"))[7:9])
ostra_ecor1c <- subset(run9, sample.name %in% species_sample_names(samples, c("Macropyxis hornei", "Macrocyprina rocas", "Macroscapha falcis"))[9:11])
ostra_ecor1d <- subset(run10, sample.name %in% species_sample_names(samples, c("Macropyxis hornei", "Macrocyprina rocas", "Macroscapha falcis"))[15])

## load reference genomes
Cyprideis_torosa <- ref.DNAseq(here("data/refgenomes/Cyprideis_torosa.raw_filtered.fasta"),
                      subselect.contigs = F)
#Cyprideis_torosa <- ref.DNAseq("/media/henrik/Science_HC/01_Science/R_InSilico/refgenomes/Cyprideis_torosa.raw_filtered.fasta", subselect.contigs = F) # reference genome saved on external drive
#Cyprideis_torosa <- ref.DNAseq("E:/01_Science/R_InSilico/refgenomes/Cyprideis_torosa.raw_filtered.fasta", subselect.contigs = F) # reference genome saved on external drive
sim_100mb <- sim.DNAseq(size=100000000, GCfreq=0.439)
sim_500mb <- sim.DNAseq(size=500000000, GCfreq=0.439)
ref_genomes <- data.frame(Cyprideis_torosa, sim_100mb, sim_500mb)

## plot
p1 <- recto_qplot_1(ostra_msp1a, 0, 50, 0.7, 0.7)
seqwidth_msp1 <- recto_digest_2(ref_genomes, msp1)
p2 <- recto_digest_ggplot(seqwidth_msp1, 100000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("ostracoda_msp1a", p3)

p1 <- recto_qplot_1(ostra_msp1b, 0, 5, 0.7, 0.7)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("ostracoda_msp1b", p3)

p1 <- recto_qplot_1(ostra_ecor1a, 0, 50, 0.7, 0.7)
seqwidth_ecor1 <- recto_digest_2(ref_genomes, ecor1)
p2 <- recto_digest_ggplot(seqwidth_ecor1, 30000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("ostracoda_ecor1a", p3)

p1 <- recto_qplot_1(ostra_ecor1b, 0, 20, 0.7, 0.7)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("ostracoda_ecor1b", p3)

p1 <- recto_qplot_1(ostra_ecor1c, 0, 20, 0.7, 0.7)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("ostracoda_ecor1c", p3)

p1 <- recto_qplot_1(ostra_ecor1d, 0, 2, 0.7, 0.7)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("ostracoda_ecor1d", p3)

## clean up
rm(ostra_msp1a, ostra_msp1b, ostra_ecor1a, ostra_ecor1b, ostra_ecor1c, ostra_ecor1d, p1, p2, p3,
   seqwidth_msp1, seqwidth_ecor1, Cyprideis_torosa, sim_100mb, sim_500mb, ref_genomes)
#####

#### Amphipoda
#####
## subset
unique(samples$species)
runs_vs_species(samples, c("Abyssorchomene gerulicorbis", "Paralicella caperesca", "Eusirus aff perdentatus"))
ager_msp1 <- subset(run8, sample.name %in% species_sample_names(samples, c("Abyssorchomene gerulicorbis", "Paralicella caperesca", "Eusirus aff perdentatus"))[1:3])
ager_ecor1 <- subset(run8, sample.name %in% species_sample_names(samples, c("Abyssorchomene gerulicorbis", "Paralicella caperesca", "Eusirus aff perdentatus"))[7:9])
pcap_msp1 <- subset(run8, sample.name %in% species_sample_names(samples, c("Abyssorchomene gerulicorbis", "Paralicella caperesca", "Eusirus aff perdentatus"))[4:6])
pcap_ecor1 <- subset(run8, sample.name %in% species_sample_names(samples, c("Abyssorchomene gerulicorbis", "Paralicella caperesca", "Eusirus aff perdentatus"))[10:11])
eper_msp1 <- subset(run9, sample.name %in% species_sample_names(samples, c("Abyssorchomene gerulicorbis", "Paralicella caperesca", "Eusirus aff perdentatus"))[12:14])
eper_ecor1a <- subset(run9, sample.name %in% species_sample_names(samples, c("Abyssorchomene gerulicorbis", "Paralicella caperesca", "Eusirus aff perdentatus"))[15:17])
eper_ecor1b <- subset(run10, sample.name %in% species_sample_names(samples, c("Abyssorchomene gerulicorbis", "Paralicella caperesca", "Eusirus aff perdentatus"))[18])

## load reference genomes
Hyalella_azteca <- ref.DNAseq(here("data/refgenomes/GCF_000764305.1_Hazt_2.0_genomic.fna"),
                      subselect.contigs = F)
#Hyalella_azteca <- ref.DNAseq("E:/01_Science/R_InSilico/refgenomes/GCF_000764305.1_Hazt_2.0_genomic.fna", subselect.contigs = F) # reference genome saved on external drive
Parhyale_hawaiensis <- ref.DNAseq(here("refgenomes/GCA_001587735.1_Phaw3.0_genomic.fna"),
                          subselect.contigs = T, prop.contigs = 0.25)
#Parhyale_hawaiensis <- ref.DNAseq("E:/01_Science/R_InSilico/refgenomes/GCA_001587735.1_Phaw3.0_genomic.fna", subselect.contigs = T, prop.contigs = 0.25) # reference genome saved on external drive
Eusirus_perdentatus <- ref.DNAseq(here("refgenomes/Eperdentatus.fasta"),
                           subselect.contigs = F)
Eusirus_perdentatus <- ref.DNAseq("E:/01_Science/R_InSilico/refgenomes/Eperdentatus.fasta", subselect.contigs = F) # reference genome saved on external drive
ref_genomes <- data.frame(Hyalella_azteca, Parhyale_hawaiensis, Eusirus_perdentatus)

## plot
p1 <- recto_qplot_1(ager_msp1, 0, 50, 0.7, 0.7)
seqwidth_msp1 <- recto_digest_2(ref_genomes, msp1)
p2 <- recto_digest_ggplot(seqwidth_msp1, 100000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("ager_msp1", p3)

p1 <- recto_qplot_1(ager_ecor1, 0, 20, 0.7, 0.7)
seqwidth_ecor1 <- recto_digest_2(ref_genomes, ecor1)
p2 <- recto_digest_ggplot(seqwidth_ecor1, 100000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("ager_ecor1", p3)

p1 <- recto_qplot_1(pcap_msp1, 0, 5, 0.7, 0.7)
p2 <- recto_digest_ggplot(seqwidth_msp1, 100000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("pcap_msp1", p3)

p1 <- recto_qplot_1(pcap_ecor1, 0, 3, 0.7, 0.7)
p2 <- recto_digest_ggplot(seqwidth_ecor1, 100000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("pcap_ecor1", p3)

p1 <- recto_qplot_1(eper_msp1, 0, 1, 0.7, 0.7)
p2 <- recto_digest_ggplot(seqwidth_msp1, 100000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("eper_msp1", p3)

p1 <- recto_qplot_1(eper_ecor1a, 0, 0.5, 0.7, 0.7)
p2 <- recto_digest_ggplot(seqwidth_ecor1, 100000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("eper_ecor1a", p3)

p1 <- recto_qplot_1(eper_ecor1b, 0, 1, 0.7, 0.7)
p2 <- recto_digest_ggplot(seqwidth_ecor1, 100000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("eper_ecor1b", p3)

# select the "best" runs from multiple trials
eper_ecor1 <- rbind(subset(eper_ecor1a, well.number %in% c (5, 6)), eper_ecor1b)
eper_ecor1$samples

p1 <- recto_qplot_1(eper_ecor1, 0, 1, 0.7, 0.7)
p2 <- recto_digest_ggplot(seqwidth_ecor1, 100000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("eper_ecor1", p3)

## clean up
rm(ostra_msp1a, ostra_msp1b, ostra_ecor1a, ostra_ecor1b, ostra_ecor1c, ostra_ecor1d, p1, p2, p3,
   seqwidth_msp1, seqwidth_ecor1, Cyprideis_torosa, sim_100mb, sim_500mb, ref_genomes)
#####

#### Bivalvia
#####
## subset
runs_vs_species(samples, "Aequiyoldia eightsii")
aeig_ecor1 <- subset(run1, sample.name %in% species_sample_names(samples, "Aequiyoldia eightsii")[1:3])
aeig_pst1 <- subset(run1, sample.name %in% species_sample_names(samples, "Aequiyoldia eightsii")[4:6])
aeig_msp1 <- subset(run2, sample.name %in% species_sample_names(samples, "Aequiyoldia eightsii")[7:9])

runs_vs_species(samples, "Laternula elliptica")
lell_ecor1 <- subset(run1, sample.name %in% species_sample_names(samples, "Laternula elliptica")[1:3])
lell_pst1 <- rbind(subset(run1, sample.name %in% species_sample_names(samples, "Laternula elliptica")[4:5]),
                   subset(run2, sample.name %in% species_sample_names(samples, "Laternula elliptica")[6]))
lell_msp1 <- subset(run2, sample.name %in% species_sample_names(samples, "Laternula elliptica")[7:9])

## load reference genomes
Pinctada_imbricata <- ref.DNAseq("E:/01_Science/R_InSilico/refgenomes/GCA_002216045.1_PinMar1.0_genomic.fna", subselect.contigs = F) # reference genome saved on external drive
Bathymodiolus_platifrons <- ref.DNAseq("E:/01_Science/R_InSilico/refgenomes/GCA_002080005.1_Bpl_v1.0_genomic.fna", subselect.contigs = F) # reference genome saved on external drive
Crassostrea_gigas <- ref.DNAseq("E:/01_Science/R_InSilico/refgenomes/GCF_000297895.1_oyster_v9_genomic.fna", subselect.contigs = F) # reference genome saved on external drive
ref_genomes <- data.frame(Pinctada_imbricata, Bathymodiolus_platifrons, Crassostrea_gigas)

## plot
p1 <- recto_qplot_1(aeig_ecor1, 0, 4, 0.7, 0.7)
seqwidth_ecor1 <- recto_digest_2(ref_genomes, ecor1)
p2 <- recto_digest_ggplot(seqwidth_ecor1, 30000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("aeig_ecor1", p3)

p1 <- recto_qplot_1(aeig_pst1, 0, 4, 0.7, 0.7)
seqwidth_pst1 <- recto_digest_2(ref_genomes, pst1)
p2 <- recto_digest_ggplot(seqwidth_pst1, 30000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("aeig_pst1", p3)

p1 <- recto_qplot_1(aeig_msp1, 0, 4, 0.7, 0.7)
seqwidth_msp1 <- recto_digest_2(ref_genomes, msp1)
p2 <- recto_digest_ggplot(seqwidth_msp1, 100000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("aeig_msp1", p3)

p1 <- recto_qplot_1(lell_msp1, 0, 10, 0.7, 0.7)
p2 <- recto_digest_ggplot(seqwidth_msp1, 100000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("lell_msp1", p3)

p1 <- recto_qplot_1(lell_pst1, 0, 4, 0.7, 0.7)
p2 <- recto_digest_ggplot(seqwidth_pst1, 30000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("lell_pst1", p3)

p1 <- recto_qplot_1(lell_ecor1, 0, 4, 0.7, 0.7)
p2 <- recto_digest_ggplot(seqwidth_ecor1, 30000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("lell_ecor1", p3)

## clean up
rm(aeig_ecor1, aeig_pst1, aeig_msp1, lell_msp1, lell_pst1, lell_ecor1, p1, p2, p3,
   seqwidth_msp1, seqwidth_ecor1, seqwidth_pst1, Pinctada_imbricata, Bathymodiolus_platifrons, Crassostrea_gigas, ref_genomes)
#####

#### Asteroidea
#####
## subset
runs_vs_species(samples, "Psilaster charcoti")
pcha_ecor1 <- subset(run3, sample.name %in% species_sample_names(samples, "Psilaster charcoti")[1:3])
pcha_pst1 <- rbind(subset(run3, sample.name %in% species_sample_names(samples, "Psilaster charcoti")[4]),
                   subset(run4, sample.name %in% species_sample_names(samples, "Psilaster charcoti")[5:6]))
pcha_msp1 <- subset(run4, sample.name %in% species_sample_names(samples, "Psilaster charcoti")[7:9])

runs_vs_species(samples, "Bathybiaster loripes")
blor_ecor1 <- subset(run3, sample.name %in% species_sample_names(samples, "Bathybiaster loripes")[1:3])
blor_pst1 <- subset(run4, sample.name %in% species_sample_names(samples, "Bathybiaster loripes")[4:6])
blor_msp1 <- rbind(subset(run4, sample.name %in% species_sample_names(samples, "Bathybiaster loripes")[7]),
                   subset(run7, sample.name %in% species_sample_names(samples, "Bathybiaster loripes")[8:9]))


## load reference genomes
Acanthaster_planci <- ref.DNAseq(here("data/refgenomes/GCF_001949145.1_OKI-Apl_1.0_genomic.fna"),
                      subselect.contigs = F)
#Acanthaster_planci <- ref.DNAseq("E:/01_Science/R_InSilico/refgenomes/GCF_001949145.1_OKI-Apl_1.0_genomic.fna", subselect.contigs = F) # reference genome saved on external drive
Patiria_miniata <- ref.DNAseq(here("data/refgenomes/GCA_000285935.1_Pmin_1.0_genomic.fna"),
                       subselect.contigs = F)
#Patiria_miniata <- ref.DNAseq("E:/01_Science/R_InSilico/refgenomes/GCA_000285935.1_Pmin_1.0_genomic.fna", subselect.contigs = F) # reference genome saved on external drive
Patiriella_regularis <- ref.DNAseq(here("data/refgenomes/GCA_900067625.1_Patiriella_regularis_genome_assembly_1.0_genomic.fna"),
                         subselect.contigs = F)
#Patiriella_regularis <- ref.DNAseq("E:/01_Science/R_InSilico/refgenomes/GCA_900067625.1_Patiriella_regularis_genome_assembly_1.0_genomic.fna", subselect.contigs = F) # reference genome saved on external drive
ref_genomes <- data.frame(Acanthaster_planci, Patiria_miniata, Patiriella_regularis)

## plot
p1 <- recto_qplot_1(pcha_ecor1, 0, 2, 0.7, 0.7)
seqwidth_ecor1 <- recto_digest_2(ref_genomes, ecor1)
p2 <- recto_digest_ggplot(seqwidth_ecor1, 30000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("pcha_ecor1", p3)

p1 <- recto_qplot_1(pcha_pst1, 0, 4, 0.7, 0.7)
seqwidth_pst1 <- recto_digest_2(ref_genomes, pst1)
p2 <- recto_digest_ggplot(seqwidth_pst1, 30000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("pcha_pst1", p3)

p1 <- recto_qplot_1(pcha_msp1, 0, 4, 0.7, 0.7)
seqwidth_msp1 <- recto_digest_2(ref_genomes, msp1)
p2 <- recto_digest_ggplot(seqwidth_msp1, 100000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("pcha_msp1", p3)

p1 <- recto_qplot_1(blor_msp1, 0, 10, 0.7, 0.7)
p2 <- recto_digest_ggplot(seqwidth_msp1, 100000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("blor_msp1", p3)

p1 <- recto_qplot_1(blor_pst1, 0, 2, 0.7, 0.7)
p2 <- recto_digest_ggplot(seqwidth_pst1, 30000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("blor_pst1", p3)

p1 <- recto_qplot_1(blor_ecor1, 0, 4, 0.7, 0.7)
p2 <- recto_digest_ggplot(seqwidth_ecor1, 30000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("blor_ecor1", p3)

## clean up
rm(pcha_ecor1, pcha_pst1, pcha_msp1, blor_msp1, blor_pst1, blor_ecor1, p1, p2, p3,
   seqwidth_msp1, seqwidth_ecor1, seqwidth_pst1, Acanthaster_planci, Patiria_miniata, Patiriella_regularis, ref_genomes)
#####

#### Actinopterygii
#####
## subset
runs_vs_species(samples, "Trematomus bernacchii")
tber_apek1 <- subset(run5, sample.name %in% species_sample_names(samples, "Trematomus bernacchii")[1:3])
tber_ecor1 <- rbind(subset(run5, sample.name %in% species_sample_names(samples, "Trematomus bernacchii")[4:5]),
                   subset(run7, sample.name %in% species_sample_names(samples, "Trematomus bernacchii")[6:7]))
tber_msp1_ecor1 <- subset(run7, sample.name %in% species_sample_names(samples, "Trematomus bernacchii")[8:10])

runs_vs_species(samples, "Trematomus loennbergii")
tloe_apek1 <- subset(run5, sample.name %in% species_sample_names(samples, "Trematomus loennbergii")[1:3])
tloe_ecor1 <- subset(run5, sample.name %in% species_sample_names(samples, "Trematomus loennbergii")[4:6])
tloe_msp1_ecor1 <- subset(run7, sample.name %in% species_sample_names(samples, "Trematomus loennbergii")[7:9])

## load reference genomes
Notothenia_coriiceps <- ref.DNAseq(here("data/refgenomes/GCF_000735185.1_NC01_genomic.fna"),
                         subselect.contigs = F)
#Notothenia_coriiceps <- ref.DNAseq("E:/01_Science/R_InSilico/refgenomes/GCF_000735185.1_NC01_genomic.fna", subselect.contigs = F) # reference genome saved on external drive
sim_100mb <- sim.DNAseq(size = 100000000, GCfreq = 0.408)
sim_180mb <- sim.DNAseq(size = 180000000, GCfreq = 0.408)
ref_genomes <- data.frame(Notothenia_coriiceps, sim_100mb, sim_180mb)

## plot
p1 <- recto_qplot_1(pcha_ecor1, 0, 2, 0.7, 0.7)
seqwidth_ecor1 <- recto_digest_2(ref_genomes, ecor1)
p2 <- recto_digest_ggplot(seqwidth_ecor1, 30000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("pcha_ecor1", p3)

p1 <- recto_qplot_1(pcha_pst1, 0, 4, 0.7, 0.7)
seqwidth_pst1 <- recto_digest_2(ref_genomes, pst1)
p2 <- recto_digest_ggplot(seqwidth_pst1, 30000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("pcha_pst1", p3)

p1 <- recto_qplot_1(pcha_msp1, 0, 4, 0.7, 0.7)
seqwidth_msp1 <- recto_digest_2(ref_genomes, msp1)
p2 <- recto_digest_ggplot(seqwidth_msp1, 100000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("pcha_msp1", p3)

p1 <- recto_qplot_1(blor_msp1, 0, 10, 0.7, 0.7)
p2 <- recto_digest_ggplot(seqwidth_msp1, 100000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("blor_msp1", p3)

p1 <- recto_qplot_1(blor_pst1, 0, 2, 0.7, 0.7)
p2 <- recto_digest_ggplot(seqwidth_pst1, 30000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("blor_pst1", p3)

p1 <- recto_qplot_1(blor_ecor1, 0, 4, 0.7, 0.7)
p2 <- recto_digest_ggplot(seqwidth_ecor1, 30000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("blor_ecor1", p3)

## clean up
rm(tber_apek1, tber_ecor1, tber_msp1_ecor1, tloe_apek1, tloe_ecor1, tloe_msp1_ecor1, p1, p2, p3,
   seqwidth_msp1, seqwidth_ecor1, seqwidth_pst1, Notothenia_coriiceps, ref_genomes)
#####

#### Aves
#####
## subset
runs_vs_species(samples, "Pagodroma nivea")
pniv_ecor1a <- rbind(subset(run3, sample.name %in% species_sample_names(samples, "Pagodroma nivea")[1:2]),
                    subset(run6, sample.name %in% species_sample_names(samples, "Pagodroma nivea")[8]))
pniv_ecor1b <- rbind(subset(run10, sample.name %in% species_sample_names(samples, "Pagodroma nivea")[12]),
                    subset(run10, sample.name %in% species_sample_names(samples, "Pagodroma nivea")[15]))
pniv_pst1 <- rbind(subset(run3, sample.name %in% species_sample_names(samples, "Pagodroma nivea")[3:4]),
                   subset(run7, sample.name %in% species_sample_names(samples, "Pagodroma nivea")[9]),
                   subset(run10, sample.name %in% species_sample_names(samples, "Pagodroma nivea")[11]),
                   subset(run10, sample.name %in% species_sample_names(samples, "Pagodroma nivea")[14]))
pniv_msp1 <- rbind(subset(run4, sample.name %in% species_sample_names(samples, "Pagodroma nivea")[5:6]),
                   subset(run6, sample.name %in% species_sample_names(samples, "Pagodroma nivea")[7]),
                   subset(run10, sample.name %in% species_sample_names(samples, "Pagodroma nivea")[10]),
                   subset(run10, sample.name %in% species_sample_names(samples, "Pagodroma nivea")[13]))


## load reference genomes
Fulmarus_glacialis <- ref.DNAseq(here("data/refgenomes/GCF_000690835.1_ASM69083v1_genomic.fna"),
                         subselect.contigs = F)
#Fulmarus_glacialis <- ref.DNAseq("E:/01_Science/R_InSilico/refgenomes/GCF_000690835.1_ASM69083v1_genomic.fna", subselect.contigs = F) # reference genome saved on external drive
sim_100mb <- sim.DNAseq(size = 100000000, GCfreq = 0.412)
sim_150mb <- sim.DNAseq(size = 150000000, GCfreq = 0.412)
ref_genomes <- data.frame(Fulmarus_glacialis, sim_100mb, sim_150mb)

## plot
pniv_ecor1a$samples
recto_qplot_1(pniv_ecor1a, 0, 2, 0.7, 0.7)
pniv_ecor1b$samples
recto_qplot_1(pniv_ecor1b, 0, 2, 0.7, 0.7)
pniv_ecor1 <- rbind(subset(run3, sample.name %in% species_sample_names(samples, "Pagodroma nivea")[1:2]),
                    subset(run10, sample.name == species_sample_names(samples, "Pagodroma nivea")[12]))
pniv_ecor1 <- subset(pniv_ecor1, well.number %in% c(1, 8, 11))
pniv_ecor1$samples
recto_qplot_1(pniv_ecor1, 0, 2, 0.7, 0.7)



p1 <- recto_qplot_1(pniv_ecor1, 0, 2, 0.7, 0.7)
p1
seqwidth_ecor1 <- recto_digest_2(ref_genomes, ecor1)
p2 <- recto_digest_ggplot(seqwidth_ecor1, 30000)
p1
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("pniva_ecor1", p3)

p1 <- recto_qplot_1(pcha_pst1, 0, 4, 0.7, 0.7)
seqwidth_pst1 <- recto_digest_2(ref_genomes, pst1)
p2 <- recto_digest_ggplot(seqwidth_pst1, 30000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("pcha_pst1", p3)

p1 <- recto_qplot_1(pcha_msp1, 0, 4, 0.7, 0.7)
seqwidth_msp1 <- recto_digest_2(ref_genomes, msp1)
p2 <- recto_digest_ggplot(seqwidth_msp1, 100000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("pcha_msp1", p3)

p1 <- recto_qplot_1(blor_msp1, 0, 10, 0.7, 0.7)
p2 <- recto_digest_ggplot(seqwidth_msp1, 100000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("blor_msp1", p3)

p1 <- recto_qplot_1(blor_pst1, 0, 2, 0.7, 0.7)
p2 <- recto_digest_ggplot(seqwidth_pst1, 30000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("blor_pst1", p3)

p1 <- recto_qplot_1(blor_ecor1, 0, 4, 0.7, 0.7)
p2 <- recto_digest_ggplot(seqwidth_ecor1, 30000)
p3 <- arrangeGrob(p1, p2, nrow = 1)
recto_ggsave("blor_ecor1", p3)

## clean up
rm(tber_apek1, tber_ecor1, tber_msp1_ecor1, tloe_apek1, tloe_ecor1, tloe_msp1_ecor1, p1, p2, p3,
   seqwidth_msp1, seqwidth_ecor1, seqwidth_pst1, Notothenia_coriiceps, ref_genomes)
#####