[![DOI](https://zenodo.org/badge/331342737.svg)](https://zenodo.org/badge/latestdoi/331342737)

# RRS pilot experiment calculations

#### By Henrik Christiansen, last update on 11/10/2021
#### Reference: Christiansen, H., Heindler, F.M., Hellemans, B. et al. Facilitating population genomics of non-model organisms through optimized experimental design for reduced representation sequencing. BMC Genomics 22, 625 (2021). [doi.org/10.1186/s12864-021-07917-3](https://doi.org/10.1186/s12864-021-07917-3)
#### All code and input files (except for the reference genomes) available on github: [github.com/notothen/radpilot](https://github.com/notothen/radpilot) 
#### Sequence data from test libraries available on SRA: [PRJNA674352](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA674352)

## Introduction

The material presented here contains R code, input data in spreadsheet files (tsv/csv), and output plots in different formats (jpg/pdf) used for a
reduced representation sequencing (RRS; including RADseq and GBS) optimization pilot experiment conducted with several Antarctic animal species.

With some adjustments the R code can possibly be used for similar analyses with other target taxa or experimental settings in mind. Below, I provide some
very general comments on how that might be achieved. However, there are no guarantees whatsoever as to the accuracy of the code and that it works in other
situations. Feel free to copy parts of the code and adjust it to your own liking.

An extended log of the project history and details on all input and output files and the R packages used can be found in the [LOG file](../main/LOG.pdf)

**If you use any of this code for your own work, don't forget to cite the R package providers (see [LOG](../main/LOG.pdf) for a list of the most important packages)!**

## Applying this code for your own setup

In principle you could just use the existing tools from the R packages:
* [**SimRAD**](https://cran.r-project.org/web/packages/SimRAD/index.html) (Lepais & Weir 2014, [Molecular Ecology Resources](https://doi.org/10.1111/1755-0998.12273))
* [**bioanalyzeR**](https://stanford.edu/~jwfoley/bioanalyzeR.html) (Foley 2021)

and the code provided here:
* [**Stacks Protocol**](https://bitbucket.org/rochette/rad-seq-genotyping-demo/src/master/demo_scripts/R_scripts/) (Rochette & Catchen 2017, [Nature Protocols](https://www.nature.com/articles/nprot.2017.123))

If you want to apply similar calculations as in Christiansen et al. 2021, then start by loading the source script ```recto_REs_and_functions.R``` from the ```scripts``` folder: https://github.com/notothen/radpilot/tree/main/scripts

This script contains R functions used for the purposes as in Christiansen et al. 2021. Load it like so:

```
library(here)
source(here("scripts/recto_REs_and_functions.R"))
```

If you don't have the **here** package yet, install it first: ```install.packages("here")```. It's a very useful package that helps making code function on different working platforms without the need to explicitly define a working directory.

Then either clone the github directory, or download the ```recto_REs_and_functions.R``` script and put it in a folder called ```scripts```.

## Part 1: *in silico* genome digestions

Start by loading these R packages (or install them first, if you don't have them installed):

```
library(here)
library(SimRAD)
library(seqinr)
library(tidyverse)
```

Now you need a reference genome to computationally digest. If you have one (e.g. from GenBank) you can put it into a folder called ```refgenomes``` and load it as such:

```
refgenome1 <- ref.DNAseq(here("refgenomes/YOUR_REFGENOME.fasta"), subselect.contigs = F)
```

You can also check its size and GC content. **Beware**, this may take some time if the genome is big. If that's the case, you may want to switch to ```subselect.contigs = T``` in the code above. See the **SimRAD** documentation for more details.

```
width(refgenome1)
GC(s2c(refgenome1))
```

Alternatively, you can simulate a genome, like so:

```
simgenome1 <- sim.DNAseq(size = 250000000, GCfreq = 0.5)
GC(s2c(simgenome1))

genome_size <- 1000000000 # genome size: 1000 Mb
ratio_simgenome1 <- genome_size/width(simgenome1)
ratio_simgenome1
```

We have now simulated a 1,000 Mb genome with 50 % GC content. In fact, we only simulated a 250 Mb genome to make it less computationally intensive. The results will be extrapolated to a 1,000 Mb genome using the ```ratio_simgenome1``` value. You can adjust this in different ways, depending on the computational power of your system.

In order to conduct *in silico* digestions with these two genomes using the restrictions enzymes and size windows as in Christiansen et al. 2021 all you need to do is run:

```
refgenome1_digest <- recto_digest(refgenome1, recto_REs, lower_size, upper_size, 1)
simgenome1_digest <- recto_digest(simgenome1, recto_REs, lower_size, upper_size, ratio_simgenome1)
```

You can also create a folder ```data/in_silico_results``` and save your results there:

```
write.csv(refgenome1_digest, file = here("data/in_silico_results/refgenome1_digest.csv"))
write.csv(simgenome1_digest, file = here("data/in_silico_results/simgenome1_digest.csv"))
```

That's it! Now you should have estimates for these enzymes:

|single digest|double digest
|---|---
|*sbf1*|*sbf1*-*sph1*
|*ecor1*|*ecor1*-*msp1*
|*sph1*|*ecor1*-*sph1*
|*pst1*|*pst1*-*msp1*
|*msp1*|*sbf1*-*msp1*
|*apek1*|
|*mse1*|

and in these size windows (plus the total number of fragments):

|lower size |upper size
|---|---
|210|260
|240|340
|0|100
|100|200
|200|300
|300|400
|400|500
|500|600
|600|700
|700|800
|800|900

You can adjust the size windows by providing any other sequence of 11 numbers as ```lower_size``` and ```upper_size```, e.g.:

```
lower_size <- c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500)
upper_size <- c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550)
```

If you want to simulate with other restriction enzymes or digest in more/less than 11 size bins, then you would need to change the respective object (```recto_REs```) or function (```recto_digest```) in the ```recto_REs_and_functions.R``` script or write your own.

You can find an extended application of these *in silico* calculations with the target taxa from Christiansen et al. 2021 in the script ```01_digests.R``` from the ```scripts``` folder: https://github.com/notothen/radpilot/tree/main/scripts

## Part 2: Comparison of *in silico* and empirical digestions

If you have done some empirical enzyme digestion tests and run them on a Bioanalyzer platform, you may want to compare those results with your *in silico* estimates. 

Start by loading these R packages (or install them first, if you don't have them installed):

```
library(here)
library(bioanalyzeR)
library(tidyverse)
library(gridExtra)
library(SimRAD)
```

If you have results from a Bioanalyzer run make sure to export these as XML files. Then you can import them like this e.g.:

```
run1 <- read.bioanalyzer(here("data/bioanalyzer_results/run01_2017-12-01_09-25-08.xml"))
```

The default names in the Bioanalyzer output are not very meaningful. You can replace them, for example like this, with the help of some manually curated metadata and a small helper function to extract the correct names:

```
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

library(plyr)
default_names <- c("sample 1", "sample 2", "sample 3", "sample 4", "sample 5", "sample 6",
                   "sample 7", "sample 8", "sample 9", "sample 10", "sample 11")
run1$samples$sample.name <- mapvalues(run1$samples$sample.name, from = default_names, to = run_sample_names(samples, 1))
detach(package:plyr) # conflicts with package:here
rm(default_names)
```

If you then want to plot the concentration over fragment size in a similar fashion as presented in Christiansen et al. 2021, you could do so like this:

```
recto_qplot_2(run1, 0, 50)
```

The two values supplied to the ```recto_qplot_2``` function are for the y axis range limits. You may need to experiment with these. You may also want to rename and subset your Bioanalyzer run data, but in principle that's all!

More details including subsetting and renaming and an extended application where bioanalyzer data is imported in R and plotted for the target taxa from Christiansen et al. 2021 can be found in the script ```02_empirical_digests.R``` from the ```scripts``` folder: https://github.com/notothen/radpilot/tree/main/scripts

## Part 3 & 4: Evaluate test library output

If you have run some test sequencing libraries and processed them following Rochette & Catchen (2017), then you may want to plot the number of loci and the number of SNPs per locus. To do this you can simply follow the code provided by these authors. Their code to run the Stacks pipeline also provides tsv outputs on the number of loci and SNPs. Put this file into a folder ```data/test_libraries``` and load it e.g. like this:

```
snps_per_loc <- read.delim(here('data/test_libraries/n_snps_per_locus.tsv'))
```

You can then just proceed to follow the R code of [Rochette & Catchen (2017)](https://www.nature.com/articles/nprot.2017.123). You can find it here:

https://bitbucket.org/rochette/rad-seq-genotyping-demo/src/master/demo_scripts/R_scripts/

A detailed adaptation of this code using the data of five test sequencing libraries as in Christiansen et al. 2021 can be found in the scripts ```03_plot_n_loci.R``` and ```04_plot_n_snps_per_locus.R``` from the ```scripts``` folder: https://github.com/notothen/radpilot/tree/main/scripts

## Part 5: Estimate marker density

In order to assess how large a portion of a given genome you are able to scan with a selected RRS setup you may want to assess marker density, i.e. the number of fragments/markers/SNPs that you are able to genotype in relation to the genome size. The function ```marker_density()``` which is included in the source script ```recto_REs_and_functions.R``` can estimate some useful parameters given a genome size, number of fragments, and information about the sequencing setup.

Start by loading these R packages (or install them first, if you don't have them installed) and the source script:

```
library(here)
library(scales)
library(tidyverse)
source(here("scripts/recto_REs_and_functions.R"))
```

Now, let us first define some constant variables. You can change these values, depending on what you plan to sequence on what kind of sequencing output you expect. Here, are the values we used:

```
reads_hiseq4000 <- 300000000 # conservative estimate
reads_hiseq2500 <- 200000000 # conservative estimate
reads_novaseq <- 1800000000 # adjust this depending on what exactly you plan to sequence on!
ind_per_lib <- 96 # adjust this depending on how many individuals you plan to pool per library!
read_length_hiseq4000 <- 150
read_length_hiseq2500 <- 125
read_length_novaseq <- 100 # adjust this depending on what exactly you plan to sequence on!
```

Next, we can use the ```marker_density``` function by supplying:

* the number of fragments you expect (e.g. as simulated in silico) or observed (as found through test sequencing libraries)
* the genome size of the target species (based on available information and/or a best guess)
* the sequencing platform (based on what we defined above either "HiSeq2500", "HiSeq4000", or "NovaSeq")
* whether you plan to do/did paired end sequencing or not
* what fraction of SNPs per base pair you expect (e.g. 1 SNP every 100 bp = 0.01)

A simple execution could look like this:

```
marker_density(fragment = 50000, genome_size = 1500000000, sequencer = "HiSeq2500", paired_end = T, 0.01)
```

The results should be pretty self-explanatory. You could go further and collect more such estimates, collect them in a dataframe (e.g. with ```rbind()```), and save the results (e.g. with ```write.csv()```). The script ```05_marker_density.R``` contains longer code to apply this function on the target taxa from Christiansen et al. 2021. The marker density results reported in Table 4 were estimated in this way. The entire output from that script is available in the file ```marker_density_results.csv``` in the ```data``` folder. The script is available in the ```scripts``` folder:
https://github.com/notothen/radpilot/tree/main/scripts

## Part 6: Plot estimated and realized numbers of loci and coverage

Finally, you may want to plot the realized number of loci and coverage for your test libraries. For this we collated a csv file with metadata and metrics about the sequenced test libraries. The file is called ```coverage_stats.csv``` and located in the ```data/test_libraries``` folder. It contains the library number (lib), the  species or family (species), specimen ID (ind), the obtained coverage when using Stacks parameter M=1 (cov_M_1), the optimal M parameter according to a parameter optimization test (opt_M), the coverage obtained with that M (cov_opt_M), the target coverage as estimated *in silico* (target_cov), the target number of loci as estimated *in silico* (target_loci), the obtained number of loci (n_loci), the number of used forward reads (n_used_fw_reads), the mean coverage (mean_cov), and a weighted mean coverage (mean_cov_ns). The latter values (n_loci, n_used_fw_reads, mean_cov, mean_cov_ns) were retrieved from the output files (gstacks.log.distribs) of the gstacks module of [**Stacks**](https://catchenlab.life.illinois.edu/stacks/). For consistency, these files (gstacks.logs.distribs) are included for each sequencing library in the repository under ```data/test_libraries``` but they are not actually needed to reproduce our results (as everything is collated in the ```coverage_stats.csv``` file. You could create a similar file with your own data.

The script ```06_plot_coverage.R``` contains code to load this metadata file and create a plot as shown in Fig. 2 of Christiansen et al. 2021. The script is also in the ```scripts``` folder:
https://github.com/notothen/radpilot/tree/main/scripts

## Conclusion

Above you have some guidelines to follow the calculations and plotting exercises as they were used for the RRS setup optimization as presented in Christiansen et al. 2021. You can copy all R script files and manipulate them to work with your own data. If you run into trouble, you could always go back to this github repository which contains everything - except for the reference genomes from GenBank - you need to reproduce the exact calculations and plots as in our paper (except for some tiny stochastic variation in the numbers) and try to see where your problem occurs. Alternatively, you could file an [issue](https://github.com/notothen/radpilot/issues) to ask for help or suggest changes.



