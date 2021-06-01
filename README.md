# RRS pilot experiment calculations

#### By Henrik Christiansen, last update on 01/06/2021
#### Reference: Christiansen et al. 2021 **"Facilitating population genomics of non-model organisms through optimized experimental design for reduced representation sequencing"** [10.1101/2021.03.30.437642](https://doi.org/10.1101/2021.03.30.437642) 
#### All code and input files (except for the reference genomes) available on github: [github.com/notothen/radpilot](https://github.com/notothen/radpilot) 
#### Sequence data from test libraries available on SRA ([reviewer link](https://dataview.ncbi.nlm.nih.gov/object/PRJNA674352?reviewer=dmj6c5816761lpqn2d1oe69qhs)) 

## Introduction

The material presented here contains R code, input data in spreadsheet files (tsv/csv), and output plots in different formats (eps/jpg/pdf/tiff) used for a
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
simgenome1 <- sim.DNAseq(size=500000000, GCfreq=0.5)
GC(s2c(simgenome1))

genome_size <- 1000000000 # genome size: 1000 Mb
ratio_simgenome1 <- genome_size/width(simgenome1)
ratio_simgenome1
```

We have now simulated a 1,000 Mb genome with 50 % GC content. In fact, we only simulated a 500 Mb genome to make it less computationally intensive. The results will be extrapolated to a 1,000 Mb genome using the ```ratio_simgenome1``` value. You can adjust this in different ways, depending on the computational power of your system.

In order to conduct *in silico* digestions with these two genomes using the restrictions enzymes and size windows as in Christiansen et al. 2021 all you need to do is run:

```
refgenome1_digest <- recto_digest(refgenome1, recto_REs, lower_size, upper_size, 1)
simgenome1_digest <- recto_digest(simgenome1, recto_REs, lower_size, upper_size, ratio_simgenome1)
```

You can also create a folder ```data\in_silico_results``` and save your results there:

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

and in these size windows:

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

You can find an extended application of these *in silico* calculations with the target taxa from Christiansen et al. 2021 in the script ```01_digests.R``` from the the ```scripts``` folder: https://github.com/notothen/radpilot/tree/main/scripts

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

If you want to plot the concentration over fragment size in a similar fashion as presented in Christiansen et al. 2021, then you could do so like this:

```
recto_qplot_2(run1, 0, 50)
```

The two values supplied to the ```recto_qplot_2``` function are for the y axis range limits. You may need to experiment with these. You may also want to rename and subset your Bioanalyzer run data, but in principle that's all!

More details including subsetting and renaming and an extended application where bioanalyzer data is imported in R and plotted for the target taxa from Christiansen et al. 2021 can be found in the script ```02_empirical_digests.R``` from the the ```scripts``` folder: https://github.com/notothen/radpilot/tree/main/scripts

## Part 3: Estimate marker density

TO BE WRITTEN

## Part 4: Evaluate test library output

If you have run some test sequencing libraries and processed them following Rochette & Catchen (2017), then you may want to plot the number of loci and the number of SNPs per locus. To do this you can simply follow the codes provided by these authors. Their codes to run the Stacks pipeline also provide tsv outputs on the number of loci and SNPs. Put this file into a folder ```data/test_libraries``` and load it e.g. like this:

```
snps_per_loc <- read.delim(here('data/test_libraries/n_snps_per_locus.tsv'))
```

You can then just proceed to follow the R code of Rochette & Catchen (2017). You can find it here:

https://bitbucket.org/rochette/rad-seq-genotyping-demo/src/master/demo_scripts/R_scripts/

A detailed adaptation of this code using the data of five test sequencing libraries as in Christiansen et al. 2021 can be found in the scripts ```03_plot_n_loci.R``` and ```04_plot_n_snps_per_locus.R``` from the ```scripts``` folder: https://github.com/notothen/radpilot/tree/main/scripts

In addition, you may want to plot the realized number of loci and coverage for your test libraries. For this we collated a csv file with all metadata about the sequenced test libraries (the library number, species, specimen ID) and ....CONTINUE HERE


