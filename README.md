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

**If you use any of this code for your own work, don't forget to cite the R package providers (see LOG for a list of the most important packages)!**

## Applying this code for your own setup

In principle you could just use the existing tools from the R packages **SimRAD** (Lepais & Weir 2014) and **bioanalyzeR** (Foley 2021).
If you want to apply similar calculations as in Christiansen et al. 2021, then start by loading the source script ```recto_REs_and_functions.R``` from the ```scripts``` folder: https://github.com/notothen/radpilot/tree/main/scripts

This script contains R functions used for the purposes as in Christiansen et al. 2021. Load it like so:

```
library(here)
source(here("scripts/recto_REs_and_functions.R"))
```

If you don't have the **here** package yet, install it first: ```install.packages("here")```. It's a very useful package that helps making code function on different working platforms without the need to explicitly define a working directory.

Then either clone the github directory, or download the ```recto_REs_and_functions.R``` script and put it in a folder called ```scripts```.

# Part 1: In silico genome digestions

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

In order to conduct in silico digestions with these two genomes using the restrictions enzymes and size windows as in Christiansen et al. 2021 all you need to do is run:

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


and in these size windows:

```
lower size: 210, 240, 0, 100, 200, 300, 400, 500, 600, 700, 800
upper size: 260, 340, 100, 200, 300, 400, 500, 600, 700, 800, 900
```

You can adjust the size windows by providing any other sequence of numbers as ```lower_size``` and ```upper_size```, e.g.:

```
lower_size <- c(210, 240, 0, 100, 200, 300, 400, 500, 600, 700, 800)
upper_size <- c(260, 340, 100, 200, 300, 400, 500, 600, 700, 800, 900)
```

You can find an extended application of these in silico calculations with the target taxa from Christiansen et al. 2021 in the script ```01_digests.R``` from the the ```scripts``` folder: https://github.com/notothen/radpilot/tree/main/scripts

# Part 2:


