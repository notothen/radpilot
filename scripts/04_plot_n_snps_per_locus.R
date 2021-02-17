#### Script to plot the number of SNPs per r80 locus
#### that can be constructed de novo with
#### varying M parameter from test RAD libraries
## 14/12/2020
## v1.6
## H. Christiansen,
## modified from Rochette & Catchen 2017:
browseURL("https://bitbucket.org/rochette/rad-seq-genotyping-demo/src/default/demo_scripts/R_scripts/")

#### load packages
library(here) # to shorten file paths
library(ggsci) # for journal type colors
library(scales) # to visualize ggsci palettes
source("~/printfig.R") # to save plots

#### lib1
#####
# lib1
snps_per_loc <- read.delim(here('data/test_libraries/lib1/n_snps_per_locus_r50.tsv'))
snps_per_loc <- subset(snps_per_loc, M==n & m==3)
snps_per_loc <- snps_per_loc[order(snps_per_loc$n_snps),]
Mn_values <- sort(unique(snps_per_loc$M))
m <- matrix(NA, nrow = length(Mn_values), ncol = max(snps_per_loc$n_snps)+1)
for(i in 1:nrow(snps_per_loc)) {
  m[snps_per_loc$M[i], snps_per_loc$n_snps[i]+1] = snps_per_loc$n_loci[i] # [n_snps+1] as column 1 is for loci with 0 SNPs
}
max_n_snps <- 10
m[, max_n_snps + 2] <- rowSums(m[, (max_n_snps + 2):ncol(m)], na.rm = T)
m <- m[, 1:(max_n_snps + 2)]
lib1 <- m / rowSums(m, na.rm = T)

# call plot
mypal <- pal_npg("nrc", alpha = 0.7)(9)
colfunc <- colorRampPalette(c(mypal[7], 'lightgray'))
colfunc(length(Mn_values))
lib1_plot <- function(x){
  barplot(lib1,
          beside = T, col = colfunc(length(Mn_values)), las = 1,
          names.arg = c(0:max_n_snps, paste('>', max_n_snps, sep = '')),
          xlab = 'number of SNPs',
          ylab = 'percentage of loci'
          )
  legend(100, 0.3, legend = Mn_values, fill = colfunc(length(Mn_values)),
         title = "M and n", bty = 'n')
}
printfig(here('figures/n_snps_lib1'), lib1_plot)
printfigpdf(here('figures/n_snps_lib1'), lib1_plot) 
rm(lib1, m, snps_per_loc, i, max_n_snps, mypal, Mn_values, colfunc, colfunc2, lib1_plot)
#####

#### lib2
#####
# lib2, lat
snps_per_loc <- read.delim(here('data/test_libraries/lib2/laternula/n_snps_per_locus_cat.tsv'))
snps_per_loc <- subset(snps_per_loc, M==n & m==3)
snps_per_loc <- snps_per_loc[order(snps_per_loc$n_snps),]
Mn_values <- sort(unique(snps_per_loc$M))

# Write the counts in a matrix.
m <- matrix(NA, nrow = length(Mn_values), ncol = max(snps_per_loc$n_snps)+1)
for(i in 1:nrow(snps_per_loc)) {
  m[snps_per_loc$M[i], snps_per_loc$n_snps[i]+1] = snps_per_loc$n_loci[i] # [n_snps+1] as column 1 is for loci with 0 SNPs
}

# Truncate the distributions.
max_n_snps <- 10
m[, max_n_snps + 2] <- rowSums(m[, (max_n_snps + 2):ncol(m)], na.rm = T)
m <- m[, 1:(max_n_snps + 2)]
lat <- m / rowSums(m, na.rm = T)

# lib2, yol
snps_per_loc <- read.delim(here('data/test_libraries/lib2/yoldia/n_snps_per_locus_cat.tsv'))
snps_per_loc <- subset(snps_per_loc, M==n & m==3)
snps_per_loc <- snps_per_loc[order(snps_per_loc$n_snps),]
Mn_values <- sort(unique(snps_per_loc$M))
m <- matrix(NA, nrow = length(Mn_values), ncol = max(snps_per_loc$n_snps)+1)
for(i in 1:nrow(snps_per_loc)) {
  m[snps_per_loc$M[i], snps_per_loc$n_snps[i]+1] = snps_per_loc$n_loci[i] # [n_snps+1] as column 1 is for loci with 0 SNPs
}
max_n_snps <- 10
m[, max_n_snps + 2] <- rowSums(m[, (max_n_snps + 2):ncol(m)], na.rm = T)
m <- m[, 1:(max_n_snps + 2)]
yol <- m / rowSums(m, na.rm = T)

# call plot
mypal <- pal_npg("nrc", alpha = 0.7)(9)
colfunc <- colorRampPalette(c(mypal[1], 'lightgray'))
colfunc(length(Mn_values))
lat_plot <- function(x){
  barplot(lat,
        beside = T, col = colfunc(length(Mn_values)), las = 1,
        names.arg = c(0:max_n_snps, paste('>', max_n_snps, sep = '')),
        xlab = 'number of SNPs',
        ylab = 'percentage of loci'
  )
  legend(100, 0.5, legend = Mn_values, fill = colfunc(length(Mn_values)),
         title = "M and n", bty = 'n')
}
# yol
colfunc2 <- colorRampPalette(c(mypal[2], 'lightgray'))
yol_plot <- function(x){
  barplot(yol,
          beside = T, col = colfunc2(length(Mn_values)), las = 1,
          names.arg = c(0:max_n_snps, paste('>', max_n_snps, sep = '')),
          xlab = 'number of SNPs',
          ylab = 'percentage of loci'
  )
  legend(100, 0.5, legend = Mn_values, fill = colfunc2(length(Mn_values)),
         title = "M and n", bty = 'n')
}
printfig(here('figures/n_snps_lib2_lat'), lat_plot)
printfigpdf(here('figures/n_snps_lib2_lat'), lat_plot) 
printfig(here('figures/n_snps_lib2_yol'), yol_plot)
printfigpdf(here('figures/n_snps_lib2_yol'), yol_plot) 
rm(lat, yol, m, snps_per_loc, i, max_n_snps, mypal, Mn_values, colfunc, colfunc2, lat_plot, yol_plot)
#####

#### lib3
#####
# lib3, bathy
snps_per_loc <- read.delim(here('data/test_libraries/lib3/bathybiaster/n_snps_per_locus_cat.tsv'))
snps_per_loc <- subset(snps_per_loc, M==n & m==3)
snps_per_loc <- snps_per_loc[order(snps_per_loc$n_snps),]
Mn_values <- sort(unique(snps_per_loc$M))
m <- matrix(NA, nrow = length(Mn_values), ncol = max(snps_per_loc$n_snps)+1)
for(i in 1:nrow(snps_per_loc)) {
  m[snps_per_loc$M[i], snps_per_loc$n_snps[i]+1] = snps_per_loc$n_loci[i] # [n_snps+1] as column 1 is for loci with 0 SNPs
}
max_n_snps <- 10
m[, max_n_snps + 2] <- rowSums(m[, (max_n_snps + 2):ncol(m)], na.rm = T)
m <- m[, 1:(max_n_snps + 2)]
bathy <- m / rowSums(m, na.rm = T)

# lib3, psi
snps_per_loc <- read.delim(here('data/test_libraries/lib3/psilaster/n_snps_per_locus_cat.tsv'))
snps_per_loc <- subset(snps_per_loc, M==n & m==3)
snps_per_loc <- snps_per_loc[order(snps_per_loc$n_snps),]
Mn_values <- sort(unique(snps_per_loc$M))
m <- matrix(NA, nrow = length(Mn_values), ncol = max(snps_per_loc$n_snps)+1)
for(i in 1:nrow(snps_per_loc)) {
  m[snps_per_loc$M[i], snps_per_loc$n_snps[i]+1] = snps_per_loc$n_loci[i] # [n_snps+1] as column 1 is for loci with 0 SNPs
}
max_n_snps <- 10
m[, max_n_snps + 2] <- rowSums(m[, (max_n_snps + 2):ncol(m)], na.rm = T)
m <- m[, 1:(max_n_snps + 2)]
psi <- m / rowSums(m, na.rm = T)

# call plot
mypal <- pal_npg("nrc", alpha = 0.7)(9)
colfunc <- colorRampPalette(c(mypal[3], 'lightgray'))
colfunc(length(Mn_values))
bathy_plot <- function(x){
  barplot(bathy,
          beside = T, col = colfunc(length(Mn_values)), las = 1,
          names.arg = c(0:max_n_snps, paste('>', max_n_snps, sep = '')),
          xlab = 'number of SNPs',
          ylab = 'percentage of loci'
  )
  legend(70, max(bathy), legend = Mn_values, fill = colfunc(length(Mn_values)),
         title = "M and n", bty = 'n')
}
# psi
colfunc2 <- colorRampPalette(c(mypal[4], 'lightgray'))
psi_plot <- function(x){
  barplot(psi,
          beside = T, col = colfunc2(length(Mn_values)), las = 1,
          names.arg = c(0:max_n_snps, paste('>', max_n_snps, sep = '')),
          xlab = 'number of SNPs',
          ylab = 'percentage of loci'
  )
  legend(70, max(psi), legend = Mn_values, fill = colfunc2(length(Mn_values)),
         title = "M and n", bty = 'n')
}
printfig(here('figures/n_snps_lib3_bathy'), bathy_plot)
printfigpdf(here('figures/n_snps_lib3_bathy'), bathy_plot) 
printfig(here('figures/n_snps_lib3_psi'), psi_plot)
printfigpdf(here('figures/n_snps_lib3_psi'), psi_plot) 
rm(bathy, psi, m, snps_per_loc, i, max_n_snps, mypal, Mn_values, colfunc, colfunc2, bathy_plot, psi_plot)
#####

#### lib4
#####
# lib4, tber
snps_per_loc <- read.delim(here('data/test_libraries/lib4/tber/c/n_snps_per_locus.tsv'))
snps_per_loc <- subset(snps_per_loc, M==n & m==3)
snps_per_loc <- snps_per_loc[order(snps_per_loc$n_snps),]
Mn_values <- sort(unique(snps_per_loc$M))
m <- matrix(NA, nrow = length(Mn_values), ncol = max(snps_per_loc$n_snps)+1)
for(i in 1:nrow(snps_per_loc)) {
  m[snps_per_loc$M[i], snps_per_loc$n_snps[i]+1] = snps_per_loc$n_loci[i] # [n_snps+1] as column 1 is for loci with 0 SNPs
}
max_n_snps <- 10
m[, max_n_snps + 2] <- rowSums(m[, (max_n_snps + 2):ncol(m)], na.rm = T)
m <- m[, 1:(max_n_snps + 2)]
tber <- m / rowSums(m, na.rm = T)

# lib4, tloe
snps_per_loc <- read.delim(here('data/test_libraries/lib4/tloe/c/n_snps_per_locus.tsv'))
snps_per_loc <- subset(snps_per_loc, M==n & m==3)
snps_per_loc <- snps_per_loc[order(snps_per_loc$n_snps),]
Mn_values <- sort(unique(snps_per_loc$M))
m <- matrix(NA, nrow = length(Mn_values), ncol = max(snps_per_loc$n_snps)+1)
for(i in 1:nrow(snps_per_loc)) {
  m[snps_per_loc$M[i], snps_per_loc$n_snps[i]+1] = snps_per_loc$n_loci[i] # [n_snps+1] as column 1 is for loci with 0 SNPs
}
max_n_snps <- 10
m[, max_n_snps + 2] <- rowSums(m[, (max_n_snps + 2):ncol(m)], na.rm = T)
m <- m[, 1:(max_n_snps + 2)]
tloe <- m / rowSums(m, na.rm = T)

# call plot
mypal <- pal_npg("nrc", alpha = 0.7)(9)
colfunc <- colorRampPalette(c(mypal[5], 'lightgray'))
colfunc(length(Mn_values))
tber_plot <- function(x){
  barplot(tber,
          beside = T, col = colfunc(length(Mn_values)), las = 1,
          names.arg = c(0:max_n_snps, paste('>', max_n_snps, sep = '')),
          xlab = 'number of SNPs',
          ylab = 'percentage of loci'
  )
  legend(100, max(tber), legend = Mn_values, fill = colfunc(length(Mn_values)),
         title = "M and n", bty = 'n')
}
# tloe
colfunc2 <- colorRampPalette(c(mypal[6], 'lightgray'))
tloe_plot <- function(x){
  barplot(tloe,
          beside = T, col = colfunc2(length(Mn_values)), las = 1,
          names.arg = c(0:max_n_snps, paste('>', max_n_snps, sep = '')),
          xlab = 'number of SNPs',
          ylab = 'percentage of loci'
  )
  legend(100, max(tloe), legend = Mn_values, fill = colfunc2(length(Mn_values)),
         title = "M and n", bty = 'n')
}
printfig(here('figures/n_snps_lib4_tber'), tber_plot)
printfigpdf(here('figures/n_snps_lib4_tber'), tber_plot) 
printfig(here('figures/n_snps_lib4_tloe'), tloe_plot)
printfigpdf(here('figures/n_snps_lib4_tloe'), tloe_plot) 
rm(tber, tloe, m, snps_per_loc, i, max_n_snps, mypal, Mn_values, colfunc, colfunc2, tber_plot, tloe_plot)
#####

#### lib5
#####
# lib5
snps_per_loc <- read.delim(here('data/test_libraries/lib5/n_snps_per_locus.tsv'))
snps_per_loc <- subset(snps_per_loc, M==n & m==3)
snps_per_loc <- snps_per_loc[order(snps_per_loc$n_snps),]
Mn_values <- sort(unique(snps_per_loc$M))
m <- matrix(NA, nrow = length(Mn_values), ncol = max(snps_per_loc$n_snps)+1)
for(i in 1:nrow(snps_per_loc)) {
  m[snps_per_loc$M[i], snps_per_loc$n_snps[i]+1] = snps_per_loc$n_loci[i] # [n_snps+1] as column 1 is for loci with 0 SNPs
}
max_n_snps <- 10
m[, max_n_snps + 2] <- rowSums(m[, (max_n_snps + 2):ncol(m)], na.rm = T)
m <- m[, 1:(max_n_snps + 2)]
lib5 <- m / rowSums(m, na.rm = T)

# call plot
mypal <- pal_npg("nrc", alpha = 0.7)(9)
colfunc <- colorRampPalette(c(mypal[8], 'lightgray'))
colfunc(length(Mn_values))
lib5_plot <- function(x){
  barplot(lib5,
          beside = T, col = colfunc(length(Mn_values)), las = 1,
          names.arg = c(0:max_n_snps, paste('>', max_n_snps, sep = '')),
          xlab = 'number of SNPs',
          ylab = 'percentage of loci'
  )
  legend(100, 0.7, legend = Mn_values, fill = colfunc(length(Mn_values)),
         title = "M and n", bty = 'n')
}
printfig(here('figures/n_snps_lib5'), lib5_plot)
printfigpdf(here('figures/n_snps_lib5'), lib5_plot) 
rm(lib5, m, snps_per_loc, i, max_n_snps, mypal, Mn_values, colfunc, lib5_plot)
#####
