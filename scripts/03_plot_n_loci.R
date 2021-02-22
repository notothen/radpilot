#### Script to plot the number of r80 loci
#### that can be constructed de novo with
#### varying M parameter from test RAD libraries
## 14/12/2019
## v1.6
## H. Christiansen,
## modified from Rochette & Catchen 2017:
browseURL("https://bitbucket.org/rochette/rad-seq-genotyping-demo/src/default/demo_scripts/R_scripts/")

#### load packages
library(here) # to shorten file paths
library(ggsci) # for journal type colors
library(scales) # to visualize ggsci palettes
source(here("scripts/printfig.R")) # to save plots

#### lib 1
#####
# read in data
snps_per_loc <- read.delim(here('data/test_libraries/lib1/n_snps_per_locus_r50.tsv'))
# Keep only M==n, m==3
snps_per_loc <- subset(snps_per_loc, M==n & m==3)
# Rename column 1
colnames(snps_per_loc)[1] <- 'par_set'

# Create a new data frame to contain the number of loci and polymorphic loci
d <- snps_per_loc[, c('par_set', 'M', 'n', 'm')]
d <- d[!duplicated(d), ]

# Compute these numbers for each parameter set, using the par_set column as an ID
rownames(d) <- d$par_set
for(p in rownames(d)) {
  s <- subset(snps_per_loc, par_set == p)
  d[p, 'n_loci'] <- sum(s$n_loci)
  s2 <- subset(s, n_snps > 0)
  d[p, 'n_loci_poly'] <- sum(s2$n_loci)
}

# Make sure the table is ordered
d <- d[order(d$M), ]

# call plot
mypal = pal_npg("nrc", alpha = 0.7)(9)
show_col(mypal)
lib1_plot <- function(x){
  plot(NULL,
       xlim = range(c(0, 9)),
       ylim = range(c(0, max(c(d$n_loci, d$n_loci))+500)),
       xlab = 'de novo assembly parameter M and n',
       ylab = 'loci shared by 50 % of samples',
       yaxt = 'n',
       xaxt = 'n',
       las = 2,
       bty = 'n'
  )
  axis(1, at = c(0:10), pos = 0)
  axis(2, at = seq(0, max(c(d$n_loci, d$n_loci))+500, by = 500),
       labels = formatC(seq(0, max(c(d$n_loci, d$n_loci))+500, by = 500), big.mark = ',', format = 'd'),
       pos = 0, las = 2)
  clip(0, 9.5, 0, max(c(d$n_loci, d$n_loci))+500)
  abline(h = 0:20*500, lty = 'dotted', col = 'grey50')
  lines(d$M, d$n_loci, lwd = 2, col = mypal[7])
  points(d$M, d$n_loci, lwd = 2, pch = 1, col = mypal[7])
  lines(d$M, d$n_loci_poly, lwd = 2, lty = 'dashed', col = mypal[7])
  points(d$M, d$n_loci_poly, lwd = 2, pch = 2, col = mypal[7])
  legend(5, 1000, c('Macrocyprididae all', 'Macrocyprididae polymorphic'),
         lty = c('solid', 'dashed'), lwd = 2, bg = 'white',
         pch = c(1, 2), col = c(mypal[7], mypal[7]))
}
printfig(here('figures/n_loci_Mn_lib1_r50'), lib1_plot)
printfigpdf(here('figures/n_loci_Mn_lib1_r50'), lib1_plot) 
dev.off()
rm(d, mypal, p, s, s2, snps_per_loc, lib1_plot)
#####

#### lib2
#####
# lib2, lat
snps_per_loc <- read.delim(here('data/test_libraries/lib2/laternula/n_snps_per_locus_cat.tsv'))
snps_per_loc <- subset(snps_per_loc, M==n & m==3)
colnames(snps_per_loc)[1] <- 'par_set'
d <- snps_per_loc[, c('par_set', 'M', 'n', 'm')]
d <- d[!duplicated(d), ]
rownames(d) <- d$par_set
for(p in rownames(d)) {
  s <- subset(snps_per_loc, par_set == p)
  d[p, 'n_loci'] <- sum(s$n_loci)
  s2 <- subset(s, n_snps > 0)
  d[p, 'n_loci_poly'] <- sum(s2$n_loci)
}
lat <- d[order(d$M), ]

# lib2, yol
snps_per_loc = read.delim(here('data/test_libraries/lib2/yoldia/n_snps_per_locus_cat.tsv'))
snps_per_loc <- subset(snps_per_loc, M==n & m==3)
colnames(snps_per_loc)[1] <- 'par_set'
d <- snps_per_loc[, c('par_set', 'M', 'n', 'm')]
d <- d[!duplicated(d), ]
rownames(d) <- d$par_set
for(p in rownames(d)) {
  s <- subset(snps_per_loc, par_set == p)
  d[p, 'n_loci'] <- sum(s$n_loci)
  s2 <- subset(s, n_snps > 0)
  d[p, 'n_loci_poly'] <- sum(s2$n_loci)
}
yol <- d[order(d$M), ]

# call plot
mypal = pal_npg("nrc", alpha = 0.7)(9)
lib2_plot <- function(x){
  plot(NULL,
       xlim = range(c(0, 9)),
       ylim = range(c(0, max(c(lat$n_loci, yol$n_loci))+5000)),
       xlab = 'de novo assembly parameter M and n',
       ylab = 'loci shared by 80 % of samples',
       yaxt = 'n',
       xaxt = 'n',
       las = 2,
       bty = 'n'
  )
  axis(1, at = c(0:10), pos = 0)
  axis(2, at = seq(0, max(c(lat$n_loci, yol$n_loci))+5000, by = 5000),
       labels = formatC(seq(0, max(c(lat$n_loci, yol$n_loci))+5000, by = 5000), big.mark = ',', format = 'd'),
       pos = 0, las = 2)
  clip(0, 9.5, 0, max(c(lat$n_loci, yol$n_loci))+5000)
  abline(h = 0:20*5000, lty = 'dotted', col = 'grey50')
  lines(lat$M, lat$n_loci, lwd = 2, col = mypal[1])
  points(lat$M, lat$n_loci, lwd = 2, pch = 1, col = mypal[1])
  lines(lat$M, lat$n_loci_poly, lwd = 2, lty = 'dashed', col = mypal[1])
  points(lat$M, lat$n_loci_poly, lwd = 2, pch = 2, col = mypal[1])
  lines(yol$M, yol$n_loci, lwd = 2, col = mypal[2])
  points(yol$M, yol$n_loci, lwd = 2, col = mypal[2])
  lines(yol$M, yol$n_loci_poly, lwd = 2, lty = 'dashed', col = mypal[2])
  points(yol$M, yol$n_loci_poly, lwd = 2, pch =2, col = mypal[2])
  legend(5, 20000, c('L. elliptica all', 'L. elliptica polymorphic', 'A. eightsii all', 'A. eightsii polymorphic'),
         lty = c('solid', 'dashed', 'solid', 'dashed'), lwd = 2, bg = 'white',
         pch = c(1, 2, 1, 2), col = c(mypal[1], mypal[1], mypal[2], mypal[2]))
}
printfig(here('figures/n_loci_Mn_lib2'), lib2_plot)
printfigpdf(here('figures/n_loci_Mn_lib2'), lib2_plot) 
dev.off()
rm(d, mypal, p, s, s2, snps_per_loc, lat, yol, lib2_plot)
#####

#### lib3
#####
# lib3, bathy
snps_per_loc <- read.delim(here('data/test_libraries/lib3/bathybiaster/n_snps_per_locus_cat.tsv'))
snps_per_loc <- subset(snps_per_loc, M==n & m==3)
colnames(snps_per_loc)[1] <- 'par_set'
d <- snps_per_loc[, c('par_set', 'M', 'n', 'm')]
d <- d[!duplicated(d), ]
rownames(d) <- d$par_set
for(p in rownames(d)) {
  s <- subset(snps_per_loc, par_set == p)
  d[p, 'n_loci'] <- sum(s$n_loci)
  s2 <- subset(s, n_snps > 0)
  d[p, 'n_loci_poly'] <- sum(s2$n_loci)
}
bathy <- d[order(d$M), ]

# lib3, psi
snps_per_loc = read.delim(here('data/test_libraries/lib3/psilaster/n_snps_per_locus_cat.tsv'))
snps_per_loc <- subset(snps_per_loc, M==n & m==3)
colnames(snps_per_loc)[1] <- 'par_set'
d <- snps_per_loc[, c('par_set', 'M', 'n', 'm')]
d <- d[!duplicated(d), ]
rownames(d) <- d$par_set
for(p in rownames(d)) {
  s <- subset(snps_per_loc, par_set == p)
  d[p, 'n_loci'] <- sum(s$n_loci)
  s2 <- subset(s, n_snps > 0)
  d[p, 'n_loci_poly'] <- sum(s2$n_loci)
}
psi <- d[order(d$M), ]

# call plot
mypal = pal_npg("nrc", alpha = 0.7)(9)
show_col(mypal)
lib3_plot <- function(x){
  plot(NULL,
       xlim = range(c(0, 9)),
       ylim = range(c(0, max(c(bathy$n_loci, psi$n_loci))+5000)),
       xlab = 'de novo assembly parameter M and n',
       ylab = 'loci shared by 80 % of samples',
       yaxt = 'n',
       xaxt = 'n',
       las = 2,
       bty = 'n'
  )
  axis(1, at = c(0:10), pos = 0)
  axis(2, at = seq(0, max(c(bathy$n_loci, psi$n_loci))+5000, by = 5000),
       labels = formatC(seq(0, max(c(bathy$n_loci, psi$n_loci))+5000, by = 5000), big.mark = ',', format = 'd'),
       pos = 0, las = 2)
  clip(0, 9.5, 0, max(c(bathy$n_loci, psi$n_loci))+5000)
  abline(h = 0:20*5000, lty = 'dotted', col = 'grey50')
  lines(bathy$M, bathy$n_loci, lwd = 2, col = mypal[3])
  points(bathy$M, bathy$n_loci, lwd = 2, pch = 1, col = mypal[3])
  lines(bathy$M, bathy$n_loci_poly, lwd = 2, lty = 'dashed', col = mypal[3])
  points(bathy$M, bathy$n_loci_poly, lwd = 2, pch = 2, col = mypal[3])
  lines(psi$M, psi$n_loci, lwd = 2, col = mypal[4])
  points(psi$M, psi$n_loci, lwd = 2, col = mypal[4])
  lines(psi$M, psi$n_loci_poly, lwd = 2, lty = 'dashed', col = mypal[4])
  points(psi$M, psi$n_loci_poly, lwd = 2, pch =2, col = mypal[4])
  legend(5, 22000, c('B. loripes all', 'B. loripes polymorphic', 'P. charcoti all', 'P. charcoti polymorphic'),
         lty = c('solid', 'dashed', 'solid', 'dashed'), lwd = 2, bg = 'white',
         pch = c(1, 2, 1, 2), col = c(mypal[3], mypal[3], mypal[4], mypal[4]))
}
printfig(here('figures/n_loci_Mn_lib3'), lib3_plot)
printfigpdf(here('figures/n_loci_Mn_lib3'), lib3_plot) 
dev.off()
rm(d, mypal, p, s, s2, snps_per_loc, bathy, psi, lib3_plot)
#####

#### lib4
#####
# lib4, tber
snps_per_loc <- read.delim(here('data/test_libraries/lib4/tber/n_snps_per_locus.tsv'))
snps_per_loc <- subset(snps_per_loc, M==n & m==3)
colnames(snps_per_loc)[1] <- 'par_set'
d <- snps_per_loc[, c('par_set', 'M', 'n', 'm')]
d <- d[!duplicated(d), ]
rownames(d) <- d$par_set
for(p in rownames(d)) {
  s <- subset(snps_per_loc, par_set == p)
  d[p, 'n_loci'] <- sum(s$n_loci)
  s2 <- subset(s, n_snps > 0)
  d[p, 'n_loci_poly'] <- sum(s2$n_loci)
}
tber <- d[order(d$M), ]

# lib4, tloe
snps_per_loc = read.delim(here('data/test_libraries/lib4/tloe/n_snps_per_locus.tsv'))
snps_per_loc <- subset(snps_per_loc, M==n & m==3)
colnames(snps_per_loc)[1] <- 'par_set'
d <- snps_per_loc[, c('par_set', 'M', 'n', 'm')]
d <- d[!duplicated(d), ]
rownames(d) <- d$par_set
for(p in rownames(d)) {
  s <- subset(snps_per_loc, par_set == p)
  d[p, 'n_loci'] <- sum(s$n_loci)
  s2 <- subset(s, n_snps > 0)
  d[p, 'n_loci_poly'] <- sum(s2$n_loci)
}
tloe <- d[order(d$M), ]

# call plot
mypal = pal_npg("nrc", alpha = 0.7)(9)
show_col(mypal)
lib4_plot <- function(x){
  plot(NULL,
     xlim = range(c(0, 9)),
     ylim = range(c(0, max(c(tber$n_loci, tloe$n_loci))+2000)),
     xlab = 'de novo assembly parameter M and n',
     ylab = 'loci shared by 80 % of samples',
     yaxt = 'n',
     xaxt = 'n',
     las = 2,
     bty = 'n'
  )
  axis(1, at = c(0:10), pos = 0)
  axis(2, at = seq(0, max(c(tber$n_loci, tloe$n_loci))+2000, by = 2000),
       labels = formatC(seq(0, max(c(tber$n_loci, tloe$n_loci))+2000, by = 2000), big.mark = ',', format = 'd'),
       pos = 0, las = 2)
  clip(0, 9.5, 0, max(c(tber$n_loci, tloe$n_loci))+2000)
  abline(h = 0:20*2000, lty = 'dotted', col = 'grey50')
  lines(tber$M, tber$n_loci, lwd = 2, col = mypal[5])
  points(tber$M, tber$n_loci, lwd = 2, pch = 1, col = mypal[5])
  lines(tber$M, tber$n_loci_poly, lwd = 2, lty = 'dashed', col = mypal[5])
  points(tber$M, tber$n_loci_poly, lwd = 2, pch = 2, col = mypal[5])
  lines(tloe$M, tloe$n_loci, lwd = 2, col = mypal[6])
  points(tloe$M, tloe$n_loci, lwd = 2, col = mypal[6])
  lines(tloe$M, tloe$n_loci_poly, lwd = 2, lty = 'dashed', col = mypal[6])
  points(tloe$M, tloe$n_loci_poly, lwd = 2, pch = 2, col = mypal[6])
  legend(5, 6000, c('T. bernacchii all', 'T. bernacchii polymorphic', 'T. loennbergii all', 'T. loennbergii polymorphic'),
         lty = c('solid', 'dashed', 'solid', 'dashed'), lwd = 2, bg = 'white',
         pch = c(1, 2, 1, 2), col = c(mypal[5], mypal[5], mypal[6], mypal[6]))
}
printfig(here('figures/n_loci_Mn_lib4'), lib4_plot)
printfigpdf(here('figures/n_loci_Mn_lib4'), lib4_plot) 
dev.off()
rm(d, mypal, p, s, s2, snps_per_loc, tber, tloe, lib4_plot)
#####

#### lib 5
#####
# read in data
snps_per_loc <- read.delim(here('data/test_libraries/lib5/n_snps_per_locus.tsv'))
# Keep only M==n, m==3
snps_per_loc <- subset(snps_per_loc, M==n & m==3)
# Rename column 1
colnames(snps_per_loc)[1] <- 'par_set'

# Create a new data frame to contain the number of loci and polymorphic loci
d <- snps_per_loc[, c('par_set', 'M', 'n', 'm')]
d <- d[!duplicated(d), ]

# Compute these numbers for each parameter set, using the par_set column as an ID
rownames(d) <- d$par_set
for(p in rownames(d)) {
  s <- subset(snps_per_loc, par_set == p)
  d[p, 'n_loci'] <- sum(s$n_loci)
  s2 <- subset(s, n_snps > 0)
  d[p, 'n_loci_poly'] <- sum(s2$n_loci)
}

# Make sure the table is ordered
d <- d[order(d$M), ]

# call plot
mypal = pal_npg("nrc", alpha = 0.7)(9)
show_col(mypal)
lib5_plot <- function(x){
  plot(NULL,
       xlim = range(c(0, 9)),
       ylim = range(c(0, max(c(d$n_loci, d$n_loci))+10000)),
       xlab = 'de novo assembly parameter M and n',
       ylab = 'loci shared by 80 % of samples',
       yaxt = 'n',
       xaxt = 'n',
       las = 2,
       bty = 'n'
  )
  axis(1, at = c(0:10), pos = 0)
  axis(2, at = seq(0, max(c(d$n_loci, d$n_loci))+10000, by = 10000),
       labels = formatC(seq(0, max(c(d$n_loci, d$n_loci))+10000, by = 10000), big.mark = ',', format = 'd'),
       pos = 0, las = 2)
  clip(0, 9.5, 0, max(c(d$n_loci, d$n_loci))+10000)
  abline(h = 0:20*10000, lty = 'dotted', col = 'grey50')
  lines(d$M, d$n_loci, lwd = 2, col = mypal[8])
  points(d$M, d$n_loci, lwd = 2, pch = 1, col = mypal[8])
  lines(d$M, d$n_loci_poly, lwd = 2, lty = 'dashed', col = mypal[8])
  points(d$M, d$n_loci_poly, lwd = 2, pch = 2, col = mypal[8])
  legend(5, 40000, c('P. nivea all', 'P. nivea polymorphic'),
         lty = c('solid', 'dashed'), lwd = 2, bg = 'white',
         pch = c(1, 2), col = c(mypal[8], mypal[8]))
}
printfig(here('figures/n_loci_Mn_lib5'), lib5_plot)
printfigpdf(here('figures/n_loci_Mn_lib5'), lib5_plot) 
dev.off()
rm(d, mypal, p, s, s2, snps_per_loc, lib5_plot)
#####

