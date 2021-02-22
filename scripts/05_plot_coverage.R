#### Script for plotting coverage results of
#### test RAD libraries
## 28/06/2019
## H. Christiansen
## v1.4

#### load packages
library(here) # to shorten file paths
library(tidyverse) # for fancy plots and data re-organization
library(ggsci) # for journal type colors
library(scales) # to visualize ggsci palettes

#### coverage, all
#####
coverage <- read.csv(here("data/test_libraries/coverage_stats.csv"), header = T)
coverage$opt_M <- as.factor(coverage$opt_M)

mypal <- pal_npg("nrc", alpha = 0.7)(9)
show_col(mypal)
mypal <- c(mypal[7], mypal[1], mypal[2], mypal[3], mypal[4], mypal[5], mypal[6], mypal[8])

p1 <-   coverage %>%
        mutate(species = fct_reorder(species, lib)) %>%
          ggplot(aes(x = factor(lib), y = cov_opt_M, fill = forcats::fct_reorder(factor(species), lib) )) +
            geom_boxplot(varwidth = T) + 
            theme_bw() +
            ylim(0, 70) +
            labs(fill = "") +
            theme(legend.position="top") +
            theme(legend.text = element_text(size = 10, face = "italic")) +
            scale_fill_manual(values = mypal) +
            ylab("Coverage") + xlab("Library") +
            theme(text = element_text(size = 20), axis.text.x = element_text(size = 20)) +
            geom_point(aes(x = lib, y = target_cov), color = "darkred", shape = 18, size =6)
p1
ggsave("Coverage.jpg", path = here("figures"),
       width = 7.25, height = 5, units = "in", dpi = 1200)
ggsave("Coverage.pdf", path = here("figures"),
       width = 7.25, height = 5, units = "in", dpi = 1200)
ggsave("Coverage.png", path = here("figures"),
       width = 7.25, height = 5, units = "in", dpi = 1200)
dev.off()
rm(coverage, mypal, p1)
#####

#### loci, all
#####
coverage <- read.csv(here("data/test_libraries/coverage_stats.csv"), header = T)

mypal <- pal_npg("nrc", alpha = 0.7)(9)
show_col(mypal)
mypal <- c(mypal[7], mypal[1], mypal[2], mypal[3], mypal[4], mypal[5], mypal[6], mypal[8])

p2 <-   coverage %>%
  mutate(species = fct_reorder(species, lib)) %>%
  ggplot(aes(x = factor(lib), y = n_loci/1000, fill = forcats::fct_reorder(factor(species), lib) )) +
  geom_boxplot(varwidth = T) + 
  theme_bw() +
  scale_fill_manual(values = mypal) +
  ylab("Loci (thousands)") + xlab("Library") +
  labs(fill = "") +
  theme(legend.position="top") +
  theme(legend.text = element_text(size = 10, face = "italic")) +
  theme(text = element_text(size = 20), axis.text.x = element_text(size = 20)) +
  geom_point(aes(x = lib, y = target_loci/1000), color = "darkred", shape = 18, size =6)
p2
ggsave("Loci.jpg", path = here("figures"),
       width = 7.25, height = 5, units = "in", dpi = 1200)
ggsave("Loci.pdf", path = here("figures"),
       width = 7.25, height = 5, units = "in", dpi = 1200)
ggsave("Loci.png", path = here("figures"),
       width = 7.25, height = 5, units = "in", dpi = 1200)
dev.off()
rm(coverage, mypal, p2)
#####
