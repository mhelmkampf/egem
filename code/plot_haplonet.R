### ===========================================================================================
### R code to plot the 12S teleo haplotype network (Figure 1a) in 
### "A plea for a conservative approach to the report of new species records from eDNA surveys"
### by Oscar Puebla and Martin Helmkampf (2022)
### ===========================================================================================

library(pegas)
library(tidyverse)

setwd("/Users/martin/Documents/Projects/2_Other/egem/3_phylo/network")


### Read in sequences and define species / populations
(sequences <- read.dna("../egem_12steleo2.fas", format = "fasta"))

(species <- str_replace(labels(sequences), ".*([a-z]{3})[a-z]{3}", "\\1"))

(pop <- str_replace(labels(sequences), ".*([a-z]{6})", "\\1") %>%
  str_replace("puemtg", "puepan") %>%
  str_replace("gemmtg", "gemflo"))


### Extract haplotypes
(haplotypes <- haplotype(sequences))


### Compute network
(network <- haploNet(haplotypes))


### Size ~ frequency
# (sz <- summary(haplotypes))
# (net.labs <- attr(network, "labels"))
# sz <- sz[net.labs]   # align order with network


### Plot basic network
plot(network, size = attr(network, "freq"), fast = FALSE)   # plot(network, size = sz, fast = FALSE)


### Define pie segments for species / population
slist <- c("flo", "gum", "gem", "may", "nig", "pue", "uni")
S <- haploFreq(sequences, fac = species, haplo = haplotypes)
(S <- S[, slist])

plist <- c("floflo", "gumpan", "gemflo", "maybel", "nigbel", "nighon", "nigpan", "puebel", "puehon", "puepan", "unibel", "unihon", "unipan")
P <- haploFreq(sequences, fac = pop, haplo = haplotypes)
(P <- P[, plist])


### Set colors and plotting options (workaround)
fs <- function(n) c("#A39CCF", "#E3A258", "dodgerblue", "#7EA7C2", "gray20", "#E17366", "#CCCCCC")

fp <- function(n) c("#A39CCF",
                    "#E3A258",
                    "dodgerblue",
                    "#7EA7C2", 
                    "gray25", "gray25", "gray25",
                    "#E17366", "#E17366", "#E17366",
                    "#CCCCCC", "#CCCCCC", "#CCCCCC")

setHaploNetOptions(pie.inner.segments.color = "gray15")


### Plot network dolor-coded by species / population
(gs <- plot(network, 
            size = attr(network, "freq"), 
            pie = S, 
            bg = fs, 
            scale.ratio = 20, 
            labels = FALSE, 
            show.mutation = 1))

(gp <- plot(network, 
            size = attr(network, "freq"), 
            pie = P, 
            bg = fp, 
            scale.ratio = 20, 
            labels = FALSE, 
            show.mutation = 1))


### Save to file (currently bugged, save manually as PDF, default settings)
ggsave(plot = gp,
       filename = "egem_network_p.pdf",
       width = 6,
       height = 5,
       device = cairo_pdf,
       bg = "transparent")


### *** Note ***
### Plotting haplotype network object creates a list that could not be converted into a ggplot object
### Haplotype network and phylogenetic tree plots were combined manually in Inkscape
