### ===========================================================================================
### R code to calculate mitochondrial genome distances for 
### "A plea for a conservative approach to the report of new species records from eDNA surveys"
### by Oscar Puebla and Martin Helmkampf (2022)
### ===========================================================================================

library(ape)
library(seqinr)

setwd("/Users/martin/Documents/Projects/2_Other/egem/3_phylo/network")


### Read in fasta file
sequences <- read.dna("../egem_mtg2_f.aln", format = "fasta")


### Calculate distance matrix
dmatrix <- as.matrix(dist.dna(sequences, model = "raw"))


### Calculate mean distance between Gulf and Caribbean clades
dclades <- dmatrix[!rownames(dmatrix) %in% c("PL17_160floflo", "ENA_gemmtg"),
                    colnames(dmatrix) %in% c("PL17_160floflo", "ENA_gemmtg")]
mean(dclades)


### Calculate mean distance within Caribbean clade
dcarib <- dmatrix[!rownames(dmatrix) %in% c("PL17_160floflo", "ENA_gemmtg"),
                  !colnames(dmatrix) %in% c("PL17_160floflo", "ENA_gemmtg")]
mean(dcarib)


### Check for bias from assemblies by removing them
# bmatrix <- dmatrix[!rownames(dmatrix) %in% c("Ref_puemtg", "ENA_gemmtg"), 
#                    !colnames(dmatrix) %in% c("Ref_puemtg", "ENA_gemmtg")]
# 
# 
# bclades <- bmatrix[rownames(bmatrix) != "PL17_160floflo",
#                    colnames(bmatrix) == "PL17_160floflo"]
# mean(bclades)
# 
# 
# bcarib <- bmatrix[rownames(bmatrix) != "PL17_160floflo",
#                   colnames(bmatrix) != "PL17_160floflo"]
# mean(bcarib)
