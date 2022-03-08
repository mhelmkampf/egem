__Description__

This repository contains code and data files that were used for the manuscript 
"*A plea for a conservative approach to the report of new species records from eDNA surveys*" by Oscar Puebla and Martin Helmkampf (submitted 2022).


__Contents__
```
.
├── README.md
├── code
│   ├── mtg_distances.R                     # calculate mt genome genetic distances
│   ├── plot_haplonet.R                     # plot 12S teleo haplotype network
│   ├── plot_tree_mtg.R                     # plot phylogenetic tree of mitochondrial genomes
│   ├── prepare_meta.R                      # prepare metadata file for genotyping
│   ├── workflow_genotyping_egem.sh         # genotyping pipeline
│   └── workflow_phylo_egem.sh              # network and phylogenetic analysis workflow
└── data
    ├── egem_12steleo2.fas                  # 12S teleo sequences
    ├── egem_labels.tsv                     # sample labels, input file for prepare_meta.R                        
    ├── egem_mtg2_f.aln                     # mitochondrial genome alignment
    └── egem_mtg2_f_gtr_500.raxml.support   # RAxML best ML tree with bootstrap supports
```
