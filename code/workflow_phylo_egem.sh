### ===========================================================================================
### bash code for haplotype network and phylogenetic analysis for the manuscript
### "A plea for a conservative approach to the report of new species records from eDNA surveys"
### by Oscar Puebla and Martin Helmkampf (2022)
### ===========================================================================================

### 0. Preparation

# Append project file structure
# egem
# ...
# |- 3_phylo
#    |- log
#    |- raxml

## Download Hypoplectrus gemma mitochondrion, complete genome (GenBank accession FJ848375)
#> Hgem_mtgenome_FJ848375.fas


### ====================================
### 1. Extraction of 12S teleo sequences

#!/bin/bash

#SBATCH --job-name=1_fasta
#SBATCH --partition=carl.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00  # D-HH:MM
#SBATCH --output=log/1_fasta_%j.out
#SBATCH --error=log/1_fasta_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-uniol-env
ml VCFtools/0.1.14

BASE_DIR=/gss/work/haex1482/2_Other/egem/3_phylo

vcftools \
    --gzvcf $BASE_DIR/../2_genotyping/out/8_geno/egem_all_LGM.vcf.gz \
    --chr LG_M \
    --from-bp 2682 \
    --to-bp 2745 \
    --recode \
    --stdout > $BASE_DIR/egem_12steleo.vcf

vcf-to-tab < \
    $BASE_DIR/egem_12steleo.vcf |
    sed -e 's/\.\/\./N\/N/g' -e 's/[ACGTN\*]\/\*/N\/N/g' -e 's/\*\/[ACGTN\*]/N\/N/g' > \
    $BASE_DIR/egem_12steleo.tab

perl ~/apps/vcf-tab-to-fasta/vcf_tab_to_fasta_alignment.pl \
    -i $BASE_DIR/egem_12steleo.tab > \
    $BASE_DIR/egem_12steleo.fas

rm $BASE_DIR/12s/egem_12steleo.tab*


# Add 12S teleo sequences of reference and gemma mitogenome manually
cat egem_12steleo.fas HpueHgem_12steleo.fas > egem_12steleo2.fas

# Confirm sequences are of equal length (aligned)
bioawk -c fastx '{ print $name, length($seq) }' egem_12steleo2.fas


### 2. Network analysis of 12S teleo

# egem_12steleo2.fas: 126 taxa and 64 sites, 5 patterns, 124 samples identical, 98.4% invariant sites

# see plot_haplonet.R


### 3. Mitogenome phylogeny: Fasta conversion

#!/bin/bash

#SBATCH --job-name=3_fasta
#SBATCH --partition=carl.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00  # D-HH:MM
#SBATCH --output=log/3_fasta_%j.out
#SBATCH --error=log/3_fasta_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-uniol-env
ml VCFtools/0.1.14

BASE_DIR=/gss/work/haex1482/2_Other/egem/3_phylo

gzip -cd $BASE_DIR/../2_genotyping/out/8_geno/egem_all_LGM.vcf.gz > $BASE_DIR/egem_mtg.vcf

sed -e 's/|/\//g' \
    $BASE_DIR/egem_mtg.vcf |
    grep -v "##" > \
    $BASE_DIR/egem_mtg_n.vcf

vcf-to-tab < \
    $BASE_DIR/egem_mtg_n.vcf |
    sed -e 's/\.\/\./N\/N/g' -e 's/[ACGTN\*]\/\*/N\/N/g' -e 's/\*\/[ACGTN\*]/N\/N/g' > \
    $BASE_DIR/egem_mtg.tab

perl ~/apps/vcf-tab-to-fasta/vcf_tab_to_fasta_alignment.pl \
    -i $BASE_DIR/egem_mtg.tab > \
    $BASE_DIR/egem_mtg.fas

rm $BASE_DIR/egem_mtg.vcf
rm $BASE_DIR/egem_mtg_n.vcf
rm $BASE_DIR/egem_mtg.tab*

# Confirm sequences are of equal length (aligned)
bioawk -c fastx '{ print $name, length($seq) }' egem_mtg.fas

# 17030 sites


### 4. Mitogenome phylogeny: Alignment and model test

# Reverse complement Hgem mitogenome (in 3_phylo)
~/apps/seqtk/seqtk seq -r Hgem_mtgenome_FJ848375.fas > Hgem_mtgenome_FJ848375_rc.fas   # from 1_rawdata/1_mt

# Blast to determine offset: Hpue 17159 (end) is Hgem 13485, manually moved sequence after that to front
#> Hgem_mtgenome_FJ848375_rcre.fas

cat egem_mtg.fas ../1_rawdata/1_mt/Hpue_LGM_01.fasta ../1_rawdata/1_mt/Hgem_mtgenome_FJ848375_rcre.fas | \
    sed -e 's/LG_M/Ref_puemtg/g' -e 's/ENA.*/ENA_gemmtg/g' > egem_mtg2.fas

#!/bin/bash

#SBATCH --job-name=4_aln
#SBATCH --partition=carl.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00  # D-HH:MM
#SBATCH --output=log/4_aln_%j.out
#SBATCH --error=log/4_aln_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml MAFFT/7.475-GCC-8.3.0-with-extensions

BASE_DIR=/gss/work/haex1482/2_Other/egem/3_phylo

mafft --auto egem_mtg2.fas > egem_mtg2_f.aln
#> FFT-NS-2 (Fast)

java -jar ~/apps/jmodeltest_2.1.10/jModelTest.jar -d egem_mtg2_f.aln -g 4 -i -f -AIC -BIC -a > egem_mtg2_f.mtest
#> GTR+I+G (AIC, BIC)

# Alignment: 126 taxa and 17174 sites, 1037 patterns, 0 samples identical, 89.8% invariant sites


### 5. Mitogenome phylogeny: phylogenetic reconstruction

#!/bin/bash

#SBATCH --job-name=5_phylo
#SBATCH --partition=carl.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=1-00:00  # D-HH:MM
#SBATCH --output=log/5_phylo_%j.out
#SBATCH --error=log/5_phylo_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

BASE_DIR=/gss/work/haex1482/2_Other/egem/3_phylo

raxml-ng \
    --all \
    --msa $BASE_DIR/egem_mtg2_f.aln \
    --model GTR+G \
    --tree pars{25},rand{25} \
    --bs-trees 500 \
    --threads 24 \
    --worker AUTO \
    --seed 123 \
    --prefix $BASE_DIR/raxml/egem_mtg2_f_gtr_500
