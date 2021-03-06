### =================================================================================================
### bash code to obtain mitochondrial genomes in VCF format from resequencing data for the manuscript
### "A plea for a conservative approach to the report of new species records from eDNA surveys"
### by Oscar Puebla and Martin Helmkampf (2022)
### =================================================================================================

### 0. Preparation

# Create project file structure
# egem
# |- 0_metadata
# |- 1_rawdata
# |- 2_genotyping
#    |- fof
#    |- log
#       |- 1_ubam
#       |- 2_adap
#       ...
#       |- 8_geno
#    |- out
#       |- (same as log)
#    |- ref
#    |- tmp

## Download raw data from ENA to 1_rawdata (124 x 2 files)

## Compile metadata file (in 1_rawdata)
for i in *.fastq.gz ; do echo $i $(zcat < $i | head -n 1) ; done
find *q.gz > ../0_metadata/raw_data.fof
# Proceed with R script: prepare_meta.R, using egem_labels.tsv and raw_data.fof
#> meta_egem.csv

## Extract and index mitogenome from reference genome assembly (GenBank accession GCA_900610375.1, in 2_genotyping/ref):
grep -A 1 '>LG_M' Hpue_genome_unmasked_01.fas > Hpue_LGM_01.fasta

ml hpc-env/8.3 GATK/4.1.9.0-GCCcore-8.3.0-Java-8 SAMtools/1.9-GCC-8.3.0 BWA/0.7.17-GCC-8.3.0

bwa index Hpue_LGM_01.fasta
gatk CreateSequenceDictionary -R Hpue_LGM_01.fasta
samtools faidx Hpue_LGM_01.fasta


### =================================
### 1. Convert Fastq to unaligned BAM

#!/bin/bash

#SBATCH --job-name=1_ubam
#SBATCH --partition=carl.p
#SBATCH --array=1-126   # 124 samples, 2 with 2 datasets
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=0-03:00  # D-HH:MM
#SBATCH --output=log/1_ubam/1_ubam_%A_%a.out
#SBATCH --error=log/1_ubam/1_ubam_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8

BASE_DIR=/gss/work/haex1482/2_Other/egem/2_genotyping
META_FILE=$BASE_DIR/../0_metadata/egem_meta.csv
LINES=$(head $META_FILE -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=";" read SAMPLE FWD REV LANEF <<< $LINES

gatk --java-options "-Xmx20G" \
    FastqToSam \
    -SM $SAMPLE \
    -F1 $BASE_DIR/../1_rawdata/$FWD \
    -F2 $BASE_DIR/../1_rawdata/$REV \
    -RG ${SAMPLE}.${LANEF} \
    -PL Illumina \
    -O $BASE_DIR/out/1_ubam/${SAMPLE}.${LANEF}_ubam.bam \
    --TMP_DIR $BASE_DIR/tmp

## Execute manually after job completion
cd $BASE_DIR/out/1_ubam ; find *_ubam.bam > $BASE_DIR/fof/1_ubam.fof


### 2. Mark adapters

#!/bin/bash

#SBATCH --job-name=2_adap
#SBATCH --partition=carl.p
#SBATCH --array=1-126
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=0-03:00  # D-HH:MM
#SBATCH --output=log/2_adap/2_adap_%A_%a.out
#SBATCH --error=log/2_adap/2_adap_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8

BASE_DIR=/gss/work/haex1482/2_Other/egem/2_genotyping
UBAM_FILE=$BASE_DIR/fof/1_ubam.fof
UBAM=$(head $UBAM_FILE -n $SLURM_ARRAY_TASK_ID | tail -n 1)
LABEL=${UBAM%_ubam.bam}

gatk --java-options "-Xmx20G" \
   MarkIlluminaAdapters \
   -I $BASE_DIR/out/1_ubam/$UBAM \
   -O $BASE_DIR/out/2_adap/${LABEL}_adap.bam \
   -M $BASE_DIR/out/2_adap/${LABEL}_adap.txt \
   --TMP_DIR $BASE_DIR/tmp

## Execute manually after job completion
cd $BASE_DIR/out/2_adap ; find *_adap.bam > $BASE_DIR/fof/2_adap.fof


### 3. Map to reference genome

#!/bin/bash

#SBATCH --job-name=3_map
#SBATCH --partition=carl.p
#SBATCH --array=1-126
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=60G
#SBATCH --time=2-00:00  # D-HH:MM
#SBATCH --output=log/3_map/3_map_%A_%a.out
#SBATCH --error=log/3_map/3_map_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8
ml BWA/0.7.17-GCC-8.3.0

BASE_DIR=/gss/work/haex1482/2_Other/egem/2_genotyping
ADAP_FILE=$BASE_DIR/fof/2_adap.fof
ADAP=$(head $ADAP_FILE -n $SLURM_ARRAY_TASK_ID | tail -n 1)
LABEL=${ADAP%_adap.bam}

gatk --java-options "-Xmx60G" \
    SamToFastq \
    -I $BASE_DIR/out/2_adap/$ADAP \
    -FASTQ /dev/stdout \
    -INTERLEAVE true \
    -NON_PF true \
    --TMP_DIR $BASE_DIR/tmp |

bwa mem -M -t 8 \
    -p $BASE_DIR/ref/Hpue_LGM_01.fasta /dev/stdin |

gatk --java-options "-Xmx60G" \
    MergeBamAlignment \
    --VALIDATION_STRINGENCY SILENT \
    --EXPECTED_ORIENTATIONS FR \
    --ATTRIBUTES_TO_RETAIN X0 \
    -ALIGNED_BAM /dev/stdin \
    -UNMAPPED_BAM $BASE_DIR/out/1_ubam/${LABEL}_ubam.bam \
    -O $BASE_DIR/out/3_map/${LABEL}_map.bam \
    -R $BASE_DIR/ref/Hpue_LGM_01.fasta \
    -PAIRED_RUN true \
    --SORT_ORDER "unsorted" \
    --IS_BISULFITE_SEQUENCE false \
    --ALIGNED_READS_ONLY false \
    --CLIP_ADAPTERS false \
    --MAX_RECORDS_IN_RAM 2000000 \
    --ADD_MATE_CIGAR true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
    --ALIGNER_PROPER_PAIR_FLAGS true \
    --UNMAP_CONTAMINANT_READS true \
    --TMP_DIR $BASE_DIR/tmp

## Execute manually after job completion
cd $BASE_DIR/out/3_map ; find *_map.bam > $BASE_DIR/fof/3_map.fof


### 4. Mark duplicates and index

#!/bin/bash

#SBATCH --job-name=4_dedup
#SBATCH --partition=carl.p
#SBATCH --array=1-126
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=0-06:00  # D-HH:MM
#SBATCH --output=log/4_dedup/4_dedup_%A_%a.out
#SBATCH --error=log/4_dedup/4_dedup_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8

BASE_DIR=/gss/work/haex1482/2_Other/egem/2_genotyping
MAP_FILE=$BASE_DIR/fof/3_map.fof
MAP=$(head $MAP_FILE -n $SLURM_ARRAY_TASK_ID | tail -n 1)
LABEL=${MAP%_map.bam}

gatk --java-options "-Xmx40G" \
    SortSam \
    -I $BASE_DIR/out/3_map/$MAP \
    -O /dev/stdout \
    --SORT_ORDER "coordinate" \
    --CREATE_INDEX false \
    --CREATE_MD5_FILE false \
    --TMP_DIR $BASE_DIR/tmp |

gatk --java-options "-Xmx40G" \
    SetNmMdAndUqTags \
    -I /dev/stdin \
    -O $BASE_DIR/out/4_dedup/${LABEL}_int.bam \
    -R $BASE_DIR/ref/Hpue_LGM_01.fasta \
    --CREATE_INDEX false \
    --CREATE_MD5_FILE false \
    --TMP_DIR $BASE_DIR/tmp

gatk --java-options "-Xmx40G" \
    MarkDuplicates \
    -I $BASE_DIR/out/4_dedup/${LABEL}_int.bam \
    -O $BASE_DIR/out/4_dedup/${LABEL}_dedup.bam \
    -M $BASE_DIR/out/4_dedup/${LABEL}_dedup.txt \
    -MAX_FILE_HANDLES 1000  \
    --TMP_DIR $BASE_DIR/tmp

gatk --java-options "-Xmx40G" \
    BuildBamIndex \
    -I $BASE_DIR/out/4_dedup/${LABEL}_dedup.bam

## Execute manually after job completion
rm $BASE_DIR/out/4_dedup/*_int.bam
cd $BASE_DIR/out/4_dedup ; find *_dedup.bam > $BASE_DIR/fof/4_dedup.fof


### 5. Coverage analysis (skipped)

# a. calculate coverage
# b. create table and histogram
# c. remove low-coverage samples


### 6. Calculate genotype likelihoods *** update to automate multilane sample processing

#!/bin/bash

#SBATCH --job-name=6_like
#SBATCH --partition=carl.p
#SBATCH --array=1-122
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=3-00:00  # D-HH:MM
#SBATCH --output=log/6_like/6_like_%A_%a.out
#SBATCH --error=log/6_like/6_like_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8

BASE_DIR=/gss/work/haex1482/2_Other/egem/2_genotyping
DEDUP_FILE=$BASE_DIR/fof/4_dedup1.fof   # 1-lane samples only, multilane samples (18903nigpan, 18429puepan) were run manually
DEDUP=$(head $DEDUP_FILE -n $SLURM_ARRAY_TASK_ID | tail -n 1)
LABEL=${DEDUP%_dedup.bam}
SAMPLE=${LABEL%.?}

gatk --java-options "-Xmx40G" \
    HaplotypeCaller  \
    -I $BASE_DIR/out/4_dedup/${DEDUP} \
    -O $BASE_DIR/out/6_like/${SAMPLE}_g.vcf.gz \
    -R $BASE_DIR/ref/Hpue_LGM_01.fasta \
    -ERC GVCF

## Execute manually after job completion
find $BASE_DIR/out/6_like/*_g.vcf.gz > $BASE_DIR/fof/6_like.fof


### 7. Combine GVCF files into cohort GVCF

#!/bin/bash

#SBATCH --job-name=7_cohort
#SBATCH --partition=carl.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH --time=2-00:00  # D-HH:MM
#SBATCH --output=log/7_cohort/7_cohort_%j.out
#SBATCH --error=log/7_cohort/7_cohort_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8

BASE_DIR=/gss/work/haex1482/2_Other/egem/2_genotyping
LIKE_FILE=$BASE_DIR/fof/6_like.fof
INPUT=$(cat $LIKE_FILE | sed -e 's/^/-V /g')

gatk --java-options "-Xmx100G" \
    CombineGVCFs \
    $INPUT \
    -O $BASE_DIR/out/7_cohort/egem_cohort.g.vcf.gz \
    -R $BASE_DIR/ref/Hpue_LGM_01.fasta \


### 8. Genotype (all sites and SNPs only)

#!/bin/bash

#SBATCH --job-name=8_geno
#SBATCH --partition=carl.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=1-00:00  # D-HH:MM
#SBATCH --output=log/8_geno/8_geno_%j.out
#SBATCH --error=log/8_geno/8_geno_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8

BASE_DIR=/gss/work/haex1482/2_Other/egem/2_genotyping

gatk --java-options "-Xmx40G" \
    GenotypeGVCFs \
    -V $BASE_DIR/out/7_cohort/egem_cohort.g.vcf.gz \
    -O $BASE_DIR/out/8_geno/egem_inter_LGM.vcf.gz \
    -R $BASE_DIR/ref/Hpue_LGM_01.fasta \
    --include-non-variant-sites true

gatk --java-options "-Xmx40G" \
    SelectVariants \
    -V $BASE_DIR/out/8_geno/egem_inter_LGM.vcf.gz \
    -O $BASE_DIR/out/8_geno/egem_all_LGM.vcf.gz \
    -R $BASE_DIR/ref/Hpue_LGM_01.fasta \
    --select-type-to-exclude INDEL

gatk --java-options "-Xmx40G" \
    SelectVariants \
    -V $BASE_DIR/out/8_geno/egem_inter_LGM.vcf.gz \
    -O $BASE_DIR/out/8_geno/egem_snps_LGM.vcf.gz \
    -R $BASE_DIR/ref/Hpue_LGM_01.fasta \
    --select-type-to-include SNP

## Execute manually after job completion
rm $BASE_DIR/out/8_geno/egem_inter_*.vcf*


### Cleanup and stats

rm -rf tmp
mkdir scripts && mv *.js scripts/

bcftools stats egem_all_LGM.vcf.gz
# SN	0	number of samples:	124
# SN	0	number of records:	17051
# SN	0	number of no-ALTs:	15334
# SN	0	number of SNPs:	1715
# SN	0	number of MNPs:	0
# SN	0	number of indels:	19
# SN	0	number of others:	4
# SN	0	number of multiallelic sites:	115
# SN	0	number of multiallelic SNP sites:	95
