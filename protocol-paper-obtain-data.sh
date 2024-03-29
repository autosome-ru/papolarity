#!/usr/bin/env bash

# Download genome and annotation

mkdir ./genome
cd ./genome
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.basic.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
find . -type f -name '*.gz' -exec bash -c 'gzip -d {}' \;


# index the genome and annotation with STAR aligner

STAR --runMode genomeGenerate --genomeDir . --genomeFastaFiles ./GRCm38.primary_assembly.genome.fa --sjdbGTFfile ./gencode.vM23.basic.annotation.gtf --runThreadN 32
cd ..

# Obtain raw data

SRR='SRR315616 SRR315617 SRR315618 SRR315619 SRR315601 SRR315602 SRR315612 SRR315613 SRR315614 SRR315615 SRR315604 SRR315605 SRR315606 SRR315607 SRR315608 SRR315609 SRR315610 SRR315611';
prefetch $SRR;
mkdir raw;
for SRR_FILE in $SRR; do fastq-dump -O ./raw --gzip $SRR_FILE --split-files; done;
mkdir fastq;
cat ./raw/SRR315616_1.fastq.gz ./raw/SRR315617_1.fastq.gz ./raw/SRR315618_1.fastq.gz ./raw/SRR315619_1.fastq.gz > ./fastq/ES_noHR_noCH_ribo.fastq.gz;
cat ./raw/SRR315601_1.fastq.gz ./raw/SRR315602_1.fastq.gz > ./fastq/ES_noHR_60sCH_ribo.fastq.gz;
cat ./raw/SRR315612_1.fastq.gz ./raw/SRR315613_1.fastq.gz ./raw/SRR315614_1.fastq.gz ./raw/SRR315615_1.fastq.gz > ./fastq/ES_90sHR_60sCH_ribo.fastq.gz;
cat ./raw/SRR315604_1.fastq.gz ./raw/SRR315605_1.fastq.gz ./raw/SRR315606_1.fastq.gz > ./fastq/ES_120sHR_60sCH_ribo.fastq.gz;
cat ./raw/SRR315607_1.fastq.gz ./raw/SRR315608_1.fastq.gz ./raw/SRR315609_1.fastq.gz > ./fastq/ES_150sHR_60sCH_ribo.fastq.gz;
cat ./raw/SRR315610_1.fastq.gz ./raw/SRR315611_1.fastq.gz > ./fastq/ES_180sHR_60sCH_ribo.fastq.gz;

# Trimming

mkdir ./trim
find ./fastq/*.fastq.gz -maxdepth 1 -type f -exec bash -c \
    'bn=$(basename {}); cutadapt -a CTGTAGGCACCATCAAT -j 32 --minimum-length 20 -q 20 --trimmed-only {} -o trim/$bn' \; \
    > cutadapt.log

# Alignment
# !!!!! ATTENTION !!!!!
# Here one should make a transcriptomic alignment, NOT a genomic alignment
# !!!!!           !!!!!
#
# NOTE! Genome-level multi-mapping reads are excluded at the alignment step.
# It might be problematic to correctly filter out multi-mapping reads later from the transcriptomic
# alignment as there are often many overlapping transcripts of a single gene

mkdir align
find ./fastq/*.fastq.gz -maxdepth 1 -type f -exec bash -c \
    's=$(basename {});bn=${s%.fastq.gz}; mkdir align/$bn' \;
find ./trim/*.fastq.gz -maxdepth 1 -type f -exec bash -c \
    's=$(basename {});bn=${s%.fastq.gz}; STAR --genomeDir ./genome/ --readFilesCommand zcat --quantMode TranscriptomeSAM --readFilesIn {} --outFileNamePrefix ./align/$bn/ --alignSJDBoverhangMin 1 --outFilterMultimapNmax 1 --runThreadN 32 --outSAMtype BAM Unsorted' \;

# Sort and index the read alignment.
# !NOTE! The earlier version of the protocol suggested filtering multi-mapping reads with samtools.
# We no longer recommend such filtering as it will exclude the reads mapped to alternative
# overlapping transcripts of a single gene and significantly reduce the effective read coverage.
# Instead, when considering only uniquely-mapped reads, we suggest to perform read mapping allowing
# only single mapping locations, see STAR parameter `--outFilterMultimapNmax 1` above.

find ./trim/*.fastq.gz -maxdepth 1 -type f -exec bash -c \
    's=$(basename {});bn=${s%.fastq.gz}; samtools sort ./align/$bn/Aligned.toTranscriptome.out.bam > ./align/$bn.bam ' \;
find ./align/*.bam -maxdepth 1 -type f -exec samtools index {} \;
