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

# Trimming

mkdir ./trim
find ./fastq/*.fastq.gz -maxdepth 1 -type f -exec bash -c \
    'bn=$(basename {}); cutadapt -a CTGTAGGCACCATCAAT -j 32 --minimum-length 20 -q 20 --trimmed-only {} -o trim/$bn' \; \
    > cutadapt.log

# Alignment

mkdir align
find ./fastq/*.fastq.gz -maxdepth 1 -type f -exec bash -c \
    's=$(basename {});bn=${s%.fastq.gz}; mkdir align/$bn' \;
find ./trim/*.fastq.gz -maxdepth 1 -type f -exec bash -c \
    's=$(basename {});bn=${s%.fastq.gz}; STAR --genomeDir ./genome/ --readFilesCommand zcat --quantMode TranscriptomeSAM --readFilesIn {} --outFileNamePrefix ./align/$bn/ --alignSJDBoverhangMin 1 --runThreadN 32 --outSAMtype BAM Unsorted' \;

# Filter uniquely mapped reads and perform sorting

find ./trim/*.fastq.gz -maxdepth 1 -type f -exec bash -c \
    's=$(basename {});bn=${s%.fastq.gz}; samtools view -b -q 255 ./align/$bn/Aligned.toTranscriptome.out.bam | samtools sort - > ./align/$bn.bam ' \;
find ./align/*.bam -maxdepth 1 -type f -exec samtools index {} \;
