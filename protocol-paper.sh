#!/usr/bin/env bash

# # Papolarity tools itself has no implicit dependencies, and no conventions about a structure of folders and file name
# #   and have a few conventions about file extensions (see README).
# # Thus all data necessary to run the command is specified in its arguments.
# # Feel free to change any part of protocol for your needs. Steps are named after paper sections.
# #
# # Nevertheless, this implementation of protocol relies on and maintains a certain file layout.
# # To start it, you should have genome annotation in GTF/GFF format (its name is specified in ANNOTATION variable)
# #   and a set of bam-files with alignment of reads in samples.
# # Alignments should be located at ./align/{SampleName}.bam
# # Sample names should be provided in SAMPLES, CONTROL and EXPERIMENTS variables.
# # If you don't want to compare samples, CONTROL and EXPERIMENTS aren't necessary
# #   (technically, steps before 3.3.4 don't need them).
# # If it's the case, specify SAMPLES explicitly and comment out steps in a section 3.3.*
# #
# # This pipeline uses several auxiliary tools: GNU parallel, csvtk and python package pasio.
# # You can find installation instructions in ./protocol-paper-installation.sh
# # To obtain data we used in a paper, follow ./protocol-paper-obtain-data.sh


CONTROL='ES_noHR_noCH_ribo';
EXPERIMENTS='ES_noHR_60sCH_ribo  ES_90sHR_60sCH_ribo  ES_120sHR_60sCH_ribo  ES_150sHR_60sCH_ribo  ES_180sHR_60sCH_ribo';
SAMPLES="$CONTROL $EXPERIMENTS";
# SAMPLES='ES_noHR_noCH_ribo  ES_noHR_60sCH_ribo  ES_90sHR_60sCH_ribo  ES_120sHR_60sCH_ribo  ES_150sHR_60sCH_ribo  ES_180sHR_60sCH_ribo';

ANNOTATION='./genome/gencode.vM23.basic.annotation.gtf'

DROP_5_FLANK=30
DROP_3_FLANK=30

# Genome is necessary only if you want to generate clipped-CDS sequences (step X.1.1)
GENOME='./genome/GRCm38.primary_assembly.genome.fa'

# 3.1. Common preprocessing
echo '=== 3.1. Common preprocessing ==='

# 3.1.1. Preprocessing transcripts annotation
echo '3.1.1. Preprocessing transcripts annotation'

mkdir -p ./genome;
papolarity cds_annotation \
    "$ANNOTATION" \
    --attr-filter transcript_type=protein_coding \
    --attr-filter gene_type=protein_coding \
    --output-file ./genome/cds_features.tsv

csvtk --tabs cut ./genome/cds_features.tsv \
                 --fields 'transcript_id,transcript_length,cds_length' \
                 --out-file ./genome/transcript_lengths.tsv

csvtk --tabs cut ./genome/cds_features.tsv \
                 --fields 'transcript_id,gene_id' \
                 --out-file ./genome/transcript2gene.tsv



# 3.1.2. Preparing coverage profiles
echo '3.1.2. Preparing coverage profiles'

mkdir -p ./coverage;
(
  for SAMPLE in $SAMPLES; do
    echo papolarity get_coverage "./align/${SAMPLE}.bam" \
                    --sort --dtype int \
                    --output-file "./coverage/${SAMPLE}.bedgraph.gz" ;
  done
) | parallel


# 3.1.3. Pooling coverage profiles
echo '3.1.3. Pooling coverage profiles'

papolarity pool_coverage ./coverage/*.bedgraph.gz --dtype int --output-file ./coverage/pooled.bedgraph.gz;


# 3.1.4. Clipping profiles withing coding segments
echo '3.1.4. Clipping profiles withing coding segments'

mkdir -p ./cds_coverage;
(
  for SAMPLE in $SAMPLES 'pooled'; do
    echo papolarity clip_cds \
                    ./genome/cds_features.tsv \
                    "./coverage/${SAMPLE}.bedgraph.gz" \
                    --drop-5-flank $DROP_5_FLANK  --drop-3-flank $DROP_3_FLANK \
                    --contig-naming original \
                    --output-file "./cds_coverage/${SAMPLE}.bedgraph.gz" ;
  done
) | parallel


#################################################

# 3.2. Polarity score estimation
echo '=== 3.2. Polarity score estimation ==='

# 3.2.1. Estimating polarity scores
echo '3.2.1. Estimating polarity scores'

mkdir -p ./coverage_features/raw;
(
  for SAMPLE in $SAMPLES 'pooled'; do
    echo papolarity coverage_features \
                    "./cds_coverage/${SAMPLE}.bedgraph.gz" \
                    --prefix "${SAMPLE}_" \
                    --output-file "./coverage_features/raw/${SAMPLE}.tsv"
  done
) | parallel


# 3.2.2. Filtering transcript lists
echo '3.2.2. Filtering transcript lists'

mkdir -p ./coverage_features/pooled
csvtk --tabs filter2 \
   "coverage_features/raw/pooled.tsv" \
   --filter '($pooled_mean_coverage >= 1) && ($pooled_coverage_q75 > 0)' \
   --out-file "coverage_features/pooled/pooled.filtered_1.tsv"


csvtk --tabs join \
    --fields transcript_id \
    "coverage_features/pooled/pooled.filtered_1.tsv" \
    "genome/transcript2gene.tsv" \
    "genome/transcript_lengths.tsv" \
    --out-file "coverage_features/pooled/pooled.filtered_1.with_gene_id.tsv"

papolarity choose_best \
    "coverage_features/pooled/pooled.filtered_1.with_gene_id.tsv" \
    pooled_mean_coverage \
    max \
    --group-by gene_id  --header \
    --output-file "coverage_features/pooled/pooled.filtered_2.tsv"

csvtk --tabs cut \
    "coverage_features/pooled/pooled.filtered_2.tsv" \
    --fields transcript_id,transcript_length,cds_length \
    --out-file ./transcripts_list.tsv


# 3.2.3. Finalizing the polarity score lists
echo '3.2.3. Finalizing the polarity score lists'

mkdir -p ./coverage_features/filtered;
for SAMPLE in $SAMPLES; do
    csvtk --tabs join \
        ./transcripts_list.tsv \
        "./coverage_features/raw/${SAMPLE}.tsv" \
        --out-file "./coverage_features/filtered/${SAMPLE}.tsv"
done

# 3.2.4. Polarity Z-score estimation
echo '3.2.4. Polarity Z-score estimation'

mkdir -p ./coverage_features/adjusted;
for SAMPLE in $SAMPLES; do
    papolarity adjust_features \
        "./coverage_features/filtered/${SAMPLE}.tsv" \
        --sort-field 'cds_length' \
        --fields "${SAMPLE}_polarity" \
        --mode z-score \
        --window 500 \
        --prefix 'zscore_' \
        --output-file "./coverage_features/adjusted/${SAMPLE}.tsv"
done

# 3.2.5. Plot per-sample polarity score distribution
echo '3.2.5. Plot per-sample polarity score distribution'

mkdir -p ./coverage_features/plot/;
for SAMPLE in $SAMPLES; do
    papolarity plot_distribution \
        "coverage_features/filtered/${SAMPLE}.tsv" \
        --fields "${SAMPLE}_polarity" \
        --no-legend \
        --title "${SAMPLE} polarity distribution" \
        --zero-line green \
        --xlim -1.0 1.0 \
        --output-file "coverage_features/plot/${SAMPLE}.png"
done

# 3.2.6. (supplementary step) Plot polarity score distribution for all samples on a single figure
echo '3.2.6. (supplementary step) Plot polarity score distribution for all samples on a single figure'

# Note: coverage_features/adjusted/all.tsv is not perfectly formatted - it has several identical columns.
# We use it for the only reason - to draw the plot.

SAMPLE_FILES_adjusted_features=$( echo $SAMPLES | xargs -n1 echo | xargs -n1 -I{} echo 'coverage_features/adjusted/{}.tsv' | tr '\n' ' ' )

csvtk --tabs join \
    ./transcripts_list.tsv \
    $SAMPLE_FILES_adjusted_features \
    --out-file coverage_features/adjusted/all.tsv;

SAMPLE_FIELDS_polarity=$( echo $SAMPLES | xargs -n1 echo | xargs -n1 -I{} echo '{}_polarity' | tr '\n' ' ' );

papolarity plot_distribution \
    "coverage_features/adjusted/all.tsv" \
    --fields $SAMPLE_FIELDS_polarity \
    --legend \
    --title "Polarity distributions" \
    --zero-line green \
    --xlim -1.0 1.0 \
    --output-file "coverage_features/plot/all.png"

#######################################################

# 3.3. Relative slope estimation
echo '=== 3.3. Relative slope estimation ==='

# 3.3.1. Segmentation of coverage profiles
echo '3.3.1. Segmentation of coverage profiles'

pasio ./coverage/pooled.bedgraph.gz --alpha 1 --beta 1 --output-file ./segmentation.bed.gz --output-mode bed

# 3.3.2. Clip segmentation
echo '3.3.2. Clip segmentation'

papolarity clip_cds \
    ./genome/cds_features.tsv  \
    ./segmentation.bed.gz  \
    --drop-5-flank $DROP_5_FLANK  --drop-3-flank $DROP_3_FLANK \
    --contig-naming original \
    --output-file ./cds_segmentation.bed.gz

# 3.3.3. (non-mandatory step) Generation of flattened coverage profiles.
echo '3.3.3. (non-mandatory step) Generation of flattened coverage profiles.'

mkdir -p ./cds_coverage_flattened;
(
  for SAMPLE in $SAMPLES 'pooled'; do
    echo papolarity flatten_coverage \
                    ./cds_segmentation.bed.gz \
                    "./cds_coverage/${SAMPLE}.bedgraph.gz" \
                    --only-matching \
                    --output-file "./cds_coverage_flattened/${SAMPLE}.bedgraph.gz";
  done
) | parallel

mkdir -p ./coverage_flattened;
(
  for SAMPLE in $SAMPLES 'pooled'; do
    echo papolarity flatten_coverage \
                    ./segmentation.bed.gz \
                    "./coverage/${SAMPLE}.bedgraph.gz" \
                    --only-matching \
                    --output-file "./coverage_flattened/${SAMPLE}.bedgraph.gz";
  done
) | parallel

# 3.3.4. Calculate slope for a pair of samples.
echo '3.3.4. Calculate slope for a pair of samples.'

mkdir -p ./comparison/raw;
(
  for EXPERIMENT in $EXPERIMENTS; do
      echo papolarity compare_coverage \
                      ./cds_segmentation.bed.gz \
                      "./cds_coverage/${CONTROL}.bedgraph.gz" \
                      "./cds_coverage/${EXPERIMENT}.bedgraph.gz" \
                      --segment-coverage-quantile 0.25 1 \
                      --prefix "${EXPERIMENT}_" \
                      --output-file "comparison/raw/${EXPERIMENT}.tsv"
  done
) | parallel

# 3.3.5. Finalizing profile comparison statistics
echo '3.3.5. Finalizing profile comparison statistics'

mkdir -p ./comparison/filtered;
for EXPERIMENT in $EXPERIMENTS; do
    csvtk --tabs join \
        ./transcripts_list.tsv \
        "comparison/raw/${EXPERIMENT}.tsv" \
        --out-file "./comparison/filtered/${EXPERIMENT}.tsv"
done

# 3.3.6. Adjust comparison statistics
echo '3.3.6. Adjust comparison statistics'

mkdir -p ./comparison/adjusted;
for EXPERIMENT in $EXPERIMENTS; do
    papolarity adjust_features \
        "comparison/filtered/${EXPERIMENT}.tsv" \
        --sort-field 'cds_length' \
        --fields "${EXPERIMENT}_slope" "${EXPERIMENT}_slopelog" "${EXPERIMENT}_l1_distance" "${EXPERIMENT}_polarity_diff" \
        --mode z-score \
        --window 500 \
        --prefix 'zscore_' \
        --output-file "./comparison/adjusted/${EXPERIMENT}.tsv"
done

# 3.3.7. Plot per-sample distributions of slope
echo '3.3.7. Plot per-sample distributions of slope'

# Note: plots are not proof-ready, typically you have to manually
# specify `--xlim min max` to set proper scaling for slope and slopelog plot

mkdir -p ./comparison/plot;
for EXPERIMENT in $EXPERIMENTS; do
    papolarity plot_distribution \
        "comparison/adjusted/${EXPERIMENT}.tsv" \
        --fields "${EXPERIMENT}_slopelog" \
        --no-legend \
        --title $'Distribution of linear regression slope\nfor normalized coverage log-ratios' \
        --zero-line green \
        --xlim -10 10 \
        --output-file "./comparison/plot/${EXPERIMENT}_slopelog.png"
done

for EXPERIMENT in $EXPERIMENTS; do
    papolarity plot_distribution \
        "comparison/adjusted/${EXPERIMENT}.tsv" \
        --fields "${EXPERIMENT}_slope" \
        --no-legend \
        --title $'Distribution of linear regression slope\nfor normalized coverage ratios' \
        --zero-line green \
        --xlim -100 100 \
        --output-file "./comparison/plot/${EXPERIMENT}_slope.png"
done

for EXPERIMENT in $EXPERIMENTS; do
    papolarity plot_distribution \
        "comparison/adjusted/${EXPERIMENT}.tsv" \
        --fields "${EXPERIMENT}_l1_distance" \
        --no-legend \
        --title $'Distribution of l1-distances\nbetween normalized coverage profiles' \
        --xlim 0 2 \
        --output-file "./comparison/plot/${EXPERIMENT}_l1_distance.png"
done


for EXPERIMENT in $EXPERIMENTS; do
    papolarity plot_distribution \
        "comparison/adjusted/${EXPERIMENT}.tsv" \
        --fields "${EXPERIMENT}_polarity_diff" \
        --no-legend \
        --title $'Distribution of polarity differences\nbetween experiment and control coverage profiles' \
        --xlim -2 2 \
        --output-file "./comparison/plot/${EXPERIMENT}_polarity_diff.png"
done

# 3.3.8. (supplementary step) Plot distributions of slopes distribution for all samples on a single figure
echo '3.3.8. (supplementary step) Plot distributions of slopes distribution for all samples on a single figure'

# Note: comparison/adjusted/all.tsv is not perfectly formatted - it has several identical columns.
# We use it for the only reason - to draw the plot.

SAMPLE_FILES_adjusted_comparison=$( echo $EXPERIMENTS | xargs -n1 echo | xargs -n1 -I{} echo 'comparison/adjusted/{}.tsv' | tr '\n' ' ' )

csvtk --tabs join \
    ./transcripts_list.tsv \
    $SAMPLE_FILES_adjusted_comparison \
    --out-file comparison/adjusted/all.tsv;

SAMPLE_FIELDS_slopelog=$( echo $EXPERIMENTS | xargs -n1 echo | xargs -n1 -I{} echo "{}_slopelog" | tr '\n' ' ' );

papolarity plot_distribution \
    "comparison/adjusted/all.tsv" \
    --fields $SAMPLE_FIELDS_slopelog \
    --labels $EXPERIMENTS \
    --legend \
    --title $'Distribution of linear regression slope\nfor normalized coverage log-ratios' \
    --xlim -10 10 \
    --zero-line green \
    --output-file "comparison/plot/all_slopelog.png";

SAMPLE_FIELDS_slope=$( echo $EXPERIMENTS | xargs -n1 echo | xargs -n1 -I{} echo "{}_slope" | tr '\n' ' ' );

papolarity plot_distribution \
    "comparison/adjusted/all.tsv" \
    --fields $SAMPLE_FIELDS_slope \
    --labels $EXPERIMENTS \
    --legend \
    --title $'Distribution of linear regression slope\nfor normalized coverage ratios' \
    --xlim -100 100 \
    --zero-line green \
    --output-file "comparison/plot/all_slope.png";


SAMPLE_FIELDS_l1=$( echo $EXPERIMENTS | xargs -n1 echo | xargs -n1 -I{} echo "{}_l1_distance" | tr '\n' ' ' );

papolarity plot_distribution \
    "comparison/adjusted/all.tsv" \
    --fields $SAMPLE_FIELDS_l1 \
    --labels $EXPERIMENTS \
    --legend \
    --title $'Distribution of l1-distances\nbetween normalized coverage profiles' \
    --xlim 0 2 \
    --output-file "comparison/plot/all_l1_distance.png";


SAMPLE_FIELDS_polarity_diff=$( echo $EXPERIMENTS | xargs -n1 echo | xargs -n1 -I{} echo "{}_polarity_diff" | tr '\n' ' ' );

papolarity plot_distribution \
    "comparison/adjusted/all.tsv" \
    --fields $SAMPLE_FIELDS_polarity_diff \
    --labels $EXPERIMENTS \
    --legend \
    --title $'Distribution of polarity differences\nbetween experiment and control coverage profiles' \
    --xlim -2 2 \
    --output-file "comparison/plot/all_polarity_diff.png";
################################################################################

# X.1.1. Supplementary step: generate CDS sequences
echo 'X.1.1. Supplementary step: generate CDS sequences'

papolarity cds_sequence "$ANNOTATION" \
                        "$GENOME" \
                        --attr-filter "transcript_type=protein_coding" \
                        --attr-filter "gene_type=protein_coding" \
                        --region-type cds \
                        --drop-5-flank $DROP_5_FLANK --drop-3-flank $DROP_3_FLANK \
                        --output-file cds_sequences.fa.gz
