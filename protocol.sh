#!/usr/bin/env bash

SAMPLE_BNS='ES_noHR_noCH_ribo  ES_noHR_60sCH_ribo  ES_90sHR_60sCH_ribo  ES_120sHR_60sCH_ribo  ES_150sHR_60sCH_ribo  ES_180sHR_60sCH_ribo';

# 3.1. Common preprocessing

# 3.1.1. Preprocessing transcripts annotation

papolarity extract_cds_annotation \
    ./genome/gencode.vM23.basic.annotation.gtf \
    --attr-filter transcript_type=protein_coding \
    --attr-filter gene_type=protein_coding \
    --output-file ./genome/gencode.vM23.cds_features.tsv


csvtk --tabs cut genome/gencode.vM23.cds_features.tsv \
                 --fields 'transcript_id,transcript_length,cds_length' \
                 --out-file genome/transcript_lengths.tsv

csvtk --tabs cut genome/gencode.vM23.cds_features.tsv \
                 --fields 'transcript_id,gene_id' \
                 --out-file genome/transcript2gene.tsv



# 3.1.2. Preparing coverage profiles

mkdir -p ./coverage/;
find ./align/ -maxdepth 1 -xtype f -name '*.bam' | xargs -n1 basename -s .bam | xargs -n1 -I{} echo 'papolarity get_coverage ./align/{}.bam --sort --dtype int --output-file ./coverage/{}.bedgraph.gz' | parallel;


# 3.1.3. Pooling coverage profiles

papolarity pool_coverages ./coverage/*.bedgraph.gz --dtype int --output-file ./coverage/pooled.bedgraph.gz;


# 3.1.4. Clipping profiles withing coding segments

mkdir -p ./cds_coverage;
find ./coverage/ -maxdepth 1 -xtype f -name '*.bedgraph.gz' \
  | xargs -n1 basename \
  | xargs -n1 -I{} echo \
    'papolarity clip_cds ./genome/gencode.vM23.cds_features.tsv  ./coverage/{} --drop-5-flank 15  --drop-3-flank 15  --contig-naming original --output-file ./cds_coverage/{}' \
  | parallel


#################################################

# 3.2. Polarity score estimation

# 3.2.1. Estimating polarity scores

mkdir -p ./coverage_features/raw;
(
for SAMPLE_BN in $SAMPLE_BNS 'pooled'; do
    echo papolarity coverage_features \
                    "./cds_coverage/${SAMPLE_BN}.bedgraph.gz" \
                    --prefix "${SAMPLE_BN}_" \
                    --output-file "./coverage_features/raw/${SAMPLE_BN}.tsv"
  done
) | parallel


# 3.2.2. Filtering transcript lists
mkdir -p ./coverage_features/pooled
csvtk --tabs filter2 \
   "coverage_features/raw/pooled.tsv" \
   --filter '$pooled_mean_coverage >= 5' \
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

mkdir -p ./coverage_features/filtered;
for SAMPLE_BN in $SAMPLE_BNS; do
    csvtk --tabs join \
        ./transcripts_list.tsv \
        "./coverage_features/raw/${SAMPLE_BN}.tsv" \
        --out-file "./coverage_features/filtered/${SAMPLE_BN}.tsv"
done

# 3.2.4. Polarity Z-score estimation

mkdir -p ./coverage_features/adjusted;
for SAMPLE_BN in $SAMPLE_BNS; do
    papolarity adjust \
        "./coverage_features/filtered/${SAMPLE_BN}.tsv" \
        --sort-field 'cds_length' \
        --fields "${SAMPLE_BN}_polarity" \
        --mode z-score \
        --window 500 \
        --prefix 'zscore_' \
        --output-file "./coverage_features/adjusted/${SAMPLE_BN}.tsv"
done

# 3.2.5. Plot per-sample polarity score distribution

mkdir -p ./coverage_features/plot/;
for SAMPLE_BN in $SAMPLE_BNS; do
    papolarity plot_distribution \
        "coverage_features/filtered/${SAMPLE_BN}.tsv" \
        --fields "${SAMPLE_BN}_polarity" \
        --no-legend \
        --title "${SAMPLE_BN} polarity distribution" \
        --zero-line green \
        --xlim -1.0 1.0 \
        --output-file "coverage_features/plot/${SAMPLE_BN}.png"
done


#######################################################


SAMPLE_FILES=$( echo $SAMPLE_BNS | xargs -n1 echo | xargs -n1 -I{} echo 'coverage_features/filtered/{}.tsv' | tr '\n' ' ' )

csvtk --tabs join \
    ./transcripts_list.tsv \
    $SAMPLE_FILES \
    --out-file coverage_features/filtered/all.tsv;

SAMPLE_FIELDS=$( echo $SAMPLE_BNS | xargs -n1 echo | xargs -n1 -I{} echo '{}_polarity' | tr '\n' ' ' );

papolarity plot_distribution \
    "coverage_features/filtered/all.tsv" \
    --fields $SAMPLE_FIELDS \
    --legend \
    --title "Polarity distributions" \
    --zero-line green \
    --xlim -1.0 1.0 \
    --output-file "coverage_features/plot/all.png"


#######################################################


