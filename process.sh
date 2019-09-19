python3 extract_annotation.py gencode.vM22.basic.annotation.gtf.gz > gencode.vM22.cds_features.tsv
python3 extract_profile.py gencode.vM22.cds_features.tsv source_data/alignments/METTL3_rna_Coots2017_m_r2.bam 

# samtools view source_data/alignments/METTL3_rna_Coots2017_m_r2.bam > METTL3.tsv
# perl gtf.pl METTL3.tsv > METTL3_coverage.tsv
