import pysam
import os
import sys
bam_fn = sys.argv[1] #"source_data/alignments/METTL3_rna_Coots2017_m_r2.bam"

index_fn = os.path.abspath(bam_fn) + '.bai'
if not os.path.isfile(index_fn):
    pysam.index(bam_fn)

samfile = pysam.AlignmentFile(bam_fn, "rb", index_filename=index_fn)
transcripts = { read.reference_name for read in samfile }
for transcript in transcripts:
    for reads_at_pos in samfile.pileup(transcript):
        # read_names = [pileup_read.alignment.query_name for pileup_read in reads_at_pos.pileups]
        print(transcript, reads_at_pos.reference_pos, reads_at_pos.get_num_aligned() , sep='\t')
