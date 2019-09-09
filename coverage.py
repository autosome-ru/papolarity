import pysam
import os
import sys
from collections import namedtuple

TranscriptAnnotation = namedtuple('TranscriptAnnotation', ['gene_id', 'transcript_id', 'len_transcript', 'len_cds', 'cds_start', 'cds_stop'])
def load_cds_features(filename):
    result = []
    with open(filename) as f:
        for line in f:
            gene_id, transcript_id, len_transcript, len_cds, cds_start, cds_stop = f.rstrip("\n").split("\t")
            result.append(TranscriptAnnotation(gene_id, transcript_id, int(len_transcript), int(len_cds), int(cds_start), int(cds_stop)))
    return result

transcripts = load_cds_features('cds_features.tsv')


bam_fn = sys.argv[1] #"source_data/alignments/METTL3_rna_Coots2017_m_r2.bam"

index_fn = os.path.abspath(bam_fn) + '.bai'
if not os.path.isfile(index_fn):
    pysam.index(bam_fn)

samfile = pysam.AlignmentFile(bam_fn, "rb")
transcripts = { read.reference_name for read in samfile }
for transcript in transcripts:
    for reads_at_pos in list(samfile.pileup(transcript)):
        # read_names = [pileup_read.alignment.query_name for pileup_read in reads_at_pos.pileups]
        print(transcript, reads_at_pos.reference_pos, reads_at_pos.get_num_aligned() , sep='\t')

# profiles = []

# for read in samfile.fetch(until_eof=True):
