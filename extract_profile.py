import pybedtools
import itertools
import numpy as np

info_by_transcript = {}
with open('cds_features.tsv') as f:
    header = f.readline()
    for line in f:
        row = line.rstrip("\n").split("\t")
        gene_id, transcript_id, len_transcript, len_cds, cds_start, cds_stop = row
        info_by_transcript[transcript_id] = {
            'gene_id': gene_id,
            'len_transcript': int(len_transcript),
            'len_cds': int(len_cds),
            'cds_start': int(cds_start),
            'cds_stop': int(cds_stop),
        }

bam = pybedtools.BedTool('source_data/alignments/METTL3_rna_Coots2017_m_r2.bam')
bg = bam.genome_coverage(bg=True)
for (transcript_id, iterator) in itertools.groupby(bg, lambda interval: interval.chrom):
    if transcript_id not in info_by_transcript:
        continue
    transcript_info = info_by_transcript[transcript_id]
    transcript_len = transcript_info['len_transcript']
    profile = np.zeros(transcript_len)
    for interval in iterator:
        coverage = int(interval[3])
        # note that coordinates are 1-based!
        profile[interval.start - 1 : interval.stop] = coverage
    cds_profile = profile[transcript_info['cds_start'] - 1: transcript_info['cds_stop']]
    # if len(cds_profile) % 3 != 0:
    #     print(transcript_id, len(cds_profile), transcript_info)
