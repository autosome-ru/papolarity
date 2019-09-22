import pybedtools
import itertools
import numpy as np
import sys
from annotation import CodingTranscriptInfo, load_transcript_cds_info

def transcript_coverages_in_file(alignment_filename, cds_info_by_transcript, sort_transcripts=False):
    bam = pybedtools.BedTool(alignment_filename)
    bg = bam.genome_coverage(bg=True)
    if sort_transcripts:
        bg = bg.sort()
    for (transcript_id, bedgraph_intervals_iter) in itertools.groupby(bg, lambda interval: interval.chrom):
        if transcript_id not in cds_info_by_transcript:
            # We don't know length of such transcripts
            continue
        transcript_info = cds_info_by_transcript[transcript_id]
        profile = profile_by_bedgraph(bedgraph_intervals_iter, transcript_info.transcript_length)
        yield (transcript_info, profile)

def profile_by_bedgraph(bedgraph_intervals_iter, transcript_length):
    # Check: is it reasonable to use real-values or integer-valued coverages?
    profile = np.zeros(transcript_length, dtype=np.int)
    for interval in bedgraph_intervals_iter:
        # Note! pybedtools use 0-based coordinates (when in integer representation).
        # See https://daler.github.io/pybedtools/3-brief-examples.html and https://daler.github.io/pybedtools/intervals.html#zero-based-coords
        coverage = int(interval[3])
        profile[interval.start : interval.stop] = coverage
    return profile

cds_annotation_fn = sys.argv[1] # 'gencode.vM22.cds_features.tsv'
alignment_filename = sys.argv[2] # 'METTL3_ribo_Coots2017_m_r2.bam'
cds_info_by_transcript = load_transcript_cds_info(cds_annotation_fn)
header = '\t'.join([CodingTranscriptInfo.header(), 'cds_total_coverage'])
print(header)
for (transcript_info, profile) in transcript_coverages_in_file(alignment_filename, cds_info_by_transcript):
    cds_profile = transcript_info.cds_profile(profile)
    print(transcript_info, ','.join(map(str, cds_profile)), sep='\t')