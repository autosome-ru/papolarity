import pybedtools
import itertools
import numpy as np
import sys
from annotation import CodingTranscriptInfo

def load_transcript_cds_info(cds_annotation_filename):
    transcript_infos = CodingTranscriptInfo.each_from_file(cds_annotation_filename)
    return {tr_info.transcript_id: tr_info  for tr_info in transcript_infos}

def transcript_coverages_in_file(alignment_filename, cds_info_by_transcript):
    bam = pybedtools.BedTool(alignment_filename)
    bg = bam.genome_coverage(bg=True)
    for (transcript_id, iterator) in itertools.groupby(bg, lambda interval: interval.chrom):
        if transcript_id not in cds_info_by_transcript:
            continue
        transcript_info = cds_info_by_transcript[transcript_id]
        profile = np.zeros(transcript_info.transcript_length, dtype=np.int)
        for interval in iterator:
            # Note! pybedtools returns coordinates of `interval` in [0-closed, 1-open) format
            coverage = int(interval[3])
            profile[interval.start : interval.stop] = coverage
        yield (transcript_info, profile)


cds_annotation_fn = sys.argv[1] # 'gencode.vM22.cds_features.tsv'
alignment_filename = sys.argv[2] # 'METTL3_ribo_Coots2017_m_r2.bam'
cds_info_by_transcript = load_transcript_cds_info(cds_annotation_fn)
header = '\t'.join([CodingTranscriptInfo.header(), 'cds_total_coverage'])
print(header)
for (transcript_info, profile) in transcript_coverages_in_file(alignment_filename, cds_info_by_transcript):
    cds_profile = transcript_info.cds_profile(profile)
    print(transcript_info, cds_profile.sum(), sep='\t')
