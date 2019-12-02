import sys
from coding_transcript_info import CodingTranscriptInfo
from coverage_profile import transcript_coverages_in_file

cds_annotation_fn = sys.argv[1] # 'gencode.vM22.cds_features.tsv'
alignment_filename = sys.argv[2] # 'METTL3_ribo_Coots2017_m_r2.bam'
cds_info_by_transcript = CodingTranscriptInfo.load_transcript_cds_info(cds_annotation_fn)
header = '\t'.join([CodingTranscriptInfo.header(), 'cds_total_coverage'])
print(header)
for (transcript_info, profile) in transcript_coverages_in_file(alignment_filename, cds_info_by_transcript):
    cds_profile = transcript_info.cds_profile(profile)
    print(transcript_info, ','.join(map(str, cds_profile)), sep='\t')
