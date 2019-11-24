from annotation import Annotation
from coding_transcript_info import CodingTranscriptInfo
import sys

gtf_annotation_fn = sys.argv[1] # 'gencode.vM22.basic.annotation.gtf.gz'

annotation = Annotation.load(gtf_annotation_fn, relevant_attributes=set(), coding_only=True)

print(CodingTranscriptInfo.header())
for transcript_id in annotation.transcript_by_id:
    transcript_info = annotation.coding_transcript_info(transcript_id)
    print(transcript_info)
