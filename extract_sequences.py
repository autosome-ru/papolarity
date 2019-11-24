from annotation import Annotation
from coding_transcript_info import CodingTranscriptInfo

gtf_annotation_fn = 'source_data/annotation/gencode.vM22.basic.annotation.gtf.gz'
assembly_fn = '/home/ilya/iogen/data/genome/mm10.fa'

annotation = Annotation.load(gtf_annotation_fn, coding_only=True)
transcript_ids_list = annotation.transcript_by_id.keys()

for transcript_id, seq in annotation.transcript_sequences(transcript_ids_list, assembly_fn, feature_type='cds'):
    print(f'>{transcript_id}')
    print(seq)
