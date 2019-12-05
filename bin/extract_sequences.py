from annotation import Annotation

gtf_annotation_fn = 'source_data/annotation/gencode.vM22.basic.annotation.gtf.gz'
assembly_fn = '/home/ilya/iogen/data/genome/mm10.fa'

annotation = Annotation.load(gtf_annotation_fn, coding_only=True, relevant_attributes=set())
transcript_ids_list = annotation.transcript_by_id.keys()

# annotation = Annotation.load(gtf_annotation_fn, coding_only=True, relevant_attributes={'level', 'transcript_support_level'})
# transcript_ids_list = annotation.transcript_by_id.keys()
# transcript_ids_list = filter(lambda tr_id: annotation.transcript_by_id[tr_id].attributes['level'] in [1, 2], transcript_ids_list)
# transcript_ids_list = filter(lambda tr_id: annotation.transcript_by_id[tr_id].attributes.get('transcript_support_level') in ['1', '2'], transcript_ids_list)

for transcript_id, seq in annotation.transcript_sequences(transcript_ids_list, assembly_fn, feature_type='cds_with_stop'):
    print(f'>{transcript_id}')
    print(seq)
