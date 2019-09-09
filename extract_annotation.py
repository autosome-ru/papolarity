from annotation import *

relevant_attributes = {'gene_id', 'transcript_id', 'gene_type', 'transcript_type'}
annotation = Annotation.load('gencode.vM22.basic.annotation.gtf.gz', relevant_attributes=relevant_attributes, coding_only=True)
with open('cds_features.tsv', 'w') as fw:
    print('gene_id', 'transcript_id', 'transcript_length', 'cds_length', 'cds_start', 'cds_stop', sep='\t', file=fw)
    for transcript_info in coding_transcript_infos(annotation):
        gene_id        = transcript_info.gene_id
        transcript_id  = transcript_info.transcript_id
        len_transcript = transcript_info.transcript_length
        len_cds        = transcript_info.cds_length
        cds_start      = transcript_info.cds_start
        cds_stop       = transcript_info.cds_stop
        print(gene_id, transcript_id, len_transcript, len_cds, cds_start, cds_stop, sep='\t', file=fw)
