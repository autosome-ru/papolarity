from collections import defaultdict
from itertools import groupby
from gtf_parser import *

def get_transcript_id(rec):
    return rec.attributes['transcript_id']

genes = {}
transcripts_by_gene = defaultdict(list)
parts_by_transcript = defaultdict(list)
for rec in parse_gtf('gencode.vM22.basic.annotation.gtf.gz'):
#for rec in parse_gtf('annotation_sample.gtf'):
    if rec.attributes['gene_type'] != 'protein_coding':
        continue

    if rec.type == 'gene':
        gene_id = rec.attributes['gene_id']
        genes[gene_id] = rec

    else:
        if rec.attributes['transcript_type'] != 'protein_coding':
            continue

        if rec.type == 'transcript':
            gene_id = rec.attributes['gene_id']
            transcripts_by_gene[gene_id].append(rec)
        else:
            transcript_id = rec.attributes['transcript_id']
            parts_by_transcript[transcript_id].append(rec)

for gene_id, gene in genes.items():
    for transcript in transcripts_by_gene[gene_id]:
        transcript_id = transcript.attributes['transcript_id']
        parts = parts_by_transcript[transcript_id]

        exons = [part for part in parts if part.type == 'exon']
        cds_segments = [part for part in parts if part.type == 'CDS']
        if len( { segment.strand for segment in (exons + cds_segments) } ) != 1:
            raise Exception('Different strands')
        strand = exons[0].strand

        if strand == '+':
            # genomic coordinates
            cds_start = min([rec.start for rec in cds_segments])
            len_exons_before_cds = sum(  (exon.end - exon.start + 1) for exon in exons  if exon.end < cds_start  )
            exons_with_cds = [exon for exon in exons  if exon.start <= cds_start <= exon.end]
            if len(exons_with_cds) > 1:
                raise Exception('Several exons contain CDS start simultaneously')
            if len(exons_with_cds) == 0:
                raise Exception('No exons contain CDS start')
            exon_with_cds = exons_with_cds[0]
            len_exon_with_cds = cds_start - exon_with_cds.start + 1 
            cds_start_transcript = len_exon_with_cds + len_exons_before_cds
            print(transcript_id, cds_start_transcript, sep='\t')

            # ToDo: посчитать cds_stop_transcript
        # ToDo: посчитать start и stop для отрицательной нитки
