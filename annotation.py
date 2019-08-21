from collections import defaultdict
from gtf_parser import *

def take_the_only(arr):
    if len(arr) > 1:
        raise Exception('Several elements when the only one is expected')
    if len(arr) == 0:
        raise Exception('No elements when one is expected')
    return arr[0]

def load_annotation_parts(filename):
    annotation_parts = {
        'genes': {},
        'geneId_by_transcript': {},
        'transcripts_by_gene': defaultdict(list),
        'parts_by_transcript': defaultdict(list),
    }
    for rec in parse_gtf(filename):
    # for rec in parse_gtf('annotation_large_sample.gtf'):
        if rec.attributes['gene_type'] != 'protein_coding':
            continue

        gene_id = rec.attributes['gene_id']
        if rec.type == 'gene':
            annotation_parts['genes'][gene_id] = rec

        else:
            if rec.attributes['transcript_type'] != 'protein_coding':
                continue

            transcript_id = rec.attributes['transcript_id']
            if rec.type == 'transcript':
                annotation_parts['geneId_by_transcript'][transcript_id] = gene_id
                annotation_parts['transcripts_by_gene'][gene_id].append(rec)
            else:
                annotation_parts['parts_by_transcript'][transcript_id].append(rec)
    return annotation_parts

def cds_annotation(annotation_parts):
    annotation = {
        'geneId_by_transcript': annotation_parts['geneId_by_transcript'],
        'len_exons_by_transcript': {},
        'len_cds_by_transcript':   {},
        'cds_start_by_transcript': {},
        'cds_stop_by_transcript':  {},
    }

    for gene_id, gene in annotation_parts['genes'].items():
        for transcript in annotation_parts['transcripts_by_gene'][gene_id]:
            transcript_id = transcript.attributes['transcript_id']
            parts = annotation_parts['parts_by_transcript'][transcript_id]

            exons = [part   for part in parts   if part.type == 'exon']
            cds_segments = [part for part in parts if part.type == 'CDS']
            if len( { segment.strand for segment in (exons + cds_segments) } ) != 1:
                raise Exception('Different strands')
            strand = exons[0].strand

            len_transcript = sum((exon.end - exon.start + 1) for exon in exons)
            len_cds_in_transcript = sum((cds.end - cds.start + 1) for cds in cds_segments)

            if strand == '+':
                # genomic coordinates
                cds_start = min([rec.start for rec in cds_segments])
                len_exons_before_cds = sum( exon.length() for exon in exons  if exon.end < cds_start  )
                exon_with_cds = take_the_only([exon for exon in exons  if exon.start <= cds_start <= exon.end])
                len_exon_with_cds = cds_start - exon_with_cds.start + 1 
                cds_start_transcript = len_exon_with_cds + len_exons_before_cds
                cds_stop_transcript = cds_start_transcript + len_cds_in_transcript - 1
            elif strand == '-':
                cds_start = max([rec.end for rec in cds_segments])
                len_exons_before_cds = sum( exon.length() for exon in exons  if cds_start < exon.start  )
                exon_with_cds = take_the_only([exon for exon in exons  if exon.start <= cds_start <= exon.end])
                len_exon_with_cds = exon_with_cds.end - cds_start + 1 
                cds_start_transcript = len_exon_with_cds + len_exons_before_cds
                cds_stop_transcript = cds_start_transcript + len_cds_in_transcript - 1
            else:
                raise Exception('Unknown strand')

            annotation['len_exons_by_transcript'][transcript_id] = len_transcript
            annotation['len_cds_by_transcript'][transcript_id] = len_cds_in_transcript
            annotation['cds_start_by_transcript'][transcript_id] = cds_start_transcript
            annotation['cds_stop_by_transcript'][transcript_id] = cds_stop_transcript
    return annotation
