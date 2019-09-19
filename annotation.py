from collections import defaultdict
from collections import namedtuple
from gtf_parser import *
from coding_transcript_info import CodingTranscriptInfo

class Annotation:
    def __init__(self):
        self.gene_by_id = {}
        self.transcript_by_id = {}
        self.geneId_by_transcript = {}
        self.transcripts_by_gene = defaultdict(list)
        self.parts_by_transcript = defaultdict(list)

    @classmethod
    def load(cls, filename, relevant_attributes=None, coding_only=False):
        annotation = cls()

        if relevant_attributes:
            attributes_filter = lambda attributes: {k: v  for (k,v) in attributes.items()  if k in relevant_attributes}
        else:
            attributes_filter = lambda x: x

        records = parse_gtf(filename, attributes_filter=attributes_filter)
        if coding_only:
            records = filter_coding(records)

        for rec in records:
            annotation.push_record(rec)
        return annotation

    def push_record(self, rec):
        gene_id = rec.attributes['gene_id']
        if rec.type == 'gene':
            self.gene_by_id[gene_id] = rec
        else:
            transcript_id = rec.attributes['transcript_id']
            if rec.type == 'transcript':
                self.transcript_by_id[transcript_id] = rec
                self.geneId_by_transcript[transcript_id] = gene_id
                self.transcripts_by_gene[gene_id].append(rec)
            else:
                self.parts_by_transcript[transcript_id].append(rec)

    def transcript_exons(self, transcript_id):
        parts = self.parts_by_transcript[transcript_id]
        return [part  for part in parts  if part.type == 'exon']

    def transcript_cds(self, transcript_id):
        parts = self.parts_by_transcript[transcript_id]
        return [part  for part in parts  if part.type == 'CDS']

    def transcript_strand(self, transcript_id):
        segment_strands = { segment.strand  for segment in  self.parts_by_transcript[transcript_id] }
        if len(segment_strands) != 1:
            raise Exception('Different strands')
        return segment_strands.pop()

    def coding_transcript_info(self, transcript_id):
        """Return coding transcript info. This method fails when called for non-coding transcript"""
        gene_id = self.geneId_by_transcript[transcript_id]
        exons = self.transcript_exons(transcript_id)
        cds_segments = self.transcript_cds(transcript_id)
        strand = self.transcript_strand(transcript_id)

        transcript_length = sum(exon.length for exon in exons)
        len_cds_in_transcript = sum(cds.length for cds in cds_segments)
        if strand == '+':
            genomic_cds_start = min([rec.start for rec in cds_segments])
            len_exons_before_cds = sum( exon.length  for exon in exons  if exon.in_upstream_of(genomic_cds_start) )
            exon_with_cds = take_the_only([exon for exon in exons  if exon.contain_position(genomic_cds_start)])
            cds_offset = genomic_cds_start - exon_with_cds.start  # offset of CDS relative to embracing exon
        elif strand == '-':
            genomic_cds_start = max([rec.end for rec in cds_segments])
            len_exons_before_cds = sum( exon.length  for exon in exons  if exon.in_upstream_of(genomic_cds_start) )
            exon_with_cds = take_the_only([exon for exon in exons  if exon.contain_position(genomic_cds_start)])
            cds_offset = exon_with_cds.end - genomic_cds_start  # offset of CDS relative to embracing exon

        cds_start = cds_offset + len_exons_before_cds
        cds_stop = cds_start + len_cds_in_transcript - 1

        return CodingTranscriptInfo(gene_id, transcript_id, transcript_length, cds_start, cds_stop)


# take only elements related to coding transcripts of coding genes
def filter_coding(iter):
    for rec in iter:
        if rec.is_coding():
            yield rec

def take_the_only(arr):
    if len(arr) > 1:
        raise Exception('Several elements when the only one is expected')
    if len(arr) == 0:
        raise Exception('No elements when one is expected')
    return arr[0]
