from collections import defaultdict
from collections import namedtuple
from gtf_parser import *
from coding_transcript_info import CodingTranscriptInfo
import itertools
import pybedtools

class Annotation:
    def __init__(self):
        self.gene_by_id = {}
        self.transcript_by_id = {}
        self.geneId_by_transcript = {}
        self.transcripts_by_gene = defaultdict(list)
        self.parts_by_transcript = defaultdict(list)

    @classmethod
    def load(cls, filename, relevant_attributes=None, coding_only=False, attr_mapping=None):
        annotation = cls()

        if relevant_attributes is not None:
            # that's the minimum list of attributes for a library to properly work
            necessary_attributes = {'gene_id', 'transcript_id', 'gene_type', 'transcript_type', 'exon_number'}
            relevant_attributes = relevant_attributes.union(necessary_attributes)
            # We should preserve attributes in order to rename them later
            if attr_mapping:
                relevant_attributes = relevant_attributes.union(attr_mapping.keys())

        records = parse_gtf(filename, relevant_attributes=relevant_attributes)
        if attr_mapping:
            records = map(lambda rec: rec.attributes_renamed(attr_mapping), records)

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
            genomic_cds_start = max([rec.stop - 1 for rec in cds_segments])
            len_exons_before_cds = sum( exon.length  for exon in exons  if exon.in_upstream_of(genomic_cds_start) )
            exon_with_cds = take_the_only([exon for exon in exons  if exon.contain_position(genomic_cds_start)])
            cds_offset = exon_with_cds.stop - 1 - genomic_cds_start  # offset of CDS relative to embracing exon

        cds_start = cds_offset + len_exons_before_cds
        cds_stop = cds_start + len_cds_in_transcript

        return CodingTranscriptInfo(gene_id, transcript_id, transcript_length, cds_start, cds_stop)

    def ordered_segments_by_type(self, transcript_id, feature_type):
        by_exon_number = lambda segment: segment.attributes['exon_number']
        if feature_type == 'exons':
            return sorted(self.transcript_exons(transcript_id), key=by_exon_number)
        elif feature_type == 'cds':
            return sorted(self.transcript_cds(transcript_id), key=by_exon_number)
        elif feature_type == 'full':
            # full, unspliced transcript
            return [self.transcript_by_id(transcript_id)]


    def segments_as_bedtool_intervals(self, segments, name='.'):
        yield from (pybedtools.Interval(s.contig, s.start, s.stop, strand=s.strand, name=name) for s in segments)

    # annotation.transcript_sequence('ENSMUST00000115529.7', './genomes/mm10.fa', feature_type='cds')
    def transcript_sequence(self, transcript_id, assembly_fasta_fn, feature_type='exons'):
        iterator = self.transcript_sequences([transcript_id], assembly_fasta_fn, feature_type)
        transcript_id, sequence = next(iterator)
        return sequence

    # feature type is one of full/exons/cds
    def transcript_sequences(self, transcript_ids, assembly_fasta_fn, feature_type='exons'):
        intervals = itertools.chain.from_iterable(
            self.segments_as_bedtool_intervals(
                self.ordered_segments_by_type(transcript_id, feature_type),
                name = transcript_id
            ) for transcript_id in transcript_ids
        )
        bed = pybedtools.BedTool(intervals)
        seq_result = bed.sequence(fi = assembly_fasta_fn, name=True)
        with open(seq_result.seqfn) as fasta_file:
            sequences = iterate_as_fasta(fasta_file)
            sequence_groups = itertools.groupby(sequences, key=lambda hdr_seq_pair: hdr_seq_pair[0])
            # we have separate sequences for each feature (such as exon) of transcript, so we should join them
            for transcript_id, hdr_seq_pairs in sequence_groups:
                sequence = ''.join(hdr_seq_pair[1] for hdr_seq_pair in hdr_seq_pairs)
                yield (transcript_id, sequence)


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

def load_transcript_cds_info(cds_annotation_filename):
    transcript_infos = CodingTranscriptInfo.each_from_file(cds_annotation_filename)
    return {tr_info.transcript_id: tr_info  for tr_info in transcript_infos}

# takes an iterator of lines in FASTA file and yield pairs (header, sequence)
# properly handles multiline FASTA
def iterate_as_fasta(fasta_lines):
    header = None
    for line in fasta_lines:
        line = line.rstrip('\n')
        if line.startswith('>'):
            if header != None:
                yield (header, ''.join(sequences))
            header = line[1:].lstrip()
            sequences = []
        else:
            sequences.append(line.strip())
    if header != None:
        yield (header, ''.join(sequences))
