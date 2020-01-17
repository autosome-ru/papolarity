from collections import defaultdict
import itertools
import pybedtools
from .gtf_parser import GTFRecord
from .utils import take_the_only
from .fasta_reader import fasta_from_file
from .dto.coding_transcript_info import CodingTranscriptInfo

class Annotation:
    def __init__(self):
        self.gene_by_id = {}
        self.transcript_by_id = {}
        self.geneId_by_transcript = {}
        self.transcripts_by_gene = defaultdict(list)
        self.parts_by_transcript = defaultdict(list)

    @classmethod
    def load(cls, filename, relevant_attributes=None, multivalue_keys=frozenset(), ignore_unknown_multivalues=False, condition=None):
        annotation = cls()
        records = GTFRecord.each_in_file(filename, multivalue_keys=multivalue_keys, ignore_unknown_multivalues=ignore_unknown_multivalues)

        if relevant_attributes is not None:
            # that's the minimum list of attributes for a library to properly work
            necessary_attributes = {'gene_id', 'transcript_id',}
            relevant_attributes = relevant_attributes.union(necessary_attributes)
            records = map(lambda rec: rec.attributes_filtered(relevant_attributes), records)
        if condition:
            records = filter(condition, records)

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

    def transcript_start_codons(self, transcript_id):
        '''
        We return a list of codons because some transcripts
        don't have an annotated start/stop codon or have
        several segments due to splicing
        '''
        parts = self.parts_by_transcript[transcript_id]
        start_codons = [part  for part in parts  if part.type == 'start_codon']
        return start_codons

    def transcript_stop_codons(self, transcript_id):
        '''
        We return a list of codons because some transcripts
        don't have an annotated start/stop codon or have
        several segments due to splicing
        '''
        parts = self.parts_by_transcript[transcript_id]
        stop_codons = [part  for part in parts  if part.type == 'stop_codon']
        return stop_codons

    def transcript_utrs(self, transcript_id):
        parts = self.parts_by_transcript[transcript_id]
        return [part  for part in parts  if part.type == 'UTR']

    def transcript_strand(self, transcript_id):
        return self.segments_strand(self.parts_by_transcript[transcript_id])

    @classmethod
    def segments_strand(cls, segments):
        strands = { segment.strand  for segment in  segments }
        if len(strands) != 1:
            raise ValueError('Different strands')
        return strands.pop()

    @classmethod
    def segments_ordered_5_to_3(cls, segments):
        strand = cls.segments_strand(segments)
        reverse_order = (strand == '-')
        by_start_coordinate = lambda segment: segment.start
        return sorted(segments, key=by_start_coordinate, reverse=reverse_order)

    def coding_transcript_info(self, transcript_id):
        """Return coding transcript info. This method fails when called for non-coding transcript"""
        gene_id = self.geneId_by_transcript[transcript_id]
        exons = self.transcript_exons(transcript_id)
        cds_segments = self.transcript_cds(transcript_id)
        strand = self.transcript_strand(transcript_id)

        transcript_length = sum(exon.length for exon in exons)

        if len(cds_segments) == 0:
            return CodingTranscriptInfo(gene_id, transcript_id, transcript_length, '', '')

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
        if feature_type == 'exons':
            return self.segments_ordered_5_to_3(self.transcript_exons(transcript_id))
        elif feature_type == 'cds':
            return self.segments_ordered_5_to_3(self.transcript_cds(transcript_id))
        elif feature_type == 'cds_with_stop':
            return self.segments_ordered_5_to_3(self.transcript_cds(transcript_id) + self.transcript_stop_codons(transcript_id))
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
        seq_result = bed.sequence(fi = assembly_fasta_fn, name=True, s=True)  # s - strandedness
        sequences = fasta_from_file(seq_result.seqfn)
        sequence_groups = itertools.groupby(sequences, key=lambda hdr_seq_pair: hdr_seq_pair[0])
        # we have separate sequences for each feature (such as exon) of transcript, so we should join them
        for transcript_id, hdr_seq_pairs in sequence_groups:
            if transcript_id[-3:] in ['(+)', '(-)']:
                transcript_id = transcript_id[0:-3]
            sequence = ''.join(hdr_seq_pair[1] for hdr_seq_pair in hdr_seq_pairs)
            yield (transcript_id, sequence)
