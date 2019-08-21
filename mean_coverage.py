from collections import defaultdict
from collections import namedtuple
from itertools import groupby
import itertools
import numpy as np
import sys
from gtf_parser import *

import sklearn.linear_model
from pasio import *

filename = sys.argv[1]

alpha = 1
beta = 1
scorer_factory = ScorerFactory(alpha, beta)
length_regularization_multiplier = 0
length_regularization_function = lambda x: x
split_number_regularization_multiplier = 0
algorithm = 'rounds'
no_split_constant = True
window_size = 2500
window_shift = 1250
num_rounds = None
split_at_gaps = False

square_splitter = SquareSplitter(scorer_factory,
    length_regularization_multiplier = 0,
    length_regularization_function = lambda x: x,
    split_number_regularization_multiplier = 0,
    split_number_regularization_function = lambda x: x)

if algorithm == 'exact':
    splitter = square_splitter
else:
    if no_split_constant:
        square_splitter = ReducerCombiner(NotConstantReducer(), square_splitter)
    else:
        square_splitter = ReducerCombiner(NotZeroReducer(), square_splitter)
    ###
    if algorithm == 'slidingwindow':
        sliding_window = SlidingWindow(window_size = window_size, window_shift = window_shift)
        reducer = SlidingWindowReducer(sliding_window = sliding_window, base_reducer = square_splitter)
        splitter = ReducerCombiner(reducer, square_splitter)
    elif algorithm == 'rounds':
        sliding_window = SlidingWindow(window_size = window_size, window_shift = window_shift)
        base_reducer = SlidingWindowReducer(sliding_window = sliding_window, base_reducer = square_splitter)
        reducer = RoundReducer(base_reducer = base_reducer, num_rounds = num_rounds)
        splitter = ReducerCombiner(reducer, NopSplitter(scorer_factory))


def take_the_only(arr):
    if len(arr) > 1:
        raise Exception('Several elements when the only one is expected')
    if len(arr) == 0:
        raise Exception('No elements when one is expected')
    return arr[0]

def get_transcript_id(rec):
    return rec.attributes['transcript_id']

CoveredPosition = namedtuple('CoveredPosition', ['contig', 'pos', 'coverage'])

def parse_row(l):
    contig, pos, coverage = l.strip().split("\t")
    return CoveredPosition(contig, int(pos), int(coverage))

def parse_coverages_file(filename):
    with open(filename) as f:
        for row in f:
            yield parse_row(row)

def coverage_list(covered_positions, contig_length=None):
    coverages_by_pos = {covered_pos.pos: covered_pos.coverage  for covered_pos in covered_positions}
    if not contig_length:
        contig_length = max(coverages_by_pos.keys(), default=0)
    return [coverages_by_pos.get(pos, 0) for pos in range(1, contig_length + 1)]


def calculate_polarity_score(m):
    total_m = 0
    numerator = 0
    for i in range(len(m)):
        total_m = total_m + m[i]
        numerator += m[i] * i
    if total_m == 0:
        return 0
    polarity_score = numerator / total_m
    normalized_score = (polarity_score * 2 / (len(m) - 1)) - 1
    return normalized_score

genes = {}
gene_by_transcript = {}
transcripts_by_gene = defaultdict(list)
parts_by_transcript = defaultdict(list)
for rec in parse_gtf('gencode.vM22.basic.annotation.gtf.gz'):
# for rec in parse_gtf('annotation_large_sample.gtf'):
    if rec.attributes['gene_type'] != 'protein_coding':
        continue

    gene_id = rec.attributes['gene_id']
    if rec.type == 'gene':
        genes[gene_id] = rec

    else:
        if rec.attributes['transcript_type'] != 'protein_coding':
            continue

        transcript_id = rec.attributes['transcript_id']
        if rec.type == 'transcript':
            gene_by_transcript[transcript_id] = gene_id
            transcripts_by_gene[gene_id].append(rec)
        else:
            parts_by_transcript[transcript_id].append(rec)


len_exons_by_transcript = {}
len_cds_by_transcript = {}
cds_start_by_transcript = {}
cds_stop_by_transcript = {}
for gene_id, gene in genes.items():
    for transcript in transcripts_by_gene[gene_id]:
        transcript_id = transcript.attributes['transcript_id']
        parts = parts_by_transcript[transcript_id]

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

            # print(transcript_id, cds_start_transcript, sep='\t')
            # ToDo: посчитать cds_stop_transcript
        # ToDo: посчитать start и stop для отрицательной нитки

        
        len_exons_by_transcript[transcript_id] = len_transcript
        len_cds_by_transcript[transcript_id] = len_cds_in_transcript
        cds_start_by_transcript[transcript_id] = cds_start_transcript
        cds_stop_by_transcript[transcript_id] = cds_stop_transcript
        # print(gene_id, transcript_id, len_transcript, len_cds_in_transcript, sep='\t')



contigs_iterator = itertools.groupby(parse_coverages_file(filename), key=lambda coverage_pos: coverage_pos.contig)

coverage_fields = ['gene_id', 'transcript_id', 'transcript_length', 'cds_start', 'cds_stop', 'read_number', 'mean_cds_coverage', 'polarity_score', 'stabilized_polarity_score', 'slope']
MeanTranscriptCoverage = namedtuple('MeanTranscriptCoverage', coverage_fields)
mean_coverages = []
num_errors = 0
for transcript_id, contig_coverages in contigs_iterator:
    if transcript_id  not in  gene_by_transcript:
        continue
    contig_coverages = list(contig_coverages)

    gene_id = gene_by_transcript[transcript_id]
    transcript_length = len_exons_by_transcript[transcript_id]
    cds_start = cds_start_by_transcript[transcript_id]
    cds_stop = cds_stop_by_transcript[transcript_id]
    read_number = len(contig_coverages)
    coverage = coverage_list(contig_coverages, contig_length = transcript_length)
    cds_coverage = coverage[cds_start - 1 : cds_stop]
    cds_length = cds_stop - cds_start + 1
    mean_cds_coverage = sum(cds_coverage) / cds_length
    
    try:
        chrom_start = 1
        counts = np.array(cds_coverage)
        split_candidates = np.arange(len(counts) + 1)
        score, splits = splitter.split(counts, split_candidates)
        scorer = splitter.scorer(counts, splits)
        sum_logfac = scorer.total_sum_logfac()
        log_likelyhood = score - sum_logfac
        stabilized_counts = np.zeros(len(counts))

        model = sklearn.linear_model.LinearRegression()
        xs = []
        ys = []
        sample_weights=[]

        for (start, stop, mean_count) in zip(splits, splits[1:], scorer.mean_counts()):
            stabilized_counts[start:stop] = mean_count
            xs.append([(start + stop) / 2])
            ys.append(mean_count)
            sample_weights.append(stop - start + 1)

        model.fit(xs, ys, sample_weight=sample_weights)
        slope = model.coef_[0]

        stabilized_cds_coverage = list(stabilized_counts)
        # print('\t'.join(map(str,list(counts))))
        # print('\t'.join(map(str,list(stabilized_counts))))

        # polarity_score = calculate_polarity_score(coverage)
        polarity_score = calculate_polarity_score(cds_coverage)
        stabilized_polarity_score = calculate_polarity_score(stabilized_cds_coverage)
        field_values = [gene_id, transcript_id, transcript_length, cds_start, cds_stop, read_number, mean_cds_coverage, polarity_score, stabilized_polarity_score, slope]
        mean_coverages.append( MeanTranscriptCoverage(*field_values) )
    except:
        num_errors += 1
        print('Errors', file = sys.stderr)
    # print(transcript_id, read_number, sep="\t")
    # print(coverage)
print(f'Num errors: {num_errors}', file = sys.stderr)

print('\t'.join(coverage_fields))
for gene_id, transcript_infos in itertools.groupby(sorted(mean_coverages, key=lambda rec: rec.gene_id), key=lambda rec: rec.gene_id):
    transcript_infos = list(transcript_infos)
    # print(list(sorted(transcript_infos, key = lambda rec: rec.mean_cds_coverage, reverse = True))[0])
        #list делает список из коллекции, sorted по умолчанию - по возрастанию. Так что reverse. [0] - выбор первой строки
    transcript_info = max(transcript_infos, key = lambda rec: rec.mean_cds_coverage)
    row = [transcript_info._asdict()[field] for field in coverage_fields]
    print('\t'.join(map(str, row)))

    # for transcript_info in sorted(transcript_infos, key = lambda rec: rec.mean_cds_coverage, reverse = True):
    #     print(gene_id, transcript_info.transcript_id, transcript_info.mean_cds_coverage, sep="\t")

