from collections import defaultdict
from collections import namedtuple
from itertools import groupby
import math
import itertools
import numpy as np
import sys

import pysam
import os

import traceback
from gtf_parser import *

import sklearn.linear_model
from pasio import *
from annotation import *
import pasio_wrapper
from pooling import join_sorted, pooling

filename_1 = sys.argv[1]
filename_2 = sys.argv[2]

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

def coverages_in_bam(bam_fn):
    index_fn = os.path.abspath(bam_fn) + '.bai'
    if not os.path.isfile(index_fn):
        pysam.index(bam_fn)

    samfile = pysam.AlignmentFile(bam_fn, "rb")
    transcripts = { read.reference_name for read in samfile }
    for transcript in transcripts:
        yield (transcript, samfile.fetch(transcript))
        # transcript, reads_at_pos.reference_pos, reads_at_pos.get_num_aligned()


def coverage_list(reads, contig_length):
    # coverages_by_pos = {covered_pos.reference_pos: covered_pos.get_num_aligned()  for covered_pos in covered_positions}
    # if not contig_length:
    #     contig_length = max(coverages_by_pos.keys(), default=0)
    # return [coverages_by_pos.get(pos, 0) for pos in range(1, contig_length + 1)]
    profile = np.zeros(contig_length, dtype=np.int)
    for read in reads:
        profile[read.reference_start - 1 : read.reference_end] += 1
    return profile


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

annotation = cds_annotation(load_annotation_parts('gencode.vM22.basic.annotation.gtf.gz'))
print('Annotation loaded', file=sys.stderr)
basic_coverage_fields = [
    'gene_id', 'transcript_id',
    'transcript_length', 'cds_start', 'cds_stop',
    'read_number', 'mean_cds_coverage',
    'polarity_score',
]
basic_coverage_difference_fields = [
    'gene_id', 'transcript_id',
    'transcript_length', 'cds_start', 'cds_stop',
    'read_number_1', 'mean_cds_coverage_1', 'polarity_score_1',
    'read_number_2', 'mean_cds_coverage_2', 'polarity_score_2',
    'slope',
]
coverage_fields = basic_coverage_fields + ['stabilized_polarity_score', 'slope']
MeanTranscriptCoverage = namedtuple('MeanTranscriptCoverage', coverage_fields)
TranscriptCoverage = namedtuple('TranscriptCoverage', basic_coverage_fields + ['coverage', 'cds_coverage'])
CoverageDifference = namedtuple('CoverageDifference', basic_coverage_difference_fields)
splitter = pasio_wrapper.pasio_splitter()

def slope_for_segmentation(filename_1, filename_2, annotation, splitter):
    coverage_diffs = []
    num_errors = 0
    contigs_iterator_1 = iterate_transcript_coverages(filename_1, annotation)
    contigs_iterator_2 = iterate_transcript_coverages(filename_2, annotation)
    for (transcript_id, transcript_info_1, transcript_info_2) in join_sorted(contigs_iterator_1, contigs_iterator_2, key=lambda x: x.transcript_id):
        pooled_cds_coverage = pooling([transcript_info_1.cds_coverage, transcript_info_2.cds_coverage])
        total_coverage_1 = sum(transcript_info_1.cds_coverage)
        total_coverage_2 = sum(transcript_info_2.cds_coverage)
        # model = sklearn.linear_model.LinearRegression()
        # xs = []
        # ys = []
        # sample_weights = []
        # try:
        #     for (start, stop, _mean_count) in pasio_wrapper.stabile_segments(pooled_cds_coverage, splitter):
        #         mean_1 = np.mean(transcript_info_1.cds_coverage[start:stop])
        #         mean_2 = np.mean(transcript_info_2.cds_coverage[start:stop])
        #         coeff = ((mean_2 + 1) / (mean_1 + 1)) * ((total_coverage_1 + 1) / (total_coverage_2 + 1))
        #         coeff = math.log(coeff)
        #         xs.append([(start + stop) / 2])
        #         ys.append(coeff)
        #         sample_weights.append(1)
        #         # sample_weights.append(stop - start + 1)
        # except:
        #     num_errors += 1
        #     print('Error', file = sys.stderr)
        #     traceback.print_exc(file = sys.stderr)
        #     continue

        # model.fit(xs, ys, sample_weight=sample_weights)
        # slope = model.coef_[0]
        slope = 0

        field_values = [
            transcript_info_1.gene_id, transcript_id,
            transcript_info_1.transcript_length, transcript_info_1.cds_start, transcript_info_1.cds_stop,
            transcript_info_1.read_number, transcript_info_1.mean_cds_coverage, transcript_info_1.polarity_score,
            transcript_info_2.read_number, transcript_info_2.mean_cds_coverage, transcript_info_2.polarity_score,
            slope,
        ]
        coverage_diffs.append( CoverageDifference(*field_values) )
    print(f'Num errors: {num_errors}', file = sys.stderr)
    return coverage_diffs

def iterate_transcript_coverages(filename, annotation):
    # contigs_iterator = itertools.groupby(parse_coverages_file(filename), key=lambda coverage_pos: coverage_pos.contig)
    contigs_iterator = coverages_in_bam(filename)
    read_number = 0
    it = 0
    for transcript_id, reads in contigs_iterator:
        it += 1
        if it == 100:
            sys.exit(0)
        # print(transcript_id, file=sys.stderr)
        # Note that `reads` is an iterator and it behaves incorrect when converted to a list (yielded objects become staled)
        if transcript_id  not in  annotation['geneId_by_transcript']:
            continue
        print(f'{transcript_id} in annotation', file=sys.stderr)

        gene_id = annotation['geneId_by_transcript'][transcript_id]
        transcript_length = annotation['len_exons_by_transcript'][transcript_id]
        cds_start = annotation['cds_start_by_transcript'][transcript_id]
        cds_stop = annotation['cds_stop_by_transcript'][transcript_id]
        read_number += 1
        coverage = coverage_list(reads, contig_length = transcript_length)
        cds_coverage = coverage[cds_start - 1 : cds_stop]
        cds_length = cds_stop - cds_start + 1
        mean_cds_coverage = sum(cds_coverage) / cds_length
        polarity_score = calculate_polarity_score(cds_coverage)

        basic_field_values = [
            gene_id, transcript_id,
            transcript_length, cds_start, cds_stop,
            read_number, mean_cds_coverage,
            polarity_score,
        ]
        # print(TranscriptCoverage(*basic_field_values, coverage, cds_coverage), file=sys.stderr)
        yield TranscriptCoverage(*basic_field_values, coverage, cds_coverage)

def calculate_coverages(filename, annotation, splitter):
    mean_coverages = []
    num_errors = 0
    for transcript_info in iterate_transcript_coverages(filename, annotation):
        cds_coverage = transcript_info.cds_coverage
        stabilized_counts = np.zeros(len(cds_coverage))

        # model = sklearn.linear_model.LinearRegression()
        # xs = []
        # ys = []
        # sample_weights = []
        # try:
        #     for (start, stop, mean_count) in pasio_wrapper.stabile_segments(cds_coverage, splitter):
        #         stabilized_counts[start:stop] = mean_count
        #         xs.append([(start + stop) / 2])
        #         ys.append(mean_count)
        #         # sample_weights.append(stop - start + 1)
        #         sample_weights.append(1)
        # except:
        #     num_errors += 1
        #     print('Error', file = sys.stderr)
        #     traceback.print_exc(file = sys.stderr)
        #     continue

        # model.fit(xs, ys, sample_weight=sample_weights)
        # slope = model.coef_[0]
        slope = 0

        stabilized_cds_coverage = list(stabilized_counts)

        # stabilized_polarity_score = calculate_polarity_score(stabilized_cds_coverage)
        stabilized_polarity_score = 0
        field_values = list(transcript_info)[0:-2] + [stabilized_polarity_score, slope]
        mean_coverages.append( MeanTranscriptCoverage(*field_values) )
    print(f'Num errors: {num_errors}', file = sys.stderr)
    return mean_coverages

# mean_coverages_1 = calculate_coverages(filename_1, annotation, splitter)
# mean_coverages_2 = calculate_coverages(filename_2, annotation, splitter)

# best_coverage_1_by_gene = {}
# best_coverage_2_by_gene = {}

# for gene_id, transcript_infos in itertools.groupby(sorted(mean_coverages_1, key=lambda rec: rec.gene_id), key=lambda rec: rec.gene_id):
#     transcript_info = max(transcript_infos, key = lambda rec: rec.mean_cds_coverage)
#     best_coverage_1_by_gene[gene_id] = transcript_info

# for gene_id, transcript_infos in itertools.groupby(sorted(mean_coverages_2, key=lambda rec: rec.gene_id), key=lambda rec: rec.gene_id):
#     transcript_info = max(transcript_infos, key = lambda rec: rec.mean_cds_coverage)
#     best_coverage_2_by_gene[gene_id] = transcript_info

# header = [f'{field}_1' for field in coverage_fields] + [f'{field}_2' for field in coverage_fields] + ['delta_angle']
# print('\t'.join(header))

# common_genes = set(best_coverage_1_by_gene.keys()).intersection(set(best_coverage_2_by_gene.keys()))

# for gene_id in common_genes:
#     transcript_info_1 = best_coverage_1_by_gene[gene_id]
#     transcript_info_2 = best_coverage_2_by_gene[gene_id]
#     if transcript_info_1.transcript_id != transcript_info_2.transcript_id:
#         continue
#     delta_angle = math.atan(transcript_info_2.slope) - math.atan(transcript_info_1.slope)
#     row = [transcript_info_1._asdict()[field] for field in coverage_fields] + \
#         [transcript_info_2._asdict()[field] for field in coverage_fields] + \
#         [delta_angle]
#     print('\t'.join(map(str, row)))

coverage_diffs = slope_for_segmentation(filename_1, filename_2, annotation, splitter)

best_coverage_1_by_gene = {}
best_coverage_2_by_gene = {}

print(*basic_coverage_difference_fields, sep='\t')
for gene_id, transcript_diff_infos in itertools.groupby(sorted(coverage_diffs, key=lambda rec: rec.gene_id), key=lambda rec: rec.gene_id):
    transcript_info = max(transcript_diff_infos, key = lambda rec: (rec.mean_cds_coverage_1 + rec.mean_cds_coverage_2) / 2)
    print(*transcript_info, sep='\t')
