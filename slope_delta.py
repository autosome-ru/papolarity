from collections import defaultdict
from collections import namedtuple
from itertools import groupby
import math
import itertools
import numpy as np
import sys
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

annotation = cds_annotation(load_annotation_parts('gencode.vM22.basic.annotation.gtf.gz'))
coverage_fields = ['gene_id', 'transcript_id', 'transcript_length', 'cds_start', 'cds_stop', 'read_number', 'mean_cds_coverage', 'polarity_score', 'stabilized_polarity_score', 'slope']
MeanTranscriptCoverage = namedtuple('MeanTranscriptCoverage', coverage_fields)

splitter = pasio_wrapper.pasio_splitter()
def calculate_coverages(filename, annotation, splitter):
    contigs_iterator = itertools.groupby(parse_coverages_file(filename), key=lambda coverage_pos: coverage_pos.contig)
    mean_coverages = []
    num_errors = 0
    for transcript_id, contig_coverages in contigs_iterator:
        if transcript_id  not in  annotation['geneId_by_transcript']:
            continue
        contig_coverages = list(contig_coverages)

        gene_id = annotation['geneId_by_transcript'][transcript_id]
        transcript_length = annotation['len_exons_by_transcript'][transcript_id]
        cds_start = annotation['cds_start_by_transcript'][transcript_id]
        cds_stop = annotation['cds_stop_by_transcript'][transcript_id]
        read_number = len(contig_coverages)
        coverage = coverage_list(contig_coverages, contig_length = transcript_length)
        cds_coverage = coverage[cds_start - 1 : cds_stop]
        cds_length = cds_stop - cds_start + 1
        mean_cds_coverage = sum(cds_coverage) / cds_length

        try:
            stabilized_counts = np.zeros(len(cds_coverage))

            model = sklearn.linear_model.LinearRegression()
            xs = []
            ys = []
            sample_weights = []

            for (start, stop, mean_count) in pasio_wrapper.stabile_segments(cds_coverage, splitter):
                stabilized_counts[start:stop] = mean_count
                xs.append([(start + stop) / 2])
                ys.append(mean_count)
                sample_weights.append(stop - start + 1)

            model.fit(xs, ys, sample_weight=sample_weights)
            slope = model.coef_[0]

            stabilized_cds_coverage = list(stabilized_counts)

            polarity_score = calculate_polarity_score(cds_coverage)
            stabilized_polarity_score = calculate_polarity_score(stabilized_cds_coverage)
            field_values = [gene_id, transcript_id, transcript_length, cds_start, cds_stop, read_number, mean_cds_coverage, polarity_score, stabilized_polarity_score, slope]
            mean_coverages.append( MeanTranscriptCoverage(*field_values) )
        except:
            num_errors += 1
            print('Error', file = sys.stderr)
            traceback.print_exc(file = sys.stderr)
    print(f'Num errors: {num_errors}', file = sys.stderr)
    return mean_coverages

mean_coverages_1 = calculate_coverages(filename_1, annotation, splitter)
mean_coverages_2 = calculate_coverages(filename_2, annotation, splitter)

best_coverage_1_by_gene = {}
best_coverage_2_by_gene = {}

for gene_id, transcript_infos in itertools.groupby(sorted(mean_coverages_1, key=lambda rec: rec.gene_id), key=lambda rec: rec.gene_id):
    transcript_info = max(transcript_infos, key = lambda rec: rec.mean_cds_coverage)
    best_coverage_1_by_gene[gene_id] = transcript_info

for gene_id, transcript_infos in itertools.groupby(sorted(mean_coverages_2, key=lambda rec: rec.gene_id), key=lambda rec: rec.gene_id):
    transcript_info = max(transcript_infos, key = lambda rec: rec.mean_cds_coverage)
    best_coverage_2_by_gene[gene_id] = transcript_info

header = [f'{field}_1' for field in coverage_fields] + [f'{field}_2' for field in coverage_fields] + ['delta_angle']
print('\t'.join(header))

common_genes = set(best_coverage_1_by_gene.keys()).intersection(set(best_coverage_2_by_gene.keys()))

for gene_id in common_genes:
    transcript_info_1 = best_coverage_1_by_gene[gene_id]
    transcript_info_2 = best_coverage_2_by_gene[gene_id]
    if transcript_info_1.transcript_id != transcript_info_2.transcript_id:
        continue
    delta_angle = math.atan(transcript_info_2.slope) - math.atan(transcript_info_1.slope)
    row = [transcript_info_1._asdict()[field] for field in coverage_fields] + \
        [transcript_info_2._asdict()[field] for field in coverage_fields] + \
        [delta_angle]
    print('\t'.join(map(str, row)))
