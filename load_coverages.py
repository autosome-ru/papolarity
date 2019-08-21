import itertools
from collections import namedtuple

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

filename = 'METTL3_coverage.tsv'
header = ['gene_id', 'transcript_id', 'mean_coverage', 'polarity_score']
print('\t'.join(header))
contigs_iterator = itertools.groupby(parse_coverages_file(filename), key=lambda coverage_pos: coverage_pos.contig)
for contig, contig_coverages in contigs_iterator:
    # print(len(list(contig_coverages)))
    contig_coverages = list(contig_coverages)
    read_number = len(contig_coverages)
    coverage = coverage_list(contig_coverages)
    print(contig, read_number, sep="\t")
    print(coverage)
