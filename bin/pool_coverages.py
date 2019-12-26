#!/usr/bin/env python
import sys
from os.path import dirname
sys.path.insert(0, dirname(dirname(__file__)))
import argparse
import gzip
from utils import align_iterators, pool_profiles, get_constant_intervals
from dto.coverage_interval import CoverageInterval
from dto.transcript_coverage import TranscriptCoverage
from gzip_utils import open_for_write

def get_argparser():
    argparser = argparse.ArgumentParser(
        prog = "pool_coverage",
        description = "Pool coverage profiles",
    )
    argparser.add_argument('coverages', nargs='*', help='Alignment coverages in bedgraph format')
    argparser.add_argument('--only-matching', action='store_true', help="Don't pool coverages of transcripts which are present not in all files")
    argparser.add_argument('--output-file', '-o', dest='output_file', help="Store results at this path")
    argparser.add_argument('--dtype', choices=['int', 'float'], default='int', help="Make int or float-valued coverage (default: int)")
    return argparser

argparser = get_argparser()
args = argparser.parse_args()

if args.dtype == 'int':
    dtype = int
elif args.dtype == 'float':
    dtype = float
else:
    raise ValueError('dtype should be either int or float')

coverage_streams = [TranscriptCoverage.each_in_file(coverage_fn, header=False, dtype=dtype) for coverage_fn in args.coverages]

with open_for_write(args.output_file) as output_stream:
    aligned_transcripts = align_iterators(coverage_streams, key=lambda transcript_coverage: transcript_coverage.transcript_id, check_sorted=True)
    for (transcript_id, coverages) in aligned_transcripts:
        if args.only_matching and not all(coverages):
            continue
        pooled_coverage_profile = pool_profiles([transcript_coverage.coverage for transcript_coverage in coverages])
        pooled_coverage_intervals = get_constant_intervals(pooled_coverage_profile)
        pooled_bedgraph = (CoverageInterval(transcript_id, *interval, dtype=dtype) for interval in pooled_coverage_intervals)
        CoverageInterval.print_tsv(pooled_bedgraph, header=False, file=output_stream)
