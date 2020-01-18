import argparse
import numpy as np
from ..utils import align_iterators, get_constant_intervals
from ..dto.coverage_interval import CoverageInterval
from ..dto.transcript_coverage import TranscriptCoverage
from ..gzip_utils import open_for_write

def configure_argparser(argparser=None):
    if not argparser:
        argparser = argparse.ArgumentParser(prog="pool_coverage", description="Pool coverage profiles")
    argparser.add_argument('coverage_profiles', nargs='*', help='Coverage profiles in bedgraph format')
    argparser.add_argument('--only-matching', action='store_true', help="Don't pool coverage profiles of transcripts which are present not in all files")
    argparser.add_argument('--output-file', '-o', dest='output_file', help="Store results at this path")
    argparser.add_argument('--output-mode', choices=['sum', 'mean'], default='sum', help="What to report")
    argparser.add_argument('--dtype', choices=['int', 'float'], default='int', help="Make int or float-valued coverage (default: %(default)s)")
    return argparser

def main():
    argparser = configure_argparser(argparser)
    args = argparser.parse_args()
    invoke(args)

def invoke(args):
    if args.dtype == 'int':
        dtype = int
    elif args.dtype == 'float':
        dtype = float
    else:
        raise ValueError('dtype should be either int or float')

    coverage_streams = [TranscriptCoverage.each_in_file(coverage_fn, header=False, dtype=dtype) for coverage_fn in args.coverage_profiles]

    with open_for_write(args.output_file) as output_stream:
        aligned_transcripts = align_iterators(coverage_streams, key=lambda transcript_coverage: transcript_coverage.transcript_id, check_sorted=True)
        for (transcript_id, transcript_coverage_profiles) in aligned_transcripts:
            if args.only_matching and not all(transcript_coverage_profiles):
                continue
            transcript_coverage_profiles = filter(lambda x: x, transcript_coverage_profiles)
            coverage_values = [transcript_coverage.coverage for transcript_coverage in transcript_coverage_profiles]
            if args.output_mode == 'sum':
                pooled_coverage_profile = np.sum(coverage_values, axis=0)
            elif args.output_mode == 'mean':
                pooled_coverage_profile = np.mean(coverage_values, axis=0)
            else:
                raise ValueError(f'Unknown output_mode `{args.output_mode}`')
            pooled_coverage_intervals = get_constant_intervals(pooled_coverage_profile)
            pooled_bedgraph = (CoverageInterval(transcript_id, *interval, dtype=dtype) for interval in pooled_coverage_intervals)
            CoverageInterval.print_tsv(pooled_bedgraph, header=False, file=output_stream)
