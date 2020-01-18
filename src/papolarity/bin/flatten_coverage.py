import argparse
import math
from ..utils import align_iterators, get_constant_intervals
from ..gzip_utils import open_for_write
from ..dto.coverage_interval import CoverageInterval
from ..dto.transcript_coverage import TranscriptCoverage
from ..segmentation import Segmentation

def configure_argparser(argparser=None):
    if not argparser:
        argparser = argparse.ArgumentParser(
            prog = "flatten_coverage",
            description = "Flatten coverage profiles by averaging data through given segments",
        )
    argparser.add_argument('segmentation', help='Segmentation of contigs in bed format')
    argparser.add_argument('coverage', help='Coverage in bedgraph format')
    argparser.add_argument('--only-matching', action='store_true', help="Don't pool coverage profiles of transcripts which are present not in all files")
    argparser.add_argument('--output-file', '-o', dest='output_file', help="Store results at this path")
    argparser.add_argument('--rounding', choices=['no', 'round', 'ceil', 'floor'], default='none', help="Rounding of float values (default: no rounding)")
    return argparser

def main():
    argparser = configure_argparser(argparser)
    args = argparser.parse_args()
    invoke(args)

def invoke(args):
    if args.rounding == 'round':
        rounding = round
    elif args.rounding == 'ceil':
        rounding = math.ceil
    elif args.rounding == 'floor':
        rounding = math.floor
    elif args.rounding in ['no', 'none']:
        rounding = lambda x: x
    else:
        raise ValueError('rounding should be either int or float')

    transcript_coverage_stream = TranscriptCoverage.each_in_file(args.coverage, header=False, dtype=float)

    segmentation_stream = Segmentation.each_in_file(args.segmentation, header=False)

    with open_for_write(args.output_file) as output_stream:
        aligned_transcripts = align_iterators([segmentation_stream, transcript_coverage_stream], key=[lambda segment: segment.chrom, lambda transcript_coverage: transcript_coverage.transcript_id], check_sorted=True)
        for (transcript_id, (segmentation, transcript_coverage)) in aligned_transcripts:
            if transcript_coverage is None:
                continue
            if (segmentation is None) and args.only_matching:
                continue
            if segmentation is not None:
                transcript_coverage = segmentation.stabilize_profile(transcript_coverage.coverage)
            coverage_intervals = get_constant_intervals(transcript_coverage)
            bedgraph = (CoverageInterval(transcript_id, *interval, dtype=float) for interval in coverage_intervals)
            CoverageInterval.print_tsv(bedgraph, header=False, file=output_stream)
