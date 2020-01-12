import argparse
from pybedtools import BedTool
from ..coverage_profile import make_coverage, coverage_intervals_from_bedgraph
from ..dto.coverage_interval import CoverageInterval

def configure_argparser(argparser=None):
    if not argparser:
        argparser = argparse.ArgumentParser(prog="get_coverage", description="Generates coverage from an alignment")
    argparser.add_argument('alignment', help='alignment in BAM format')
    argparser.add_argument('--sort', action='store_true', help="Sort resulting alignments by transcript name")
    argparser.add_argument('--output-file', '-o', dest='output_file', help="Store results at this path")
    argparser.add_argument('--dtype', choices=['int', 'float'], default='int', help="Make int or float-valued coverage (default: int)")
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

    alignment = BedTool(args.alignment)
    bedgraph = make_coverage(alignment, sort_transcripts=args.sort, stream=True, dtype=dtype)

    if args.output_file:
        bedgraph.saveas(args.output_file)
    else:
        intervals = coverage_intervals_from_bedgraph(bedgraph, dtype=dtype)
        CoverageInterval.print_tsv(intervals, header=False)
