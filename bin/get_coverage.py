#!/usr/bin/env python
import sys
from os.path import dirname
sys.path.insert(0, dirname(dirname(__file__)))
import argparse
import coverage_profile
from pybedtools import BedTool

def get_argparser():
    argparser = argparse.ArgumentParser(
        # prog = "get_coverage",
        description = "Generates coverage from an alignment",
    )
    argparser.add_argument('alignment', help='alignment in BAM format')
    argparser.add_argument('--sort', action='store_true', help="Sort resulting alignments by transcript name")
    argparser.add_argument('--output-file', help="Store bedgraph at this path")
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

alignment = BedTool(args.alignment)
bedgraph = coverage_profile.make_coverage(alignment, sort_transcripts=args.sort, stream=True, dtype=dtype)

if args.output_file:
    bedgraph.saveas(args.output_file)
else:
    for interval in bedgraph:
        print('\t'.join(interval[0:4]))
