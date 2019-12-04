#!/usr/bin/env python
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
    return argparser

argparser = get_argparser()
args = argparser.parse_args()

alignment = BedTool(args.alignment)
bedgraph = coverage_profile.make_coverage(alignment, sort_transcripts=args.sort, stream=True)

if args.output_file:
    bedgraph.saveas(args.output_file)
else:
    for interval in bedgraph:
        print('\t'.join(interval[0:4]))
