import sys
from os.path import dirname
sys.path.insert(0, dirname(dirname(__file__)))
import argparse
import numpy as np
from gzip_utils import open_for_write
from dto.transcript_coverage import TranscriptCoverage
from polarity_score import polarity_score

def get_argparser():
    argparser = argparse.ArgumentParser(
        prog = "coverage_properties",
        description = "Calculate coverage profile properties",
    )
    argparser.add_argument('coverage', metavar='coverage.bedgraph', help='Coverage data')
    argparser.add_argument('--output-file', '-o', dest='output_file', help="Store results at this path")
    return argparser

argparser = get_argparser()
args = argparser.parse_args()

coverages = TranscriptCoverage.each_in_file(args.coverage, header=False, dtype=int)

with open_for_write(args.output_file) as output_stream:
    header = ['transcript_id', 'transcript_length', 'mean_coverage', 'total_coverage', 'polarity']
    print('\t'.join(header), file=output_stream)
    for transcript_coverage in coverages:
        coverage = transcript_coverage.coverage

        transcript_id = transcript_coverage.transcript_id
        transcript_length = len(coverage)
        mean_coverage = np.mean(coverage)
        total_coverage = np.sum(coverage)
        polarity = polarity_score(coverage)

        info = [transcript_id, transcript_length, mean_coverage, total_coverage, polarity]
        print('\t'.join(map(str, info)), file=output_stream)
