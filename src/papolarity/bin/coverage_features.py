import argparse
import numpy as np
from ..gzip_utils import open_for_write
from ..dto.transcript_coverage import TranscriptCoverage
from ..polarity_score import polarity_score

def configure_argparser(argparser=None):
    if not argparser:
        argparser = argparse.ArgumentParser(prog="coverage_features", description="Calculate coverage profile features")
    argparser.add_argument('coverage', metavar='coverage.bedgraph', help='Coverage data')
    argparser.add_argument('--prefix', default='', help='Prefix of feature columns (to distinguish samples)')
    argparser.add_argument('--output-file', '-o', dest='output_file', help="Store results at this path")
    return argparser

def main():
    argparser = configure_argparser()
    args = argparser.parse_args()
    invoke(args)

def invoke(args):
    coverages = TranscriptCoverage.each_in_file(args.coverage, header=False, dtype=int)
    with open_for_write(args.output_file) as output_stream:
        feature_names = ['mean_coverage', 'total_coverage', 'polarity']
        prefixed_feature_names = [f'{args.prefix}{name}' for name in feature_names]
        header = ['transcript_id', *prefixed_feature_names]
        print('\t'.join(header), file=output_stream)
        for transcript_coverage in coverages:
            coverage = transcript_coverage.coverage

            transcript_id = transcript_coverage.transcript_id
            mean_coverage = np.mean(coverage)
            total_coverage = np.sum(coverage)
            polarity = polarity_score(coverage)

            info = [transcript_id, mean_coverage, total_coverage, polarity]
            print('\t'.join(map(str, info)), file=output_stream)
