import argparse
import numpy as np
from ..gzip_utils import open_for_write
from ..utils import tsv_string_empty_none
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
    coverage_profiles = TranscriptCoverage.each_in_file(args.coverage, header=False, dtype=int)
    with open_for_write(args.output_file) as output_stream:
        # Note: 'q50' etc goes as the last part of name because csvtk-0.19.1
        # filter2 function had some problems with column names containing digits in the middle of the name.
        # The problem was resolved in 0.19.2 version.
        # (see https://github.com/shenwei356/csvtk/issues/44 for details)
        feature_names = ['mean_coverage', 'coverage_q25', 'coverage_q50', 'coverage_q75', 'total_coverage', 'polarity']
        prefixed_feature_names = [f'{args.prefix}{name}' for name in feature_names]
        header = ['transcript_id', *prefixed_feature_names]
        print('\t'.join(header), file=output_stream)
        for transcript_coverage in coverage_profiles:
            coverage = transcript_coverage.coverage

            transcript_id = transcript_coverage.transcript_id
            mean_coverage = np.mean(coverage)
            coverage_q25 = np.quantile(coverage, 0.25)
            coverage_q50 = np.quantile(coverage, 0.5)
            coverage_q75 = np.quantile(coverage, 0.75)
            total_coverage = np.sum(coverage)
            polarity = polarity_score(coverage)

            info = [transcript_id, mean_coverage, coverage_q25, coverage_q50, coverage_q75, total_coverage, polarity]
            print(tsv_string_empty_none(info), file=output_stream)
