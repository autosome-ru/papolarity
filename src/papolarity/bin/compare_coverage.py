import argparse
from ..utils import tsv_string_empty_none
from ..gzip_utils import open_for_write
from ..dto.transcript_coverage import TranscriptCoverage
from ..segmentation import Segmentation
from ..profile_comparison import compare_coverage_streams

def configure_argparser(argparser=None):
    if not argparser:
        argparser = argparse.ArgumentParser(
            prog = "compare_coverage",
            description = "Coverage profile comparison",
        )
    argparser.add_argument('segmentation', metavar='segmentation.bed', help='Segmentation')
    argparser.add_argument('coverage_control', metavar='control.bedgraph', help='Coverage for control data')
    argparser.add_argument('coverage_experiment', metavar='experiment.bedgraph', help='Coverage for experiment data')
    argparser.add_argument('--segment-coverage-quantile', nargs=2, metavar=('<quantile>', '<threshold>'), default=['0.5', '0'],
                           help='Filter transcripts with too low value of quantile over segments.\n'
                                'E.g. `--segment-coverage-quantile  0.5  1.0` filters out all transcripts with median '
                                'segment coverage less than 1.0')
    argparser.add_argument('--prefix', default='', help='Prefix of feature columns (to distinguish samples)')
    argparser.add_argument('--output-file', '-o', dest='output_file', help="Store results at this path")
    return argparser

def main():
    argparser = configure_argparser()
    args = argparser.parse_args()
    invoke(args)

def invoke(args):
    segmentation_stream = Segmentation.each_in_file(args.segmentation, header=False)

    control_coverage_profiles = TranscriptCoverage.each_in_file(args.coverage_control, header=False, dtype=int)
    experiment_coverage_profiles = TranscriptCoverage.each_in_file(args.coverage_experiment, header=False, dtype=int)

    quantile_q, quantile_threshold = [float(x) for x in args.segment_coverage_quantile]

    with open_for_write(args.output_file) as output_stream:
        feature_names = ['slope', 'slopelog', 'l1_distance', 'polarity_diff']
        prefixed_feature_names = [f'{args.prefix}{name}' for name in feature_names]
        header = ['transcript_id', *prefixed_feature_names]
        print('\t'.join(header), file=output_stream)
        for rec in compare_coverage_streams(segmentation_stream, control_coverage_profiles, experiment_coverage_profiles,
                                            quantile_q=quantile_q, quantile_threshold=quantile_threshold):
            info = [rec[field] for field in ['transcript_id', *feature_names]]
            print(tsv_string_empty_none(info), file=output_stream)
