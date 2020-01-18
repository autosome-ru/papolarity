import argparse
from ..utils import common_subsequence
from ..gzip_utils import open_for_write
from ..dto.transcript_coverage import TranscriptCoverage
from ..segmentation import Segmentation
from ..dto.coverage_comparison_stats import CoverageComparisonStats

def compare_coverage_streams(segmentations, control_coverage_profiles, experiment_coverage_profiles):
    segmentation_and_coverage_profiles = [
        segmentations,
        control_coverage_profiles,
        experiment_coverage_profiles,
    ]

    key_extractors = [
        lambda segment: segment.chrom,
        lambda transcript_coverage: transcript_coverage.transcript_id,
        lambda transcript_coverage: transcript_coverage.transcript_id,
    ]

    transcript_stream = common_subsequence(segmentation_and_coverage_profiles, key=key_extractors, check_sorted=True)
    for (transcript_id, (segmentation, control_coverage, experiment_coverage)) in transcript_stream:
        yield CoverageComparisonStats.make_from_profiles(transcript_id, control_coverage.coverage, experiment_coverage.coverage, segmentation)

def configure_argparser(argparser=None):
    if not argparser:
        argparser = argparse.ArgumentParser(
            prog = "compare_coverage",
            description = "Coverage profile comparison",
        )
    argparser.add_argument('segmentation', metavar='segmentation.bed', help='Segmentation')
    argparser.add_argument('coverage_control', metavar='control.bedgraph', help='Coverage for control data')
    argparser.add_argument('coverage_experiment', metavar='experiment.bedgraph', help='Coverage for experiment data')
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

    with open_for_write(args.output_file) as output_stream:
        feature_names = ['slope', 'slopelog', 'l1_distance']
        prefixed_feature_names = [f'{args.prefix}{name}' for name in feature_names]
        header = ['transcript_id', *prefixed_feature_names]
        print('\t'.join(header), file=output_stream)
        for rec in compare_coverage_streams(segmentation_stream, control_coverage_profiles, experiment_coverage_profiles):
            info = [rec.transcript_id, rec.slope, rec.slopelog, rec.l1_distance]
            print('\t'.join(map(str, info)), file=output_stream)
