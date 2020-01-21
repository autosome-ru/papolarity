import argparse
import numpy as np
from math import log2
from ..utils import tsv_string_empty_none
from ..gzip_utils import open_for_write
from ..dto.transcript_coverage import TranscriptCoverage
from ..segmentation import Segmentation
from ..dto.coverage_comparison_stats import CoverageComparisonStats
from ..profile_comparison import align_profile_streams_to_segmentation, segmentwise_sums

def configure_argparser(argparser=None):
    if not argparser:
        argparser = argparse.ArgumentParser(
            prog = "compare_coverage_extract",
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

    feature_names = [
        'segment_start', 'segment_stop', 'coord', 'rel_coord', \
        'control_segment_sum', 'experiment_segment_sum', \
        'mean_control', 'mean_experiment', \
        'ratio', 'log2ratio',
    ]
    with open_for_write(args.output_file) as output_stream:
        prefixed_feature_names = [f'{args.prefix}{name}' for name in feature_names]
        header = ['transcript_id', *prefixed_feature_names]
        print('\t'.join(header), file=output_stream)
        aligned_stream = align_profile_streams_to_segmentation(
            segmentation_stream, 
            [control_coverage_profiles, experiment_coverage_profiles]
        )
        for (transcript_id, (segmentation, control_coverage, experiment_coverage)) in aligned_stream:
            control_sums = segmentwise_sums(segmentation, control_coverage.coverage)
            experiment_sums = segmentwise_sums(segmentation, experiment_coverage.coverage)

            for row in slope_segments(control_sums, experiment_sums, segmentation):
                print(tsv_string_empty_none([transcript_id, *row]), file=output_stream)

def slope_segments(control_sums, experiment_sums, segmentation):
    assert len(segmentation.segments) == len(control_sums) == len(experiment_sums)
    profile_len = segmentation.segmentation_length
    num_segments = segmentation.num_segments

    total_coverage_control = np.sum(control_sums)
    total_coverage_experiment = np.sum(experiment_sums)
    if (total_coverage_control == 0) or (total_coverage_experiment ==0):
        pass # yield nothing for uncovered transcripts

    for (segment, control_segment_sum, experiment_segment_sum) in zip(segmentation.segments, control_sums, experiment_sums):
        coord = (segment.start + segment.stop - 1) / 2
        rel_coord = coord / profile_len

        if (control_segment_sum == 0) and (experiment_segment_sum == 0):
            yield (
                segment.start, segment.stop, coord, rel_coord, \
                control_segment_sum, experiment_segment_sum, \
                0, 0, \
                None, None,
            )
            continue

        # We add pseudocount of 1 to each total coverage of each segment
        mean_control = (control_segment_sum + 1) / (total_coverage_control + num_segments)
        mean_experiment = (experiment_segment_sum + 1) / (total_coverage_experiment + num_segments)
        detrended_profile = mean_experiment / mean_control

        yield (
            segment.start, segment.stop, coord, rel_coord, \
            control_segment_sum, experiment_segment_sum, \
            mean_control, mean_experiment, \
            detrended_profile, log2(detrended_profile),
        )
