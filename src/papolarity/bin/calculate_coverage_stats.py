import argparse
from ..utils import common_subsequence
from ..dto.transcript_coverage import TranscriptCoverage
from ..segmentation import Segmentation
from ..dto.coding_transcript_info import CodingTranscriptInfo
from ..dto.coverage_comparison_stats import CoverageComparisonStats

def compare_coverage_streams(cds_info_by_transcript, segmentations, control_coverages, experiment_coverages):
    segmentation_and_coverages = [
        segmentations,
        control_coverages,
        experiment_coverages,
    ]

    key_extractors = [
        lambda segment: segment.chrom,
        lambda transcript_coverage: transcript_coverage.transcript_id,
        lambda transcript_coverage: transcript_coverage.transcript_id,
    ]

    transcript_stream = common_subsequence(segmentation_and_coverages, key=key_extractors, check_sorted=True)
    for (transcript_id, (segmentation, control_coverage, experiment_coverage)) in transcript_stream:
        if transcript_id not in cds_info_by_transcript:
            continue
        cds_info = cds_info_by_transcript[transcript_id]
        yield CoverageComparisonStats.make_from_profiles(cds_info, control_coverage.coverage, experiment_coverage.coverage, segmentation)

def configure_argparser(argparser=None):
    if not argparser:
        argparser = argparse.ArgumentParser(
            prog = "calculate_coverage_stats",
            description = "Coverage profile comparison",
        )
    argparser.add_argument('cds_annotation', metavar='cds_annotation.tsv', help='CDS annotation') # 'gencode.vM22.cds_features.tsv'
    argparser.add_argument('segmentation', metavar='segmentation.bed', help='Segmentation')
    argparser.add_argument('coverage_control', metavar='control.bedgraph', help='Coverage for control data')
    argparser.add_argument('coverage_experiment', metavar='experiment.bedgraph', help='Coverage for experiment data')
    return argparser

def main():
    argparser = configure_argparser()
    args = argparser.parse_args()
    invoke(args)

def invoke(args):
    segmentation_stream = Segmentation.each_in_file(args.segmentation, header=False)

    cds_info_by_transcript = CodingTranscriptInfo.load_transcript_cds_info(args.cds_annotation)
    control_coverages = TranscriptCoverage.each_in_file(args.coverage_control, header=False, dtype=int)
    experiment_coverages = TranscriptCoverage.each_in_file(args.coverage_experiment, header=False, dtype=int)

    transcript_comparison_infos = compare_coverage_streams(cds_info_by_transcript, segmentation_stream, control_coverages, experiment_coverages)
    transcript_comparison_infos = list(transcript_comparison_infos)
    CoverageComparisonStats.print(transcript_comparison_infos, extended=True)
