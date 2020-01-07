import sys
from os.path import dirname
sys.path.insert(0, dirname(dirname(__file__)))
import argparse
from dto.transcript_coverage import TranscriptCoverage
from segmentation import Segmentation
from dto.coding_transcript_info import CodingTranscriptInfo
from dto.coverage_comparison_stats import CoverageComparisonStats
from transcript_comparator import TranscriptComparator

def get_argparser():
    argparser = argparse.ArgumentParser(
        prog = "calculate_coverage_stats",
        description = "Coverage profile comparison",
    )
    argparser.add_argument('cds_annotation', metavar='cds_annotation.tsv', help='CDS annotation') # 'gencode.vM22.cds_features.tsv'
    argparser.add_argument('segmentation', metavar='segmentation.bed', help='Segmentation')
    argparser.add_argument('coverage_control', metavar='control.bedgraph', help='Coverage for control data')
    argparser.add_argument('coverage_experiment', metavar='experiment.bedgraph', help='Coverage for experiment data')
    return argparser

argparser = get_argparser()
args = argparser.parse_args()

segmentation_stream = Segmentation.each_in_file(args.segmentation, header=False)

cds_info_by_transcript = CodingTranscriptInfo.load_transcript_cds_info(args.cds_annotation)
control_coverages = TranscriptCoverage.each_in_file(args.coverage_control, header=False, dtype=int)
experiment_coverages = TranscriptCoverage.each_in_file(args.coverage_experiment, header=False, dtype=int)

comparator = TranscriptComparator(cds_info_by_transcript, drop_start_flank=15, drop_stop_flank=15)
transcript_comparison_infos = comparator.compare_coverage_streams(segmentation_stream, control_coverages, experiment_coverages)
transcript_comparison_infos = list(transcript_comparison_infos)
CoverageComparisonStats.print(transcript_comparison_infos, extended=True)
