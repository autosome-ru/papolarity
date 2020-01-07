import sys
from os.path import dirname
sys.path.insert(0, dirname(dirname(__file__)))
import argparse
from pybedtools import BedTool
from dto.coding_transcript_info import CodingTranscriptInfo
from dto.coverage_comparison_stats import CoverageComparisonStats
from transcript_comparator import TranscriptComparator
import pasio
import logging

def get_argparser():
    argparser = argparse.ArgumentParser(
        prog = "calculate_coverage_stats",
        description = "Coverage profile comparison",
    )
    argparser.add_argument('cds_annotation', help='CDS annotation') # 'gencode.vM22.cds_features.tsv'
    argparser.add_argument('alignment_control', help='Alignment for control data')
    argparser.add_argument('alignment_experiment', help='Alignment for experiment data')
    return argparser

logger = logging.getLogger('pasio')
logger.setLevel(logging.WARNING)
splitter = pasio.configure_splitter(alpha=1, beta=1, algorithm='rounds', window_size=2500, window_shift=1250, num_rounds=None, no_split_constant=True)

argparser = get_argparser()
args = argparser.parse_args()

cds_info_by_transcript = CodingTranscriptInfo.load_transcript_cds_info(args.cds_annotation)
alignment_control = BedTool(args.alignment_control)
alignment_experiment = BedTool(args.alignment_experiment)

comparator = TranscriptComparator(cds_info_by_transcript, splitter, drop_start_flank=15, drop_stop_flank=15)
transcript_comparison_infos = comparator.compare_multiple_alignments(alignment_control, alignment_experiment)
transcript_comparison_infos = list(transcript_comparison_infos)
CoverageComparisonStats.print(transcript_comparison_infos, extended=True)
