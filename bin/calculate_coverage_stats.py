import sys
from pybedtools import BedTool
from dto.coding_transcript_info import CodingTranscriptInfo
from dto.coverage_comparison_stats import CoverageComparisonStats
from transcript_comparator import TranscriptComparator

import pasio
import logging

logger = logging.getLogger('pasio')
logger.setLevel(logging.WARNING)
splitter = pasio.configure_splitter(alpha=1, beta=1, algorithm='rounds', window_size=2500, window_shift=1250, num_rounds=None, no_split_constant=True)

cds_annotation_fn = sys.argv[1] # 'gencode.vM22.cds_features.tsv'
cds_info_by_transcript = CodingTranscriptInfo.load_transcript_cds_info(cds_annotation_fn)

alignment_control_fn = sys.argv[2]
alignment_experiment_fn = sys.argv[3]

alignment_control = BedTool(alignment_control_fn)
alignment_experiment = BedTool(alignment_experiment_fn)

comparator = TranscriptComparator(cds_info_by_transcript, splitter, drop_start_flank=15, drop_stop_flank=15)
transcript_comparison_infos = comparator.compare_multiple_alignments(alignment_control, alignment_experiment)
transcript_comparison_infos = list(transcript_comparison_infos)
CoverageComparisonStats.print(transcript_comparison_infos, extended=True)