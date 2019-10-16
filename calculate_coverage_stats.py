import sys
import pasio_wrapper
from pybedtools import BedTool
from annotation import load_transcript_cds_info
from slope import CoverageComparisonStats, TranscriptComparator

splitter = pasio_wrapper.pasio_splitter()

cds_annotation_fn = sys.argv[1] # 'gencode.vM22.cds_features.tsv'
cds_info_by_transcript = load_transcript_cds_info(cds_annotation_fn)

alignment_control_fn = sys.argv[2]
alignment_experiment_fn = sys.argv[3]

alignment_control = BedTool(alignment_control_fn)
alignment_experiment = BedTool(alignment_experiment_fn)

transcript_comparator = TranscriptComparator(cds_info_by_transcript, splitter, drop_start_flank=15, drop_stop_flank=15)
infos = transcript_comparator.compare_multiple_alignments(alignment_control, alignment_experiment)
CoverageComparisonStats.print(infos)
