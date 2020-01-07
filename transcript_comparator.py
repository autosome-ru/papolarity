from coverage_profile import transcript_coverages_from_alignment
from utils import common_subsequence
from segmentation import make_joint_segmentation
from dto.coverage_comparison_stats import CoverageComparisonStats

class TranscriptComparator:
    # start and stop codon can have piles of reads, so we usually want to drop them
    # flank lengths to drop are of 15nt ~= half-ribosome (~half of riboseq footprint length)
    def __init__(self, cds_info_by_transcript, drop_start_flank=15, drop_stop_flank=15):
        self.cds_info_by_transcript = cds_info_by_transcript
        self.drop_start_flank = drop_start_flank
        self.drop_stop_flank = drop_stop_flank

    def compare_coverage_streams(self, segmentations, control_coverages, experiment_coverages):
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
            if transcript_id not in self.cds_info_by_transcript:
                continue
            cds_info = self.cds_info_by_transcript[transcript_id]

            stable_control_coverage = segmentation.stabilize_profile(control_coverage.coverage)
            stable_experiment_coverage = segmentation.stabilize_profile(experiment_coverage.coverage)
            cds_profile_control = self.clipped_cds_profile(cds_info, stable_control_coverage)
            cds_profile_experiment = self.clipped_cds_profile(cds_info, stable_experiment_coverage)
            cds_segmentation = self.clipped_segmentation(cds_info, segmentation)

            if cds_segmentation.segmentation_length == 0:
                continue
            yield CoverageComparisonStats.make_from_profiles(cds_info, cds_profile_control, cds_profile_experiment, cds_segmentation)

    def clipped_segmentation(self, cds_info, segmentation):
        return segmentation.clip_to_window(cds_info.cds_start + self.drop_start_flank,
                                           cds_info.cds_stop - self.drop_stop_flank)

    def clipped_cds_profile(self, cds_info, coverage):
        profile = cds_info.cds_profile(coverage)
        # start and stop codon can have piles of reads, so we usually want to drop them
        clipped_profile = profile[self.drop_start_flank : -self.drop_stop_flank]
        return clipped_profile
