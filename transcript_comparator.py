from coverage_profile import transcript_coverages_from_alignment
from utils import common_subsequence
from segmentation import make_joint_segmentation
from dto.coverage_comparison_stats import CoverageComparisonStats

class TranscriptComparator:
    # start and stop codon can have piles of reads, so we usually want to drop them
    # flank lengths to drop are of 15nt ~= half-ribosome (~half of riboseq footprint length)
    def __init__(self, cds_info_by_transcript, splitter, drop_start_flank=15, drop_stop_flank=15):
        self.cds_info_by_transcript = cds_info_by_transcript
        self.splitter = splitter
        self.drop_start_flank = drop_start_flank
        self.drop_stop_flank = drop_stop_flank

    def compare_coverage_streams(self, coverages_control, coverages_experiment):
        get_transcript_id = lambda transcript_coverage: transcript_coverage.transcript_id
        common_transcripts = common_subsequence([coverages_control, coverages_experiment], key=get_transcript_id)
        for (transcript_id, (coverage_control, coverage_experiment)) in common_transcripts:
            if transcript_id not in self.cds_info_by_transcript:
                continue
            cds_info = self.cds_info_by_transcript[transcript_id]
            info = self.compare_profiles(cds_info, coverage_control.coverage, coverage_experiment.coverage)
            if info:
                yield info

    def clipped_cds_profile(self, cds_info, coverage):
        profile = cds_info.cds_profile(coverage)
        # start and stop codon can have piles of reads, so we usually want to drop them
        clipped_profile = profile[self.drop_start_flank : -self.drop_stop_flank]
        return clipped_profile

    def compare_profiles(self, cds_info, coverage_control, coverage_experiment):
        cds_profile_control = self.clipped_cds_profile(cds_info, coverage_control)
        cds_profile_experiment = self.clipped_cds_profile(cds_info, coverage_experiment)
        segments = make_joint_segmentation([cds_profile_control, cds_profile_experiment], self.splitter)
        if segments == None:
            return None
        # stable_cds_profile_control = stabilize_profile(cds_profile_control, segments)
        # stable_cds_profile_experiment = stabilize_profile(cds_profile_experiment, segments)
        return CoverageComparisonStats.make_from_profiles(cds_info, cds_profile_control, cds_profile_experiment, segments)
