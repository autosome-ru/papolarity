import dataclasses
from .dataclass_tsv_serializable import DataclassTsvSerializable
from ..profile_comparison import l1_distance_by_segment_counts, slope_by_segment_counts, segmentwise_sums

@dataclasses.dataclass(frozen=True)
class CoverageComparisonStats(DataclassTsvSerializable):
    transcript_id: str
    slope: float
    slopelog: float
    l1_distance: float

    @classmethod
    def make_from_profiles(cls, transcript_id, control_profile, experiment_profile, segmentation):
        assert len(control_profile) == len(experiment_profile) == segmentation.segmentation_length
        control_sums = segmentwise_sums(segmentation, control_profile)
        experiment_sums = segmentwise_sums(segmentation, experiment_profile)

        slopes = slope_by_segment_counts(control_sums, experiment_sums, segmentation, log_mode=False)

        # slope of logarithms
        slopelogs = slope_by_segment_counts(control_sums, experiment_sums, segmentation, log_mode=True)
        info = {
            'transcript_id': transcript_id,
            'slope': slopes['center'],
            'slopelog': slopelogs['center'],
            'l1_distance': l1_distance_by_segment_counts(control_sums, experiment_sums),
        }
        return cls(**info)
