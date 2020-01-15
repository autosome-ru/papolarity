import dataclasses
from .dataclass_tsv_serializable import DataclassTsvSerializable
from ..profile_comparison import profile_difference, slope_by_profiles

@dataclasses.dataclass(frozen=True)
class CoverageComparisonStats(DataclassTsvSerializable):
    transcript_id: str
    slope: float
    multipoint_slope: float
    profile_difference: float

    @classmethod
    def make_from_profiles(cls, transcript_id, profile_control, profile_experiment, segmentation):
        info = {
            'transcript_id': transcript_id,
            'slope': slope_by_profiles(profile_control, profile_experiment, segmentation, mode='center'),
            'multipoint_slope': slope_by_profiles(profile_control, profile_experiment, segmentation, mode='every_point'),
            'profile_difference': profile_difference(profile_control, profile_experiment, segmentation),
        }
        return cls(**info)
