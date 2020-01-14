import sys
import dataclasses
import numpy as np
from .dataclass_tsv_serializable import DataclassTsvSerializable
from ..profile_comparison import profile_difference, slope_by_profiles

@dataclasses.dataclass(frozen=True)
class CoverageComparisonStats(DataclassTsvSerializable):
    transcript_id: str
    slope: float
    control_mean_coverage: float
    experiment_mean_coverage: float
    control_total_coverage: float
    experiment_total_coverage: float
    num_segments: int
    multipoint_slope: float
    profile_difference: float

    @classmethod
    def make_from_profiles(cls, transcript_id, cds_profile_control, cds_profile_experiment, segmentation):
        info = {
            'transcript_id': transcript_id,
            'slope': slope_by_profiles(cds_profile_control, cds_profile_experiment, segmentation, mode='center'),
            'control_mean_coverage': np.mean(cds_profile_control), 'experiment_mean_coverage': np.mean(cds_profile_experiment),
            'control_total_coverage': np.sum(cds_profile_control), 'experiment_total_coverage': np.sum(cds_profile_experiment),
            'num_segments': segmentation.num_segments,
            'multipoint_slope': slope_by_profiles(cds_profile_control, cds_profile_experiment, segmentation, mode='every_point'),
            'profile_difference': profile_difference(cds_profile_control, cds_profile_experiment, segmentation),
        }
        return cls(**info)

    @classmethod
    def print(cls, infos, file=sys.stdout, extended=False):
        if extended:
            print('\t'.join([cls.header(), 'geom_mean_coverage']), file=file)
            for info in infos:
                print('\t'.join(map(str, [info, info.geom_mean_coverage,])), file=file)
        else:
            print(cls.header(), file=file)
            for info in infos:
                print(info, file=file)

    @property
    def geom_mean_coverage(self):
        return (self.control_mean_coverage * self.experiment_mean_coverage) ** 0.5
