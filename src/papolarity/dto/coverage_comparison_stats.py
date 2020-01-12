import sys
import dataclasses
from itertools import groupby
import numpy as np
from .dataclass_tsv_serializable import DataclassTsvSerializable
from .coding_transcript_info import CodingTranscriptInfo
from ..polarity_score import polarity_score
from ..profile_comparison import profile_difference, slope_by_profiles

@dataclasses.dataclass(frozen=True)
class CoverageComparisonStats(DataclassTsvSerializable):
    gene_id: str
    transcript_id: str
    transcript_length: int
    cds_start: int
    cds_stop: int
    slope: float
    control_mean_coverage: float
    experiment_mean_coverage: float
    control_total_coverage: float
    experiment_total_coverage: float
    control_polarity_score: float
    experiment_polarity_score: float
    num_segments: int
    multipoint_slope: float
    profile_difference: float

    @property
    def cds_info(self):
        return CodingTranscriptInfo(self.gene_id, self.transcript_id, self.transcript_length, self.cds_start, self.cds_stop)

    @classmethod
    def make_from_profiles(cls, cds_info, cds_profile_control, cds_profile_experiment, segmentation):
        info = {
            **dataclasses.asdict(cds_info),
            'slope': slope_by_profiles(cds_profile_control, cds_profile_experiment, segmentation, mode='center'),
            'control_mean_coverage': np.mean(cds_profile_control), 'experiment_mean_coverage': np.mean(cds_profile_experiment),
            'control_total_coverage': np.sum(cds_profile_control), 'experiment_total_coverage': np.sum(cds_profile_experiment),
            'control_polarity_score': polarity_score(cds_profile_control), 'experiment_polarity_score': polarity_score(cds_profile_experiment),
            'num_segments': segmentation.num_segments,
            'multipoint_slope': slope_by_profiles(cds_profile_control, cds_profile_experiment, segmentation, mode='every_point'),
            'profile_difference': profile_difference(cds_profile_control, cds_profile_experiment, segmentation),
        }
        return cls(**info)

    @classmethod
    def print(cls, infos, file=sys.stdout, extended=False):
        if extended:
            infos = list(infos)
            best_transcript_ids = set(info.transcript_id for info in cls.choose_best_transcript(infos))

            print('\t'.join([cls.header(), 'geom_mean_coverage', 'polarity_difference', 'is_best_transcript']), file=file)
            for info in infos:
                if info.transcript_id in best_transcript_ids:
                    is_best_transcript = '+'
                else:
                    is_best_transcript = '-'
                print('\t'.join(map(str, [info, info.geom_mean_coverage(), info.polarity_difference(), is_best_transcript])), file=file)
        else:
            print(cls.header(), file=file)
            for info in infos:
                print(info, file=file)

    def geom_mean_coverage(self):
        return (self.control_mean_coverage * self.experiment_mean_coverage) ** 0.5

    def polarity_difference(self):
        return self.experiment_polarity_score - self.control_polarity_score

    # by default takes max by geometric mean of cds-mean between experiments
    @classmethod
    def choose_best_transcript(cls, slope_data, key=lambda info: info.geom_mean_coverage()):
        by_gene = lambda info: info.gene_id
        slope_data = sorted(slope_data, key=by_gene)
        for (gene_id, infos) in groupby(slope_data, by_gene):
            best_transcript_info = max(infos, key=key)
            yield best_transcript_info
