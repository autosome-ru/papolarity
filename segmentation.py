from utils import pool_profiles
import itertools
import pasio
import numpy as np
import dataclasses
from dto.interval import Interval
from typing import List

def make_joint_segmentation(coverages, splitter):
    pooled_coverage = pool_profiles(coverages)
    if len(pooled_coverage) == 0:
        return None
    return list(pasio.segments_with_scores(pooled_coverage, splitter))

def stabilize_profile(profile, segments):
    stable_profile = np.zeros_like(profile)
    for segment in segments:
        start = segment.start
        stop = segment.stop
        stable_profile[start:stop] = np.mean(profile[start:stop])
    return stable_profile

@dataclasses.dataclass(frozen=True)
class Segmentation:
    chrom: str
    segments: List[Interval] # bed-coordinates [a, b)

    def stabilize_profile(self, profile):
        stable_profile = np.zeros_like(profile)
        for segment in self.segments:
            start = segment.start
            stop = segment.stop
            stable_profile[start:stop] = np.mean(profile[start:stop])
        return stable_profile

    @classmethod
    def each_in_file(cls, filename, force_gzip=None, header=False):
        segment_stream = Interval.each_in_file(filename, header=header, force_gzip=force_gzip)
        for (chrom, segments_iter) in itertools.groupby(segment_stream, key=lambda segment: segment.chrom):
            segments = list(segments_iter)
            yield cls(chrom, segments)
