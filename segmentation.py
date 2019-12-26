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

    def __post_init__(self):
        object.__setattr__(self, 'segments', sorted(self.segments)) # assign frozen attribute.
        if any(self.chrom != segment.chrom for segment in self.segments):
            raise ValueError(f"Segmentation intervals all should be on specified chromosome/contig `{self.chrom}`")
        if len(self.segments) > 0 and self.segments[0].start != 0:
            raise ValueError(f"Segmentation should cover entire interval (but first segment {self.segments[0]} doesn't start from zero)")
        if any(self.segments[i].stop != self.segments[i + 1].start for i in range(len(self.segments) - 1)):
            raise ValueError(f"Segments should be contigious but weren't: {self.segments}")

    def __len__(self):
        return self.segments[-1].stop

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
