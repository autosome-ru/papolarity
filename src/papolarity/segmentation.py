import itertools
import dataclasses
from typing import List
import numpy as np
from .dto.interval import Interval

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

    @property
    def segmentation_length(self):
        return self.segments[-1].stop

    @property
    def num_segments(self):
        return len(self.segments)

    def stabilize_profile(self, profile):
        stable_profile = np.zeros_like(profile)
        for segment in self.segments:
            start = segment.start
            stop = segment.stop
            stable_profile[start:stop] = np.mean(profile[start:stop])
        return stable_profile

    # clip flanks so that region is in specified window
    def clip_to_window(self, start, stop):
        return self.clip_flanks(start, self.segmentation_length - stop)

    # clip flanks of specified length
    def clip_flanks(self, flank_5, flank_3):
        clip_len = self.segmentation_length - flank_3 - flank_5
        clipped_segments = []
        for segment in self.segments:
            segment_start = max(0, segment.start - flank_5)
            segment_stop = min(clip_len, segment.stop - flank_5)
            if (segment_stop <= 0) or (clip_len <= segment_start):
                continue
            clipped_segment = Interval(self.chrom, segment_start, segment_stop)
            clipped_segments.append(clipped_segment)
        return Segmentation(self.chrom, clipped_segments)

    @classmethod
    def each_in_file(cls, filename, force_gzip=None, header=False):
        segment_stream = Interval.each_in_file(filename, header=header, force_gzip=force_gzip)
        for (chrom, segments_iter) in itertools.groupby(segment_stream, key=lambda segment: segment.chrom):
            segments = list(segments_iter)
            yield cls(chrom, segments)
