import dataclasses
import itertools
import numpy as np
from .coverage_interval import CoverageInterval

@dataclasses.dataclass
class TranscriptCoverage:
    transcript_id: str
    coverage: np.array

    @classmethod
    def each_in_file(cls, filename, header=False, dtype=float):
        bedgraph_stream = CoverageInterval.each_in_file(filename, header=False)
        yield from cls.each_in_bedgraph(bedgraph_stream, dtype=dtype)

    @classmethod
    def each_in_bedgraph(cls, bedgraph_stream, dtype=float):
        for (transcript_id, bedgraph_intervals_iter) in itertools.groupby(bedgraph_stream, lambda interval: interval.chrom):
            bedgraph_intervals = list(bedgraph_intervals_iter)
            transcript_length = max(interval.stop for interval in bedgraph_intervals)
            profile = np.zeros(transcript_length, dtype=dtype)
            for interval in bedgraph_intervals:
                # Note! pybedtools use 0-based coordinates (when in integer representation).
                # See https://daler.github.io/pybedtools/3-brief-examples.html and https://daler.github.io/pybedtools/intervals.html#zero-based-coords
                profile[interval.start : interval.stop] = interval.coverage
            yield cls(transcript_id, profile)
