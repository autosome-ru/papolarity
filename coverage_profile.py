from decimal import Decimal
import itertools
import numpy as np
import coreutils_sort
from dto.transcript_coverage import TranscriptCoverage
from dto.coverage_interval import CoverageInterval

def transcript_coverages_from_alignment(alignment, sort_transcripts=False, dtype=float):
    bedgraph = make_coverage(alignment, sort_transcripts=sort_transcripts, dtype=dtype)
    coverage_intervals = coverage_intervals_from_bedgraph(bedgraph, dtype=dtype)
    yield from transcript_coverages_from_bedgraph(coverage_intervals, dtype=dtype)

def make_coverage(alignment, sort_transcripts=False, stream=True, dtype=float):
    if sort_transcripts:
        # sort by chromosome (transcript_id in case of transcriptomic alignments)
        bedgraph = alignment.genome_coverage(bga=True).coreutils_sort(key=['1,1', '2,2n'], stream=stream)
    else:
        bedgraph = alignment.genome_coverage(bga=True, stream=stream)
    return bedgraph

def coverage_intervals_from_bedgraph(bedgraph, dtype=float):
    for interval in bedgraph:
        yield CoverageInterval(*interval[0:4], dtype=dtype)

def transcript_coverages_from_bedgraph(bedgraph, dtype=float):
    for (transcript_id, bedgraph_intervals_iter) in itertools.groupby(bedgraph, lambda interval: interval.chrom):
        bedgraph_intervals = list(bedgraph_intervals_iter)
        transcript_length = max(interval.stop for interval in bedgraph_intervals)
        profile = np.zeros(transcript_length, dtype=dtype)
        for interval in bedgraph_intervals:
            # Note! pybedtools use 0-based coordinates (when in integer representation).
            # See https://daler.github.io/pybedtools/3-brief-examples.html and https://daler.github.io/pybedtools/intervals.html#zero-based-coords
            profile[interval.start : interval.stop] = interval.coverage
        yield TranscriptCoverage(transcript_id, profile)
