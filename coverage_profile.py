from pybedtools import BedTool
from decimal import Decimal
import itertools
import sys
import numpy as np
import coreutils_sort

def transcript_coverages_from_file(alignment_fn, sort_transcripts=False, dtype=float):
    alignment = BedTool(alignment_fn)
    yield from transcript_coverages_from_alignment(alignment, sort_transcripts=sort_transcripts, dtype=dtype)

def transcript_coverages_from_alignment(alignment, sort_transcripts=False, dtype=float):
    bedgraph = make_coverage(alignment, sort_transcripts=sort_transcripts)
    yield from transcript_coverages_from_bedgraph(bedgraph, dtype=dtype)

def make_coverage(alignment, sort_transcripts=False):
    if sort_transcripts:
        # sort by chromosome (transcript_id in case of transcriptomic alignments)
        bedgraph = alignment.genome_coverage(bga=True).coreutils_sort(key=['1,1'], stream=True)
    else:
        bedgraph = alignment.genome_coverage(bga=True, stream=True)
    return bedgraph

def transcript_coverages_from_bedgraph(bedgraph, dtype=float):
    for (transcript_id, bedgraph_intervals_iter) in itertools.groupby(bedgraph, lambda interval: interval.chrom):
        bedgraph_intervals = list(bedgraph_intervals_iter)
        transcript_length = max(interval.stop for interval in bedgraph_intervals)
        profile = profile_by_bedgraph(bedgraph_intervals, transcript_length, dtype=dtype)
        yield (transcript_id, profile)

def profile_by_bedgraph(bedgraph_intervals_iter, transcript_length, dtype=float):
    profile = np.zeros(transcript_length, dtype=dtype)
    for interval in bedgraph_intervals_iter:
        # Note! pybedtools use 0-based coordinates (when in integer representation).
        # See https://daler.github.io/pybedtools/3-brief-examples.html and https://daler.github.io/pybedtools/intervals.html#zero-based-coords
        try:
            coverage = dtype(interval[3])
        except:
            # if dtype is integer, string "1.2345e6" can't be directly converted to int
            coverage = dtype(Decimal(interval[3]))
        profile[interval.start : interval.stop] = coverage
    return profile
