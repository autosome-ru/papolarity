import pybedtools
import itertools
import sys
import numpy as np
import coreutils_sort

def transcript_coverages_in_file(alignment, cds_info_by_transcript, sort_transcripts=False):
    if sort_transcripts:
        bedgraph = alignment.genome_coverage(bg=True).coreutils_sort(key=['1,1'], stream=True)
    else:
        bedgraph = alignment.genome_coverage(bg=True, stream=True)
    for (transcript_id, bedgraph_intervals_iter) in itertools.groupby(bedgraph, lambda interval: interval.chrom):
        if transcript_id not in cds_info_by_transcript:
            # We don't know length of such transcripts
            print(f'Warning: there are reads from unknown transcript `{transcript_id}` in `{alignment.fn}`', file=sys.stderr)
            continue
        transcript_info = cds_info_by_transcript[transcript_id]
        profile = profile_by_bedgraph(bedgraph_intervals_iter, transcript_info.transcript_length)
        yield (transcript_info, profile)

def profile_by_bedgraph(bedgraph_intervals_iter, transcript_length):
    # Check: is it reasonable to use real-values or integer-valued coverages?
    profile = np.zeros(transcript_length, dtype=np.int)
    for interval in bedgraph_intervals_iter:
        # Note! pybedtools use 0-based coordinates (when in integer representation).
        # See https://daler.github.io/pybedtools/3-brief-examples.html and https://daler.github.io/pybedtools/intervals.html#zero-based-coords
        coverage = int(interval[3])
        profile[interval.start : interval.stop] = coverage
    return profile
