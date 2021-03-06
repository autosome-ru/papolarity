from . import coreutils_sort
from .dto.transcript_coverage import TranscriptCoverage
from .dto.coverage_interval import CoverageInterval

def make_coverage(alignment, sort_transcripts='no', stream=True, dtype=float):
    if sort_transcripts == 'case-sensitive':
        # sort by chromosome (transcript_id in case of transcriptomic alignments)
        opts = {'key': ['1,1', '2,2n'], 'ignore-case': False}
        bedgraph = alignment.genome_coverage(bga=True).coreutils_sort(**opts, stream=stream)
    elif sort_transcripts == 'case-insensitive':
        opts = {'key': ['1,1', '2,2n'], 'ignore-case': True}
        bedgraph = alignment.genome_coverage(bga=True).coreutils_sort(**opts, stream=stream)
    elif sort_transcripts == 'no':
        bedgraph = alignment.genome_coverage(bga=True, stream=stream)
    else:
        raise ValueError('sort_transcripts should be one of no/case-sensitive/case-insensitive.')
    return bedgraph

def coverage_intervals_from_bedgraph(bedgraph, dtype=float):
    for interval in bedgraph:
        yield CoverageInterval(*interval[0:4], dtype=dtype)
