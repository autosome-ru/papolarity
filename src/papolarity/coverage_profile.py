from . import coreutils_sort
from .dto.transcript_coverage import TranscriptCoverage
from .dto.coverage_interval import CoverageInterval

def transcript_coverages_from_alignment(alignment, sort_transcripts=False, dtype=float):
    bedgraph = make_coverage(alignment, sort_transcripts=sort_transcripts, dtype=dtype)
    coverage_intervals = coverage_intervals_from_bedgraph(bedgraph, dtype=dtype)
    yield from TranscriptCoverage.each_in_bedgraph(coverage_intervals, dtype=dtype)

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
