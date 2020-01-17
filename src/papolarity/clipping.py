import dataclasses
from itertools import groupby
from .dto.interval import Interval

def segments_clipped_to_window(segments, window_start, window_stop, contig_name):
    for segment in segments:
        segment_start = max(0, segment.start - window_start)
        segment_stop = min(window_stop, segment.stop) - window_start
        if segment_stop - segment_start > 0:
            yield Interval(contig_name, segment_start, segment_stop, segment.rest)

@dataclasses.dataclass
class Clipper:
    contig_naming_mode: str = 'window'
    drop_5_flank: int = 0
    drop_3_flank: int = 0

    def __post_init__(self):
        if self.contig_naming_mode not in ['original', 'window']:
            raise ValueError(f'Unknown contig naming mode `{self.contig_naming_mode}`')

    def customize_contig_name(self, contig_name, window_start, window_stop):
        if self.contig_naming_mode == 'original':
            return contig_name
        elif self.contig_naming_mode == 'window':
            return f'{contig_name}:{window_start}-{window_stop}'
        else:
            raise ValueError(f'Unknown contig naming mode `{self.contig_naming_mode}`')

    def segments_clipped_to_cds(self, segments, cds_info, contig_name):
        window_start = cds_info.cds_start + self.drop_5_flank
        window_stop = cds_info.cds_stop - self.drop_3_flank
        customized_contig_name = self.customize_contig_name(contig_name, window_start, window_stop)
        yield from segments_clipped_to_window(segments, window_start, window_stop, customized_contig_name)

    def bedfile_clipped_to_cds(self, bed_stream, cds_info_by_transcript, allow_non_matching=False):
        '''
        Attention! if `contig_naming_mode` is not "original" and `allow_non_matching` is True,
        then contig names will be in inconsistent formats:
        transcripts with have matching cds_info, will have customized contig names,
        while transcripts without matching cds_info will have original contig names.
        '''
        for (contig_name, segments) in groupby(bed_stream, lambda segment: segment.chrom):
            if contig_name in cds_info_by_transcript:
                cds_info = cds_info_by_transcript[contig_name]
                if cds_info.is_coding:
                    yield from self.segments_clipped_to_cds(segments, cds_info, contig_name)
                elif allow_non_matching:
                    yield from segments
            else:
                if allow_non_matching:
                    yield from segments
