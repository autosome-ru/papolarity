import sys
from os.path import dirname
sys.path.insert(0, dirname(dirname(__file__)))
import argparse
import itertools
from gzip_utils import open_for_write
from dto.interval import Interval
from dto.coding_transcript_info import CodingTranscriptInfo
import dataclasses

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
        for (contig_name, segments) in itertools.groupby(bed_stream, lambda segment: segment.chrom):
            if contig_name in cds_info_by_transcript:
                cds_info = cds_info_by_transcript[contig_name]
                yield from self.segments_clipped_to_cds(segments, cds_info, contig_name)
            else:
                if allow_non_matching:
                    yield from segments

def get_argparser():
    argparser = argparse.ArgumentParser(
        prog = "clip_cds",
        description = "Clip any bed file in transcriptomic coordinates to CDS-region",
    )
    argparser.add_argument('cds_annotation', metavar='cds_annotation.tsv', help='CDS annotation') # 'gencode.vM22.cds_features.tsv'
    argparser.add_argument('bedfile', metavar='bedfile.bed', help = 'Coverage or segmentation file in bed format (3 standard columns + any number of non-standard)')

    # start and stop codon can have piles of reads, so we usually want to drop them
    # flank lengths to drop are of 15nt ~= half-ribosome (~half of riboseq footprint length)
    argparser.add_argument('--drop-5-flank', metavar='N', type=int, default=0, help="Clip N additional nucleotides from transcript start (5'-end)")
    argparser.add_argument('--drop-3-flank', metavar='N', type=int, default=0, help="Clip N additional nucleotides from transcript end (3'-end)")

    argparser.add_argument('--output-file', '-o', dest='output_file', help="Store results at this path")
    argparser.add_argument('--allow-non-matching', action='store_true', help="Allow transcripts which are not present in CDS-annotation (they are not clipped)")
    argparser.add_argument('--contig-naming', dest='contig_naming_mode', choices=['original', 'window'], default='window', help="Use original (chr1) or modified (chr1:23-45) contig name for resulting intervals")
    return argparser

argparser = get_argparser()
args = argparser.parse_args()

if (args.allow_non_matching) and (args.contig_naming_mode != 'original'):
    print('Attention! When `--allow-non-matching` is set, only `--contig-naming original` will give consistent contig names', file=sys.stderr)

cds_info_by_transcript = CodingTranscriptInfo.load_transcript_cds_info(args.cds_annotation)
bed_stream = Interval.each_in_file(args.bedfile)
clipper = Clipper(contig_naming_mode=args.contig_naming_mode,
                  drop_5_flank=args.drop_5_flank,
                  drop_3_flank=args.drop_3_flank)
with open_for_write(args.output_file) as output_stream:
    for interval in clipper.bedfile_clipped_to_cds(bed_stream, cds_info_by_transcript, allow_non_matching=args.allow_non_matching):
        print(interval, file=output_stream)
