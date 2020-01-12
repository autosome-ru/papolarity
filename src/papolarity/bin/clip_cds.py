import sys
import argparse
from ..gzip_utils import open_for_write
from ..dto.interval import Interval
from ..dto.coding_transcript_info import CodingTranscriptInfo
from ..clipping import Clipper

def configure_argparser(argparser=None):
    if not argparser:
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

def main():
    argparser = configure_argparser()
    args = argparser.parse_args()
    invoke(args)

def invoke(args):
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
