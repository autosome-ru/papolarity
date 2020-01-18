import argparse
from ..gzip_utils import open_for_write
from ..annotation import Annotation
from ..annotation_filter import parse_condition, create_record_filter

def clip_sequence(sequence, drop_5_flank, drop_3_flank):
    return sequence[drop_5_flank : (len(sequence) - drop_3_flank)]

def configure_argparser(argparser=None):
    if not argparser:
        argparser = argparse.ArgumentParser(
            prog = "cds_sequence",
            description = "Extract CDS sequences from GTF annotation and genome assembly",
        )
    argparser.add_argument('gtf_annotation', metavar='annotation.gtf', help='Genomic annotation in GTF format')
    argparser.add_argument('assembly', metavar='assembly.fa', help='Genome assembly in FASTA format')
    argparser.add_argument('--region-type', choices=['cds', 'cds_with_stop', 'exons', 'full'], default='cds',
                           help='Genome assembly in FASTA format')
    argparser.add_argument('--drop-5-flank', metavar='N', type=int, default=0, help="Clip N additional nucleotides from transcript start (5'-end)")
    argparser.add_argument('--drop-3-flank', metavar='N', type=int, default=0, help="Clip N additional nucleotides from transcript end (3'-end)")

    argparser.add_argument('--output-file', '-o', dest='output_file', help="Store results at this path")
    argparser.add_argument('--attr-filter', action='append', dest='filters', default=[], 
                                            help="Filter records so that attributes has one of specified values.\n"
                                                 "Format: `attribute=value_1,value_2,...`")
    # Known issue: Right now all multivalue keys are ignored, so filters won't work with them
    # ToDo: add ability to specify multivalue keys and add filters like `--filter-has tag value`
    return argparser

def main():
    argparser = configure_argparser()
    args = argparser.parse_args()
    invoke(args)

def invoke(args):
    condition_configs = [parse_condition(condition_str) for condition_str in args.filters]
    filters = [create_record_filter(condition_config) for condition_config in condition_configs]
    filter_conjuction = lambda rec: all(f(rec) for f in filters)
    relevant_attributes = {k for (k,v) in condition_configs}

    annotation = Annotation.load(
        args.gtf_annotation,
        relevant_attributes=relevant_attributes,
        multivalue_keys=set(),
        ignore_unknown_multivalues=True,
        condition=filter_conjuction,
    )
    transcript_ids_list = annotation.transcript_by_id.keys()

    with open_for_write(args.output_file) as output_stream:
        for transcript_id, sequence in annotation.transcript_sequences(transcript_ids_list, args.assembly, feature_type=args.region_type):
            clipped_sequence = clip_sequence(sequence, drop_5_flank=args.drop_5_flank, drop_3_flank=args.drop_3_flank)
            print(f'>{transcript_id}', file=output_stream)
            print(clipped_sequence, file=output_stream)
