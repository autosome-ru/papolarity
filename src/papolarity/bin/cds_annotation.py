import argparse
from ..gzip_utils import open_for_write
from ..annotation import Annotation
from ..dto.coding_transcript_info import CodingTranscriptInfo
from ..annotation_filter import parse_condition, create_record_filter

def configure_argparser(argparser=None):
    if not argparser:
        argparser = argparse.ArgumentParser(prog="cds_annotation", description = "Extract CDS annotation in transcriptomic coordinates from genomic annotation")
    argparser.add_argument('gtf_annotation', metavar='annotation.gtf', help='Genomic annotation in GTF-format')
    argparser.add_argument('--output-file', '-o', dest='output_file', help="Store results at this path")
    argparser.add_argument('--attr-filter', action='append', dest='filters', default=[], 
                                            help="Filter records so that attributes has one of specified values.\n"
                                                 "Format: `attribute=value_1,value_2,...`")
    # Known issue: Right now all multivalue keys are ignored, so filters won't work with them
    # ToDo: add ability to specify multivalue keys and add filters like `--filter-has tag value`
    return argparser

def main():
    argparser = configure_argparser(argparser)
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
    with open_for_write(args.output_file) as output_stream:
        print(CodingTranscriptInfo.header(), file=output_stream)
        for transcript_id in annotation.transcript_by_id:
            cds_info = annotation.coding_transcript_info(transcript_id)
            print(cds_info, file=output_stream)
