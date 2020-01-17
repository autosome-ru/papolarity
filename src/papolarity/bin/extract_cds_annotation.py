import argparse
from ..gzip_utils import open_for_write
from ..annotation import Annotation
from ..dto.coding_transcript_info import CodingTranscriptInfo

# gene:attr=value,value,value
# transcript:attr=value,value,value
# *:attr=value,value,value
def parse_condition(condition_str):
    feature_types, attr_condition_str = condition_str.split(':', maxsplit=1)
    k, vs = attr_condition_str.split('=', maxsplit=1)
    return {'attr': k, 'acceptable_values': vs.split(','), 'feature_types': feature_types.split(',')}

# Known issue: filters treat all attribute values as strings
def create_attr_filter(condition_config):
    feature_types = condition_config['feature_types']
    k = condition_config['attr']
    vs = condition_config['acceptable_values']
    if '*' in feature_types:
        return lambda rec: str(rec.attributes.get(k, '')) in vs
    else:
        return lambda rec: ((str(rec.attributes.get(k, '')) in vs) and (rec.type in feature_types)) or (rec.type not in feature_types)

def configure_argparser(argparser=None):
    if not argparser:
        argparser = argparse.ArgumentParser(prog="extract_cds_annotation", description = "Extract CDS annotation in transcriptomic coordinates from genomic annotation")
    argparser.add_argument('gtf_annotation', metavar='annotation.gtf', help='Genomic annotation in GTF-format')
    argparser.add_argument('--output-file', '-o', dest='output_file', help="Store results at this path")
    argparser.add_argument('--attr-filter', action='append', dest='attr_filters', default=[],
                                            help="Filter records so that attributes has one of specified values.\n"
                                                 "Filter is applied only to records of specified feature type (`*` for all feature types).\n"
                                                 "Format: `type1,type2:tag=value1,value2,...`")
    # Known issue: Right now all multivalue keys are ignored, so filters won't work with them
    # ToDo: add ability to specify multivalue keys and add filters like `--filter-has tag value`
    return argparser

def main():
    argparser = configure_argparser(argparser)
    args = argparser.parse_args()
    invoke(args)

def invoke(args):
    condition_configs = [parse_condition(condition_str) for condition_str in args.attr_filters]
    filters = [create_attr_filter(condition_config) for condition_config in condition_configs]
    filter_conjuction = lambda rec: all(f(rec) for f in filters)
    relevant_attributes = {config['attr'] for config in condition_configs}

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
