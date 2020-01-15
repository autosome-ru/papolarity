import argparse
from .version import __version__
from .bin import get_coverage, pool_coverages, extract_cds_annotation, clip_cds, \
                 coverage_properties, join_tables, choose_best, \
                 compare_coverages, plot_slope_distributions, \
                 shrinkage, plot_slope_distributions, \
                 flatten_coverage

def configure_argparser(argparser=None):
    if not argparser:
        argparser = argparse.ArgumentParser(prog="papolarity", description = "Main entrypoint of papolarity")
    argparser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    subparsers = argparser.add_subparsers(metavar='subcommand')
    
    subparser_configs = [
        {'cmd': 'get_coverage', 'namespace': get_coverage, 'help': 'Generates coverage from an alignment'},
        {'cmd': 'pool_coverages', 'namespace': pool_coverages, 'help': 'Pool coverage profiles'},
        {'cmd': 'extract_cds_annotation', 'namespace': extract_cds_annotation, 'help': 'Extract CDS annotation in transcriptomic coordinates from genomic annotation'},
        {'cmd': 'clip_cds', 'namespace': clip_cds, 'help': 'Clip any bed file in transcriptomic coordinates to CDS-region'},
        {'cmd': 'coverage_properties', 'namespace': coverage_properties, 'help': 'Calculate coverage profile properties'},
        {'cmd': 'join_tables', 'namespace': join_tables, 'help': 'Join the second table to the first one'},
        {'cmd': 'choose_best', 'namespace': choose_best, 'help': 'Choose best element from each group (e.g. best transcipt for each gene)'},
        {'cmd': 'compare_coverages', 'namespace': compare_coverages, 'help': 'Coverage profile comparison'},
        {'cmd': 'plot_slope_distributions', 'namespace': plot_slope_distributions, 'help': 'Plot distributions of polarities and slopes over genes'},
        {'cmd': 'shrinkage', 'namespace': shrinkage, 'help': 'Make length-dependend shrinkage of slope properties'},
        {'cmd': 'flatten_coverage', 'namespace': flatten_coverage, 'help': 'Flatten coverage profiles by averaging data through given segments'},
    ]
    for subparser_config in subparser_configs:
        invocation_fn = subparser_config['namespace'].invoke
        configurator = subparser_config['namespace'].configure_argparser
        subparser = subparsers.add_parser(subparser_config['cmd'], help=subparser_config['help'])
        subparser.set_defaults(invocation_fn=invocation_fn)
        configurator(subparser)
    return argparser

def main():
    argparser = configure_argparser()
    args = argparser.parse_args()
    args.invocation_fn(args)

if __name__ == '__main__':
    main()
