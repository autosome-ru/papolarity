import argparse
from .version import __version__
from .bin import get_coverage, pool_coverage, cds_annotation, clip_cds, \
                 coverage_features, choose_best, \
                 compare_coverage, \
                 adjust_features, plot_distribution, \
                 flatten_coverage, cds_sequence

def configure_argparser(argparser=None):
    if not argparser:
        argparser = argparse.ArgumentParser(prog="papolarity", description = "Main entrypoint of papolarity")
    argparser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    subparsers = argparser.add_subparsers(metavar='subcommand')
    
    subparser_configs = [
        {'cmd': 'get_coverage', 'namespace': get_coverage, 'help': 'Generates coverage from an alignment'},
        {'cmd': 'pool_coverage', 'namespace': pool_coverage, 'help': 'Pool coverage profiles'},
        {'cmd': 'cds_annotation', 'namespace': cds_annotation, 'help': 'Extract CDS annotation in transcriptomic coordinates from genomic annotation'},
        {'cmd': 'clip_cds', 'namespace': clip_cds, 'help': 'Clip any bed file in transcriptomic coordinates to CDS-region'},
        {'cmd': 'coverage_features', 'namespace': coverage_features, 'help': 'Calculate coverage profile features'},
        {'cmd': 'choose_best', 'namespace': choose_best, 'help': 'Choose best element from each group (e.g. best transcipt for each gene)'},
        {'cmd': 'compare_coverage', 'namespace': compare_coverage, 'help': 'Coverage profile comparison'},
        {'cmd': 'plot_distribution', 'namespace': plot_distribution, 'help': 'Plot distributions of features'},
        {'cmd': 'adjust_features', 'namespace': adjust_features, 'help': 'Make length-dependend adjustment of features'},
        {'cmd': 'flatten_coverage', 'namespace': flatten_coverage, 'help': 'Flatten coverage profiles by averaging data through given segments'},
        {'cmd': 'cds_sequence', 'namespace': cds_sequence, 'help': 'Extract CDS sequences from GTF annotation and genome assembly'},
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
