import argparse
from ..gzip_utils import open_for_write
from ..utils import flatten
from ..tsv_reader import stream_table_column_highlighted

def configure_argparser(argparser=None):
    if not argparser:
        argparser = argparse.ArgumentParser(prog="join_tables", description = "Join the second table to the first one")
    argparser.add_argument('reference_filename', metavar='reference_table.tsv', help='Table to be expanded in a tab-separated format')
    argparser.add_argument('reference_column', help='Name or indices of column in reference table to join')
    argparser.add_argument('supplementary_filenames', metavar='supplementary_tables.tsv', nargs='*', help='Supplementary tables in a tab-separated format')
    argparser.add_argument('--columns', nargs='*', dest='supplementary_columns', help='Names or indices of columns in supplementary tables to join')

    has_header_group = argparser.add_mutually_exclusive_group(required=True)
    has_header_group.add_argument('--header', action='store_true', dest='has_header', help="Table has header")
    has_header_group.add_argument('--no-header', action='store_false', dest='has_header', help="Tables doesn't have header")
    
    only_matching_group = argparser.add_mutually_exclusive_group(required=True)
    only_matching_group.add_argument('--only-matching', action='store_false', dest='allow_non_matching', help="Don't take lines which are present not in all files")
    only_matching_group.add_argument('--allow-non-matching', action='store_true', dest='allow_non_matching', help="Allow transcripts which are not present in CDS-annotation (they are not clipped)")

    argparser.add_argument('--output-file', '-o', dest='output_file', help="Store results at this path")
    return argparser

def main():
    argparser = configure_argparser(argparser)
    args = argparser.parse_args()
    invoke(args)

def invoke(args):
    if args.supplementary_columns:
        supplementary_columns = args.supplementary_columns
    else:
        supplementary_columns = [args.reference_column] * len(args.supplementary_filenames)

    if len(args.supplementary_filenames) != len(supplementary_columns):
        raise ValueError(f'Number of columns ({args.supplementary_columns}) should match number of tables ({args.supplementary_filenames})')

    supplementary_tables = []
    for (table_filename, column) in zip(args.supplementary_filenames, supplementary_columns):
        table_stream = stream_table_column_highlighted(table_filename, [column], has_header=args.has_header, pop_column=True)
        if args.has_header:
            _, header = next(table_stream)
        mapping = dict(table_stream)
        supplementary_tables.append({'header': header, 'mapping': mapping, 'padding': [''] * len(header)})

    with open_for_write(args.output_file) as output_stream:
        table_stream = stream_table_column_highlighted(args.reference_filename, [args.reference_column], has_header=args.has_header, pop_column=False)
        if args.has_header:
            _, reference_header = next(table_stream)
            header = reference_header + flatten(tbl['header'] for tbl in supplementary_tables)
            print('\t'.join(header), file=output_stream)
        for (key, reference_row) in table_stream:
            if args.allow_non_matching:
                supplementary_rows = [tbl['mapping'].get(key, tbl['padding']) for tbl in supplementary_tables]
            else:
                supplementary_rows = [tbl['mapping'].get(key, None) for tbl in supplementary_tables]
                if not all(supplementary_rows):
                    continue
            row = reference_row + flatten(supplementary_rows)
            print('\t'.join(row), file=output_stream)
