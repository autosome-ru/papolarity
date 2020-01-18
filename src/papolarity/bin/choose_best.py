import argparse
import itertools
from ..gzip_utils import open_for_write
from ..tsv_reader import stream_table_column_highlighted

def configure_argparser(argparser=None):
    if not argparser:
        argparser = argparse.ArgumentParser(
            prog = "choose_best",
            description = "Choose best element from each group (e.g. best transcipt for each gene)",
        )
    argparser.add_argument('table', metavar='table.tsv', help='Input table in tab-separated format')
    argparser.add_argument('column', help="Take the biggest or the lowest transcript according to this column")
    argparser.add_argument('criteria', choices=['max', 'min'], help="Criteria to choose transcript: `max` takes the biggest, `min` - the lowest")
    argparser.add_argument('--group-by', default='gene_id', dest='group_by_column', help="Column to group transcripts")
    
    header_group = argparser.add_mutually_exclusive_group(required=True)
    header_group.add_argument('--header', action='store_true', dest='has_header', help="Table has header")
    header_group.add_argument('--no-header', action='store_false', dest='has_header', help="Tables doesn't have header")

    argparser.add_argument('--output-file', '-o', dest='output_file', help="Store results at this path")
    return argparser

def main():
    argparser = configure_argparser()
    args = argparser.parse_args()
    invoke(args)

def invoke(args):
    dtype = float

    get_group = lambda row: row[0][0]
    get_value = lambda row: dtype(row[0][1])
    get_row = lambda row: row[1]

    with open_for_write(args.output_file) as output_stream:
        table_stream = stream_table_column_highlighted(args.table, [args.group_by_column, args.column], has_header=args.has_header, pop_column=False)
        if args.has_header:
            _, header = next(table_stream)
            print('\t'.join(header), file=output_stream)
        
        data = list(table_stream)
        data.sort(key=get_group)
        for group, rows in itertools.groupby(data, key=get_group):
            if args.criteria == 'max':
                row = get_row(max(rows, key=get_value))
            elif args.criteria == 'min':
                row = get_row(min(rows, key=get_value))
            else:
                raise ValueError(f'Unknown criteria `{criteria}` (only min/max allowed)')
            print('\t'.join(row), file=output_stream)
