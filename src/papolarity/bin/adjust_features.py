import argparse
import numpy as np
from ..gzip_utils import open_for_write
from ..tsv_reader import each_in_tsv

def sliding_row_with_window(arr, size):
    for idx in range(0, size // 2):
        yield (arr[idx], arr[0:size])
    for idx in range(size // 2, len(arr) - size // 2):
        yield (arr[idx], arr[(idx - size // 2):(idx - size // 2 + size)])
    for idx in range(len(arr) - size // 2, len(arr)):
        yield (arr[idx], arr[-size:])

# stddev --> 1
def standardize_stddev(val, val_mean, val_stddev):
    return val_mean + (val - val_mean) / val_stddev

# mean --> 0
def standardize_mean(val, val_mean, val_stddev):
    return val - val_mean

# mean --> 0
# stddev --> 1
# (aka z-score)
def standardize_zscore(val, val_mean, val_stddev):
    return (val - val_mean) / val_stddev

def configure_argparser(argparser=None):
    if not argparser:
        argparser = argparse.ArgumentParser(prog="adjust_features", description = "Make length-dependend adjustment of features")
    argparser.add_argument('table', metavar='table.tsv', help='Table in tab-separated format')
    argparser.add_argument('--sort-field', dest='sorting_field', required=True, help='Field to sort a table')
    argparser.add_argument('--fields', nargs='*', dest='fields_to_correct', default=[], help='List of fields to correct')
    argparser.add_argument('--prefix', default='adjusted_', help='Prefix for corrected column name')
    argparser.add_argument('--window', metavar='SIZE', dest='window_size', type=int, default=100, help='Size of sliding window (default: %(default)s)')
    argparser.add_argument('--mode', choices=['zero_mean', 'unit_stddev', 'z-score'], default='z-score', help='How to standardize statistics (default: %(default)s)')
    argparser.add_argument('--output-file', '-o', dest='output_file', help="Store results at this path")
    return argparser

def main():
    argparser = configure_argparser(argparser)
    args = argparser.parse_args()
    invoke(args)

def invoke(args):
    if args.mode == 'zero_mean':
        standardization = standardize_mean
    elif args.mode == 'unit_stddev':
        standardization = standardize_stddev
    elif args.mode == 'z-score':
        standardization = standardize_zscore
    else:
        raise ValueError(f'Unknown mode `{args.mode}`')

    data = list(each_in_tsv(args.table))

    fields = list(data[0].keys())
    TRANSFORMED_FIELDS_KEY = object() # unique object not to clash with column name
    for row in data:
        row[TRANSFORMED_FIELDS_KEY] = {}
        for field in [args.sorting_field, *args.fields_to_correct]:
            # we don't want to modify original values because type conversion can screw integer values during output
            row[TRANSFORMED_FIELDS_KEY][field] = float(row[field])
    data.sort(key=lambda info: info[TRANSFORMED_FIELDS_KEY][args.sorting_field])

    output_fields = fields[:]
    for field in args.fields_to_correct:
        output_fields.append(f'{args.prefix}{field}')

    with open_for_write(args.output_file) as output_stream:
        print('\t'.join(output_fields), file=output_stream)
        for row, window in sliding_row_with_window(data, args.window_size):
            modified_row = dict(row)
            for field in args.fields_to_correct:
                field_vals = [window_row[TRANSFORMED_FIELDS_KEY][field] for window_row in window]
                modified_row[f'{args.prefix}{field}'] = standardization(row[TRANSFORMED_FIELDS_KEY][field], val_mean=np.mean(field_vals), val_stddev=np.std(field_vals))
            print('\t'.join([str(modified_row[field]) for field in output_fields]), file=output_stream)
