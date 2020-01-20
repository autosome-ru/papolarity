import argparse
import numpy as np
from ..gzip_utils import open_for_write
from ..tsv_reader import each_in_tsv
from .. import utils

def window_around_idx(arr, idx, window_size, drop_none=False):
    if drop_none:
        left_part = utils.drop_none(arr[:idx])
        element = arr[idx]
        right_part = utils.drop_none(arr[idx + 1:])

        if element is not None:
            arr = left_part + [element] + right_part
        else:
            arr = left_part + right_part
        idx = len(left_part)

    half_window_size = window_size // 2
    if idx < half_window_size:
        return arr[0:window_size]
    elif idx - half_window_size + window_size <= len(arr):
        return arr[(idx - half_window_size):(idx - half_window_size + window_size)]
    else:
        return arr[-window_size:]

def standardize_values(values, standardization, window_size, drop_none=True):
    standardized_values = []
    for (idx, val) in enumerate(values):
        window_values = window_around_idx(values, idx, window_size, drop_none=drop_none)
        mean = np.mean(window_values)
        stddev = np.std(window_values)
        if val is not None:
            standardized_value = standardization(val, val_mean=mean, val_stddev=stddev)
            standardized_values.append(standardized_value)
        else:
            standardized_values.append(None)
    return standardized_values

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

    dtype = lambda x: float(x) if x != '' else None
    data = list(each_in_tsv(args.table))

    fields = list(data[0].keys())
    TRANSFORMED_FIELDS_KEY = object() # unique object not to clash with column name
    for row in data:
        row[TRANSFORMED_FIELDS_KEY] = {}
        for field in [args.sorting_field, *args.fields_to_correct]:
            # we don't want to modify original values because type conversion can screw integer values during output
            row[TRANSFORMED_FIELDS_KEY][field] = dtype(row[field])
    data.sort(key=lambda info: dtype(info[args.sorting_field]))

    modified_data = [row.copy() for row in data]
    output_field_names = fields[:]

    for field in args.fields_to_correct:
        output_field_name = f'{args.prefix}{field}'
        output_field_names.append(output_field_name)
        field_values = [dtype(row[field]) for row in data]
        standardized_field_values = standardize_values(field_values, standardization, args.window_size, drop_none=True)
        for idx in range(len(data)):
            modified_data[idx][output_field_name] = standardized_field_values[idx]

    with open_for_write(args.output_file) as output_stream:
        print('\t'.join(output_field_names), file=output_stream)
        for row in modified_data:
            output_values = [row[field] for field in output_field_names]
            print(utils.tsv_string_empty_none(output_values), file=output_stream)
