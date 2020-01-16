import csv
from .gzip_utils import open_for_read

def stream_table_column_highlighted(filename, columns, has_header, pop_column=True):
    with open_for_read(filename) as file:
        if has_header:
            header = next(file)
            header = header.rstrip("\n").split("\t")
            try:
                column_indices = [int(column) for column in columns]
            except ValueError:
                column_indices = [header.index(column) for column in columns]
            column_names = tuple([header[column_idx] for column_idx in column_indices])
            if pop_column:
                for column_idx in column_indices:
                    header.pop(column_idx)
            yield (column_names, header)
        else:
            column_idx = int(column)

        for line in file:
            row = line.rstrip("\n").split("\t")
            key = tuple([row[column_idx] for column_idx in column_indices])
            if pop_column:
                for column_idx in column_indices:
                    row.pop(column_idx)
            yield (key, row)

def each_in_tsv(filename):
    with open_for_read(filename) as input_stream:
        reader = csv.DictReader(input_stream, delimiter='\t')
        for row in reader:
            yield row
