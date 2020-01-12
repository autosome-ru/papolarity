from gzip_utils import open_for_read

def stream_table_column_highlighted(filename, column, has_header, pop_column=True):
    with open_for_read(filename) as file:
        if has_header:
            header = next(file)
            header = header.rstrip("\n").split("\t")
            try:
                column_idx = int(column)
            except ValueError:
                column_idx = header.index(column)
            column_name = header[column_idx]
            if pop_column:
                header.pop(column_idx)
            yield (column_name, header)
        else:
            column_idx = int(column)

        for line in file:
            row = line.rstrip("\n").split("\t")
            key = row[column_idx]
            if pop_column:
                row.pop(column_idx)
            yield (key, row)
