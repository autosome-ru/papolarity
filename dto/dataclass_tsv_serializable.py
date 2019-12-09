import sys
import dataclasses
import gzip

@dataclasses.dataclass
class DataclassTsvSerializable:
    @classmethod
    def header(cls):
        fields = [field.name for field in dataclasses.fields(cls)]
        return '\t'.join(fields)

    def __str__(self):
        return self.tsv_string()

    def tsv_string(self):
        fields = [getattr(self, field.name) for field in dataclasses.fields(self)]
        return '\t'.join(map(str, fields))

    @classmethod
    def from_string(cls, line, **kwargs):
        row = line.rstrip('\n').split('\t')
        attrs = {}
        for field, value in zip(dataclasses.fields(cls), row):
            if field.metadata.get('skip_conversion', False):
                attrs[field.name] = value
            else:
                attrs[field.name] = field.type(value)
        return cls(**attrs, **kwargs)

    @classmethod
    def _open_func(gzip, filename):
        '''
        If `gzip` is True or False, use corresponding open function.
        If `gzip` is None, it's guessed based on filename extension
        '''
        if (gzip is True) or ((gzip is None) and filename.endswith('.gz')):
            open_func = gzip.open
        else (gzip is False) or ((gzip is None) and not filename.endswith('.gz')):
            open_func = open
        else:
            raise ValueError("gzip should be one of True/False/None")

    @classmethod
    def each_in_file(cls, filename, header=True, gzip=None):
        open_func = cls._open_func(gzip=gzip, filename=file)
        with open_func(filename, 'rt') as f:
            if header:
                f.readline() # skip header
            for line in f:
                yield cls.from_string(line)

    @classmethod
    def print_tsv(cls, collection, file=sys.stdout, header=True):
        if header:
            print(cls.header(), file=file)
        for record in collection:
            print(record.tsv_string(), file=file)

    @classmethod
    def store_tsv(cls, collection, filename, header=True, gzip=None):
        open_func = cls._open_func(gzip=gzip, filename=file)
        with open_func(filename, 'wt') as f:
            cls.print_tsv(collection, f, header=header)
    