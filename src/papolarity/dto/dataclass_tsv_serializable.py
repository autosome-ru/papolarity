import sys
import dataclasses
from ..gzip_utils import open_for_read, open_for_write

@dataclasses.dataclass
class DataclassTsvSerializable:
    # Additional columns (goes after normal ones) which are to be printed/stored but not loaded'''
    computable_properties = [] # [('column_name', 'property_name'), ...]

    @classmethod
    def header(cls):
        fields = [field.name for field in dataclasses.fields(cls)] + [name for (name, prop) in cls.computable_properties]
        return '\t'.join(fields)

    def __str__(self):
        return self.tsv_string()

    def tsv_string(self):
        fields = [getattr(self, field.name) for field in dataclasses.fields(self)] + [getattr(self, prop) for (name, prop) in self.computable_properties]
        return '\t'.join([(str(f) if (f is not None) else '') for f in fields])

    @classmethod
    def from_string(cls, line, **kwargs):
        row = line.rstrip('\n').split('\t')
        attrs = {}
        for field, value in zip(dataclasses.fields(cls), row): # skips computable (and other auxiliary) properties
            if field.metadata.get('skip_conversion', False):
                attrs[field.name] = value
            else:
                converter = field.metadata.get('converter', lambda value: field.type(value))
                attrs[field.name] = converter(value)
        return cls(**attrs, **kwargs)

    @classmethod
    def each_in_file(cls, filename, header=True, force_gzip=None):
        with open_for_read(filename, force_gzip=force_gzip) as f:
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
    def store_tsv(cls, collection, filename, header=True, force_gzip=None):
        with open_for_write(filename, force_gzip=force_gzip) as f:
            cls.print_tsv(collection, f, header=header)
