import dataclasses

@dataclasses.dataclass
class DataclassTsvSerializable:
    @classmethod
    def header(cls):
        fields = [field.name for field in dataclasses.fields(cls)]
        return '\t'.join(fields)

    def __str__(self):
        fields = [getattr(self, field.name) for field in dataclasses.fields(self)]
        return '\t'.join(map(str, fields))

    @classmethod
    def from_string(cls, line):
        row = line.rstrip('\n').split('\t')
        return cls(*[field.type(value) for field, value in zip(dataclasses.fields(cls), row)])

    @classmethod
    def each_in_file(cls, filename):
        with open(filename) as f:
            header = f.readline()
            for line in f:
                yield cls.from_string(line)
