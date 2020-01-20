import dataclasses
from typing import List, Any
from .dataclass_tsv_serializable import DataclassTsvSerializable
from ..utils import tsv_string_empty_none

@dataclasses.dataclass(order=True, frozen=True)
class Interval(DataclassTsvSerializable):
    '''bed-coordinates [a, b)'''
    chrom: str
    start: int
    stop: int
    rest: List[Any] = dataclasses.field(default_factory=list, metadata={'skip_conversion': True})

    def __post_init__(self):
        if not isinstance(self.chrom, str):
            raise TypeError(f'Chromosome should be a string but was `{self.chrom}`')
        if not isinstance(self.start, int):
            raise TypeError(f'Interval start should be integer but was `{self.start}`')
        if not isinstance(self.stop, int):
            raise TypeError(f'Interval stop should be integer but was `{self.stop}`')
        if self.start >= self.stop:
            raise ValueError(f'Interval start={self.start} should be less than stop={self.stop}')

    @classmethod
    def from_string(cls, line, **kwargs):
        row = line.rstrip('\n').split('\t')
        chrom, start, stop, *rest = row
        attrs = {'chrom': chrom, 'start': int(start), 'stop': int(stop), 'rest': rest}
        return cls(**attrs, **kwargs)

    def tsv_string(self):
        fields = [self.chrom, self.start, self.stop, *self.rest]
        return tsv_string_empty_none(fields)

    @property
    def length(self):
        return self.stop - self.start
