import dataclasses
from dto.dataclass_tsv_serializable import DataclassTsvSerializable

@dataclasses.dataclass(order=True)
class Interval(DataclassTsvSerializable):
    '''bed-coordinates [a, b)'''
    chrom: str
    start: int
    stop: int

    def __post_init__(self):
        if not isinstance(self.chrom, str):
            raise TypeError(f'Chromosome should be a string but was `{self.chrom}`')
        if not isinstance(self.start, int):
            raise TypeError(f'Interval start should be integer but was `{self.start}`')
        if not isinstance(self.stop, int):
            raise TypeError(f'Interval stop should be integer but was `{self.stop}`')
        if self.start >= self.stop:
            raise ValueError(f'Interval start={self.start} should be less than stop={self.stop}')
