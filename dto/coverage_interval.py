import dataclasses
from dto.dataclass_tsv_serializable import DataclassTsvSerializable


@dataclasses.dataclass
class Interval(DataclassTsvSerializable):
    chrom: str
    start:int
    stop: int

@dataclasses.dataclass
class CoverageInterval(Interval):
    coverage: int

@dataclasses.dataclass
class FloatCoverageInterval(Interval):
    coverage: float
