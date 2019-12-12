import dataclasses
from dto.dataclass_tsv_serializable import DataclassTsvSerializable

@dataclasses.dataclass
class Interval(DataclassTsvSerializable):
    chrom: str
    start:int
    stop: int
