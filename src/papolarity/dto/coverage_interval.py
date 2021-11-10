import dataclasses
from typing import Union
from decimal import Decimal
from .dataclass_tsv_serializable import DataclassTsvSerializable

@dataclasses.dataclass(frozen=True)
class CoverageInterval(DataclassTsvSerializable):
    chrom: str
    start:int
    stop: int
    coverage: Union[float, int] = dataclasses.field(metadata={'skip_conversion': True})
    dtype: dataclasses.InitVar[type] = float
    def __post_init__(self, dtype):
        try:
            # object.__setattr_ is used to assign frozen attribute
            object.__setattr__(self, 'coverage', dtype(self.coverage))
        except:
            # if dtype is integer, string "1.2345e6" can't be directly converted to int
            object.__setattr__(self, 'coverage', dtype(Decimal(self.coverage)))
