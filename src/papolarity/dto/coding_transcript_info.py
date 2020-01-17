import dataclasses
from .dataclass_tsv_serializable import DataclassTsvSerializable

# cds_start and cds_stop are 0-based related to transcript
@dataclasses.dataclass(frozen=True)
class CodingTranscriptInfo(DataclassTsvSerializable):
    gene_id: str
    transcript_id: str
    transcript_length: int
    cds_start: int = dataclasses.field(metadata={'converter': lambda x: int(x) if x != '' else None})
    cds_stop: int = dataclasses.field(metadata={'converter': lambda x: int(x) if x != '' else None})

    computable_properties = [('cds_length', 'cds_length'), ]

    @classmethod
    def load_transcript_cds_info(cls, cds_annotation_filename):
        transcript_infos = cls.each_in_file(cds_annotation_filename)
        return {tr_info.transcript_id: tr_info  for tr_info in transcript_infos}

    @property
    def cds_length(self):
        if self.is_coding:
            return self.cds_stop - self.cds_start
        else:
            return None
    
    def cds_profile(self, profile):
        return profile[self.cds_start : self.cds_stop]

    @property
    def is_coding(self):
        return self.cds_start and self.cds_stop
