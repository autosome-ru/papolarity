import dataclasses
from .dataclass_tsv_serializable import DataclassTsvSerializable

# cds_start and cds_stop are 0-based related to transcript
@dataclasses.dataclass(frozen=True)
class CodingTranscriptInfo(DataclassTsvSerializable):
    gene_id: str
    transcript_id: str
    transcript_length: int
    cds_start: int
    cds_stop: int

    @classmethod
    def load_transcript_cds_info(cls, cds_annotation_filename):
        transcript_infos = cls.each_in_file(cds_annotation_filename)
        return {tr_info.transcript_id: tr_info  for tr_info in transcript_infos}

    @property
    def cds_length(self):
        return self.cds_stop - self.cds_start
    
    def cds_profile(self, profile):
        return profile[self.cds_start : self.cds_stop]
