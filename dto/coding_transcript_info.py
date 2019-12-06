from collections import namedtuple
import dataclasses

# cds_start and cds_stop are 0-based related to transcript
@dataclasses.dataclass
class CodingTranscriptInfo:
    gene_id: str
    transcript_id: str
    transcript_length: int
    cds_start: int
    cds_stop: int

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

    @classmethod
    def load_transcript_cds_info(cls, cds_annotation_filename):
        transcript_infos = cls.each_in_file(cds_annotation_filename)
        return {tr_info.transcript_id: tr_info  for tr_info in transcript_infos}

    @property
    def cds_length(self):
        return self.cds_stop - self.cds_start
    
    def cds_profile(self, profile):
        return profile[self.cds_start : self.cds_stop]
