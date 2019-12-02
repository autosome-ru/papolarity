from collections import namedtuple

# cds_start and cds_stop are 0-based related to transcript
class CodingTranscriptInfo(namedtuple("CodingTranscriptInfo", ['gene_id', 'transcript_id', 'transcript_length', 'cds_start', 'cds_stop'])):
    @classmethod
    def header(cls):
        return '\t'.join(['gene_id', 'transcript_id', 'transcript_length', 'cds_start', 'cds_stop'])

    def __repr__(self):
        return '\t'.join(map(str, [self.gene_id, self.transcript_id, self.transcript_length, self.cds_start, self.cds_stop]))

    @classmethod
    def from_string(cls, line):
        row = line.rstrip("\n").split("\t")
        gene_id, transcript_id, transcript_length, cds_start, cds_stop = row
        return CodingTranscriptInfo(gene_id, transcript_id, int(transcript_length), int(cds_start), int(cds_stop))

    @classmethod
    def each_from_file(cls, filename):
        with open(filename) as f:
            header = f.readline()
            for line in f:
                yield CodingTranscriptInfo.from_string(line)

    @classmethod
    def load_transcript_cds_info(cls, cds_annotation_filename):
        transcript_infos = cls.each_from_file(cds_annotation_filename)
        return {tr_info.transcript_id: tr_info  for tr_info in transcript_infos}

    @property
    def cds_length(self):
        return self.cds_stop - self.cds_start
    
    def cds_profile(self, profile):
        return profile[self.cds_start : self.cds_stop]
