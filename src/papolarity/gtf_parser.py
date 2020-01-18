# GFF parser adapted from https://techoverflow.net/2013/11/30/a-simple-gff3-parser-in-python/

from collections import namedtuple
import json
from .gzip_utils import open_for_read

_gff_info_fields = ["contig", "source", "type", "start", "stop", "score", "strand", "phase", "attributes"]
class GTFRecord(namedtuple("GTFRecord", _gff_info_fields)):
    @property
    def length(self):
        return self.stop - self.start

    def __repr__(self):
        row = [self.contig, self.source, self.type, self.start, self.stop, self.score, self.strand, self.phase, self.encoded_attributes()]
        row = [elem if elem else '.'  for elem in row]
        return '\t'.join(map(str, row))

    def contain_position(self, pos):
        return self.start <= pos < self.stop

    def in_upstream_of(self, pos):
        if self.strand == '+':
            return self.stop <= pos
        elif self.strand == '-':
            return pos < self.start

    def in_downstream_of(self, pos):
        if self.strand == '+':
            return pos < self.start
        elif self.strand == '-':
            return self.stop <= pos

    def attributes_filtered(self, relevant_attributes):
        '''
        Use to drop unnecessary attributes and make a program more memory-efficient:
          record.attributes_filtered({'attr_1', 'attr_2', ...})
        Returns a new record with modified attributes
        '''
        attributes_filtered = {k: v  for (k,v) in self.attributes.items()  if k in relevant_attributes}
        return GTFRecord(self.contig, self.source, self.type, self.start, self.stop, self.score, self.strand, self.phase, attributes_filtered)

    def encoded_attributes(self):
        return self.encode_gtf_attributes(self.attributes)

    @classmethod
    def encode_gtf_attributes(cls, attributes):
        if not attributes:
            return '.'
        attr_strings = [f'{k} {json.dumps(v)};' for k,v in attributes.items()]
        return ' '.join(attr_strings)

    @classmethod
    def each_in_file(cls, filename, multivalue_keys=None, ignore_unknown_multivalues=False):
        """
        A minimalistic GTF format parser.
        Yields objects that contain info about a single GTF feature.

        Supports transparent gzip decompression.
        """
        with open_for_read(filename, encoding='utf-8') as infile:
            for line in infile:
                if line.startswith("#"): continue
                parts = line.strip().split("\t")
                # If this fails, the file format is not standard-compatible
                assert len(parts) == len(_gff_info_fields)
                assert parts[0] != '.'  # contig
                assert parts[2] != '.'  # type
                assert parts[3] != '.'  # start
                assert parts[4] != '.'  # stop
                assert parts[6] in {'+', '-'}

                # Normalize data
                normalized_info = {
                    "contig": parts[0],
                    "source": None if parts[1] == "." else parts[1],
                    "type":   parts[2],
                    "start":  int(parts[3]) - 1, # 0-based, included
                    "stop":   int(parts[4]),     # 0-based, excluded
                    "score": None if parts[5] == "." else float(parts[5]),
                    "strand": parts[6],
                    "phase": None if parts[7] == "." else parts[7],
                    "attributes": cls.parse_gtf_attributes(parts[8], multivalue_keys=multivalue_keys, ignore_unknown_multivalues=ignore_unknown_multivalues),
                }
                yield GTFRecord(**normalized_info)

    @classmethod
    def parse_gtf_attributes(cls, attribute_string, multivalue_keys=None, ignore_unknown_multivalues=False):
        """Parse the GTF attribute column and return a dict"""
        if attribute_string == ".":
            return {}
        ret = {}
        for attribute in attribute_string.strip().rstrip(";").split(";"):
            key, value = attribute.strip().split(" ", maxsplit=1)

            if value[0] == '"':
                val = value[1:-1]
            else:
                val = int(value)

            # Some records has several attributes with the same key
            # (e.g. `tag tag_1; tag tag_2;`)
            # We treat values for such keys as lists
            if key not in multivalue_keys:
                if key not in ret:
                    ret[key] = val
                else:
                    if not ignore_unknown_multivalues:
                        raise Exception(f'Key `{key}` already in attributes.\n'
                                         'Probably you should add this attribute to a set of `multivalue_keys`')
            else:
                if key not in ret:
                    ret[key] = []
                ret[key].append(val)
        return ret
