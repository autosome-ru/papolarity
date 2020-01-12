# GFF parser adapted from https://techoverflow.net/2013/11/30/a-simple-gff3-parser-in-python/

from collections import namedtuple
import gzip
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

    def attributes_renamed(self, attr_mapping):
        '''
        Use to rename incorrectly named attributes:
          record.attributes_renamed({'wrong_attr_name': 'correct_attr_name', ...})
        Returns a new record with modified attributes; unspecified attributes are left unmodified
        '''
        attributes_modified = {}
        for (orig_attr_name, attr_value) in self.attributes.items():
            if orig_attr_name in attr_mapping:
                attr_name = attr_mapping[orig_attr_name]
            else:
                attr_name = orig_attr_name
            attributes_modified[attr_name] = attr_value
        return GTFRecord(self.contig, self.source, self.type, self.start, self.stop, self.score, self.strand, self.phase, attributes_modified)

    def encoded_attributes(self):
        return self.encode_gtf_attributes(self.attributes)

    @classmethod
    def encode_gtf_attributes(cls, attributes):
        if not attributes:
            return '.'
        attr_strings = [f'{k} {json.dumps(v)};' for k,v in attributes.items()]
        return ' '.join(attr_strings)

    @classmethod
    def each_in_file(cls, filename, multivalue_keys=None):
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
                    "attributes": cls.parse_gtf_attributes(parts[8], multivalue_keys=multivalue_keys),
                }
                yield GTFRecord(**normalized_info)

    @classmethod
    def parse_gtf_attributes(cls, attribute_string, multivalue_keys=None):
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
                    raise Exception(f'Key `{key}` already in attributes.\n'
                                     'Probably you should add this attribute to a set of `multivalue_keys`')
            else:
                if key not in ret:
                    ret[key] = []
                ret[key].append(val)
        return ret

# deprecated
def parse_gtf(filename, multivalue_keys=None):
    yield from GTFRecord.each_in_file(filename, multivalue_keys=multivalue_keys)
