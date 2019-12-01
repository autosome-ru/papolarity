# GFF parser adapted from https://techoverflow.net/2013/11/30/a-simple-gff3-parser-in-python/

from collections import namedtuple
import gzip
import json
# import urllib.request, urllib.parse, urllib.error

gff_info_fields = ["contig", "source", "type", "start", "stop", "score", "strand", "phase", "attributes"]
class GTFRecord(namedtuple("GTFRecord", gff_info_fields)):
    @property
    def length(self):
        return self.stop - self.start

    def __repr__(self):
        row = [self.contig, self.source, self.type, self.start, self.stop, self.score, self.strand, self.phase, encode_gtf_attributes(self.attributes)]
        row = [elem if elem else '.'  for elem in row]
        return '\t'.join(map(str, row))

    def is_coding(self):
        if self.attributes['gene_type'] != 'protein_coding':
            return False
        if (self.type != 'gene') and (self.attributes['transcript_type'] != 'protein_coding'):
            return False
        return True

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

def encode_gtf_attributes(attributes):
    if not attributes:
        return '.'
    attr_strings = [f'{k} {json.dumps(v)};' for k,v in attributes.items()]
    return ' '.join(attr_strings)

def parse_gtf(filename, attributes_filter=lambda x: x):
    """
    A minimalistic GTF format parser.
    Yields objects that contain info about a single GTF feature.

    Supports transparent gzip decompression.
    """
    # Parse with transparent decompression
    open_func = gzip.open if filename.endswith(".gz") else open
    with open_func(filename, "rt", encoding='utf-8') as infile:
        for line in infile:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            # If this fails, the file format is not standard-compatible
            assert len(parts) == len(gff_info_fields)
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
                "attributes": attributes_filter(parse_gtf_attributes(parts[8]))
            }
            # Alternatively, you can emit the dictionary here, if you need mutability:
            #    yield normalizedInfo
            yield GTFRecord(**normalized_info)

def parse_gtf_attributes(attribute_string):
    """Parse the GTF attribute column and return a dict"""
    if attribute_string == ".":
        return {}
    ret = {}
    # Some records has several attributes with the same key (like `tag`), we treat values for such keys as lists
    multivalue_keys = {'tag', 'ont'}
    for attribute in attribute_string.strip().rstrip(";").split(";"):
        key, value = attribute.strip().split(" ", maxsplit=1)

        if value[0] == '"':
            val = value[1:-1]
        else:
            val = int(value)

        if key not in multivalue_keys:
            if key not in ret:
                ret[key] = val
            else:
                raise Exception(f'Key `{key}` already in attributes')
        else:
            if key not in ret:
                ret[key] = []
            ret[key].append(val)
    return ret
