# GFF parser adapted from https://techoverflow.net/2013/11/30/a-simple-gff3-parser-in-python/

from collections import namedtuple
import gzip
import json
# import urllib.request, urllib.parse, urllib.error

gff_info_fields = ["contig", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
class GTFRecord(namedtuple("GFFRecord", gff_info_fields)):
    def length(self):
        return self.end - self.start + 1

    def __repr__(self):
        row = [self.contig, self.source, self.type, self.start, self.end, self.score, self.strand, self.phase, encode_gtf_attributes(self.attributes)]
        row = [elem if elem else '.'  for elem in row]
        return '\t'.join(map(str, row))

def encode_gtf_attributes(attributes):
    if not attributes:
        return '.'
    attr_strings = [f'{k} {json.dumps(v)};' for k,v in attributes.items()]
    return ' '.join(attr_strings)

def parse_gtf(filename):
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
            # Normalize data
            normalized_info = {
                "contig": None if parts[0] == "." else parts[0],
                "source": None if parts[1] == "." else parts[1],
                "type": None if parts[2] == "." else parts[2],
                "start": None if parts[3] == "." else int(parts[3]), # 1-based
                "end": None if parts[4] == "." else int(parts[4]), # 1-based
                "score": None if parts[5] == "." else float(parts[5]),
                "strand": None if parts[6] == "." else parts[6],
                "phase": None if parts[7] == "." else parts[7],
                "attributes": parse_gtf_attributes(parts[8])
            }
            # Alternatively, you can emit the dictionary here, if you need mutability:
            #    yield normalizedInfo
            yield GTFRecord(**normalized_info)

def parse_gtf_attributes(attribute_string):
    """Parse the GTF attribute column and return a dict"""
    if attribute_string == ".":
        return {}
    ret = {}
    for attribute in attribute_string.strip().rstrip(";").split(";"):
        key, value = attribute.strip().split(" ")
        if value[0] == '"':
            ret[key] = value[1:-1]
        else:
            ret[key] = int(value)
    return ret
