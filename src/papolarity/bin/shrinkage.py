import argparse
import numpy as np
from ..dto.coverage_comparison_stats import CoverageComparisonStats
from ..gzip_utils import open_for_write

def sliding_window(arr, size):
    try:
        window = []
        iterator = iter(arr)
        for _ in range(size):
            window.append(next(iterator))
        while True:
            yield window[:]
            window.pop(0)
            window.append( next(iterator) )
    except StopIteration:
        pass

def configure_argparser(argparser=None):
    if not argparser:
        argparser = argparse.ArgumentParser(prog="shrinkage", description = "Make length-dependend shrinkage of slope properties")
    argparser.add_argument('coverage_properties', metavar='coverage_properties.tsv', help='Genomic annotation in GTF-format')
    argparser.add_argument('--output-file', '-o', dest='output_file', help="Store results at this path")
    return argparser

def main():
    argparser = configure_argparser(argparser)
    args = argparser.parse_args()
    invoke(args)

def invoke(args):
    slope_data = list(CoverageComparisonStats.each_in_file(args.coverage_properties))
    slope_data = [info for info in slope_data  if info.geom_mean_coverage() >= 10]
    slope_data = list(CoverageComparisonStats.choose_best_transcript(slope_data))

    sorted_slopes = sorted(slope_data, key=lambda info: info.cds_stop - info.cds_start)

    header = ['length_mean', 'slope_mean', 'slope_median', 'slope_q90', 'slope_max', 'slope_stddev']
    with open_for_write(args.output_file) as output_stream:
        print('\t'.join(header), file=output_stream)
        for slopes_in_window in sliding_window(sorted_slopes, 500):
            length_mean  = np.mean([info.cds_stop - info.cds_start for info in slopes_in_window])
            slope_mean   = np.mean([info.slope for info in slopes_in_window])
            slope_median   = np.median([info.slope for info in slopes_in_window])
            slope_q90   = np.quantile([info.slope for info in slopes_in_window], 0.9)
            slope_max   = np.max([info.slope for info in slopes_in_window])
            slope_stddev = np.std([info.slope for info in slopes_in_window])
            output_info = [length_mean, slope_mean, slope_median, slope_q90, slope_max, slope_stddev]
            print('\t'.join([str(round(x, 3)) for x in output_info]), file=output_stream)
