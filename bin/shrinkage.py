import sys
from os.path import dirname
sys.path.insert(0, dirname(dirname(__file__)))
import numpy as np
from dto.coverage_comparison_stats import CoverageComparisonStats

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

filename = sys.argv[1] 
slope_data = list(CoverageComparisonStats.each_in_file(filename))
slope_data = [info for info in slope_data  if info.geom_mean_coverage() >= 10]
slope_data = list(CoverageComparisonStats.choose_best_transcript(slope_data))

sorted_slopes = sorted(slope_data, key=lambda info: info.cds_stop - info.cds_start)

header = ['length_mean', 'slope_mean', 'slope_median', 'slope_q90', 'slope_max', 'slope_stddev']
print('\t'.join(header))
for slopes_in_window in sliding_window(sorted_slopes, 500):
    length_mean  = np.mean([info.cds_stop - info.cds_start for info in slopes_in_window])
    slope_mean   = np.mean([info.slope for info in slopes_in_window])
    slope_median   = np.median([info.slope for info in slopes_in_window])
    slope_q90   = np.quantile([info.slope for info in slopes_in_window], 0.9)
    slope_max   = np.max([info.slope for info in slopes_in_window])
    slope_stddev = np.std([info.slope for info in slopes_in_window])
    print('\t'.join([str(round(x, 3)) for x in [length_mean, slope_mean, slope_median, slope_q90, slope_max, slope_stddev]]))
