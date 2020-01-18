import numpy as np

def polarity_score(coverage):
    coverage = np.array(coverage)
    normalized_coverage = coverage / coverage.sum()
    positions = np.linspace(-1, 1, len(coverage))
    return np.dot(normalized_coverage, positions)
