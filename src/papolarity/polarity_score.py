import numpy as np

def polarity_score(coverage):
    coverage = np.array(coverage)
    total_coverage = coverage.sum()
    if total_coverage == 0:
        return None
    normalized_coverage = coverage / total_coverage
    positions = np.linspace(-1, 1, len(coverage))
    return np.dot(normalized_coverage, positions)
