def polarity_score(coverage):
    total_coverage = 0
    numerator = 0
    for i in range(len(coverage)):
        total_coverage = total_coverage + coverage[i]
        numerator += coverage[i] * i
    if total_coverage == 0:
        return 0
    polarity_score = numerator / total_coverage
    normalized_score = (polarity_score * 2 / (len(coverage) - 1)) - 1
    return normalized_score
