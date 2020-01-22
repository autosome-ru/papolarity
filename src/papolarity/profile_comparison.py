from sklearn.linear_model import LinearRegression
import numpy as np
from math import log2
from .polarity_score import polarity_score
from .utils import common_subsequence, drop_none

def random_probabilities(num_bins):
    values = np.random.rand(num_bins)
    return values / values.sum()

def comparison_infos(transcript_id, control_profile, experiment_profile, segmentation):
    assert len(control_profile) == len(experiment_profile) == segmentation.segmentation_length
    control_sums = segmentwise_sums(segmentation, control_profile)
    experiment_sums = segmentwise_sums(segmentation, experiment_profile)
    control_total_coverage = np.sum(control_sums)
    experiment_total_coverage = np.sum(experiment_sums)
    bootstrap_fraction = 0.1
    num_iterations = 30
    slope_samples = []
    slopelog_samples = []
    for _ in range(num_iterations):
        control_sums_corrected = control_sums + bootstrap_fraction * control_total_coverage * random_probabilities(len(control_sums))
        experiment_sums_corrected = experiment_sums + bootstrap_fraction * experiment_total_coverage * random_probabilities(len(experiment_sums))
        slope_sample = slope_by_segment_counts(control_sums_corrected, experiment_sums_corrected, segmentation, log_mode=False)
        slopelog_sample = slope_by_segment_counts(control_sums_corrected, experiment_sums_corrected, segmentation, log_mode=True)
        slope_samples.append(slope_sample)
        slopelog_samples.append(slopelog_sample)

    slope_samples = drop_none(slope_samples)
    if len(slope_samples) > 0:
        slope_stats = {'slope_min': np.min(slope_samples), 'slope_max': np.max(slope_samples), 'slope_median': np.median(slope_samples), 'slope_mean': np.mean(slope_samples), 'slope_stddev': np.std(slope_samples), 'slope_rel_stddev': np.abs(np.std(slope_samples) / np.mean(slope_samples))}
    else:
        slope_stats = {'slope_min': None, 'slope_max': None, 'slope_mean': None, 'slope_median': None, 'slope_stddev': None, 'slope_rel_stddev': None,}

    slopelog_samples = drop_none(slopelog_samples)
    if len(slopelog_samples) > 0:
        slopelog_stats = {'slopelog_min': np.min(slopelog_samples), 'slopelog_max': np.max(slopelog_samples), 'slopelog_median': np.median(slopelog_samples), 'slopelog_mean': np.mean(slopelog_samples), 'slopelog_stddev': np.std(slopelog_samples), 'slopelog_rel_stddev': np.abs(np.std(slopelog_samples) / np.mean(slopelog_samples)),}
    else:
        slopelog_stats = {'slopelog_min': None, 'slopelog_max': None, 'slopelog_median': None, 'slopelog_mean': None, 'slopelog_stddev': None, 'slopelog_rel_stddev': None,}


    info = {
        'transcript_id': transcript_id,
        'control_median': np.median(control_profile), 'experiment_median': np.median(experiment_profile),
        'control_median_segments': np.median(control_sums), 'experiment_median_segments': np.median(experiment_sums),
        'slope': slope_by_segment_counts(control_sums, experiment_sums, segmentation, log_mode=False),
        'slopelog': slope_by_segment_counts(control_sums, experiment_sums, segmentation, log_mode=True),
        'num_segments': segmentation.num_segments,'control_total_coverage': control_total_coverage,'experiment_total_coverage': experiment_total_coverage,
        **slope_stats, **slopelog_stats,
        'l1_distance': l1_distance_by_segment_counts(control_sums, experiment_sums),
        'polarity_diff': polarity_diff(control_profile, experiment_profile),
    }
    return info

def polarity_diff(control_profile, experiment_profile):
    control_polarity = polarity_score(control_profile)
    exp_polarity = polarity_score(experiment_profile)
    if (exp_polarity is None) or (control_polarity is None):
        return None
    return exp_polarity - control_polarity

def segmentation_stops(segmentation):
    stops = [s.start for s in segmentation.segments]
    stops.append(segmentation.segments[-1].stop)
    return stops

def segmentwise_sums(segmentation, profile):
    profile_cumsum = np.hstack([0, np.cumsum(profile)])
    return np.diff(profile_cumsum[segmentation_stops(segmentation)])

    control_sums = segmentwise_sums(segmentation, control_profile)
    experiment_sums = segmentwise_sums(segmentation, experiment_profile)


def slope_by_profiles(control_profile, experiment_profile, segmentation, log_mode=False):
    control_sums = segmentwise_sums(segmentation, control_profile)
    experiment_sums = segmentwise_sums(segmentation, experiment_profile)
    return slope_by_segment_counts(control_sums, experiment_sums, segmentation, log_mode=log_mode)

def slope_by_segment_counts(control_sums, experiment_sums, segmentation, log_mode=False):
    assert len(segmentation.segments) == len(control_sums) == len(experiment_sums)
    profile_len = segmentation.segmentation_length
    num_segments = segmentation.num_segments

    total_coverage_control = np.sum(control_sums)
    total_coverage_experiment = np.sum(experiment_sums)
    if (total_coverage_control == 0) or (total_coverage_experiment ==0):
        return None

    if log_mode:
        rate_transform = log2
    else:
        rate_transform = lambda x: x

    xs, ys = [], []
    for (segment, control_segment_sum, experiment_segment_sum) in zip(segmentation.segments, control_sums, experiment_sums):
        if (control_segment_sum == 0) and (experiment_segment_sum == 0):
            continue

        # We add pseudocount of 1 to each total coverage of each segment
        mean_control = (control_segment_sum + 1) / (total_coverage_control + num_segments)
        mean_experiment = (experiment_segment_sum + 1) / (total_coverage_experiment + num_segments)
        detrended_profile = mean_experiment / mean_control

        value = rate_transform(detrended_profile)

        coord = (segment.start + segment.stop - 1) / 2
        rel_coord = coord / profile_len

        # mode: 'center'
        xs.append(rel_coord)
        ys.append(value)

    return slope_by_points(xs, ys)

def slope_by_points(xs, ys):
    if len(xs) < 2:
        return None
    xs = np.array(xs)
    model = LinearRegression()
    model.fit(xs.reshape(-1, 1), ys)
    return model.coef_[0]

def l1_distance(control_profile, experiment_profile, segmentation):
    control_sums = segmentwise_sums(segmentation, control_profile)
    experiment_sums = segmentwise_sums(segmentation, experiment_profile)
    assert len(control_sums) == len(experiment_sums) == len(segmentation.segments)
    return l1_distance_by_segment_counts(control_sums, experiment_sums)

def l1_distance_by_segment_counts(control_sums, experiment_sums):
    assert len(control_sums) == len(experiment_sums)
    control_profile_sum = np.sum(control_sums)
    experiment_profile_sum = np.sum(experiment_sums)
    if (control_profile_sum == 0) or (experiment_profile_sum == 0):
        return None
    control_normed = control_sums / control_profile_sum
    experiment_normed = experiment_sums / experiment_profile_sum
    return np.sum(np.abs(control_normed - experiment_normed))

def _coverage_contig_fn(transcript_coverage):
    return transcript_coverage.transcript_id

def _segmentation_contig_fn(segment):
    return segment.chrom

def align_profile_streams(profile_streams):
    '''
    yields tuples: (transcript_id, profiles)
    '''
    keys = [_coverage_contig_fn] * len(profile_streams)
    yield from common_subsequence(profile_streams, key=keys, check_sorted=True)

def align_profile_streams_to_segmentation(segmentation_stream, profile_streams):
    '''
    yields tuples: (transcript_id, (segmentation, *profiles))
    '''
    streams = [segmentation_stream, *profile_streams]
    keys = [_segmentation_contig_fn] + [_coverage_contig_fn] * len(profile_streams)
    yield from common_subsequence(streams, key=keys, check_sorted=True)

def compare_coverage_streams(segmentation_stream, control_coverage_profiles, experiment_coverage_profiles):
    aligned_stream = align_profile_streams_to_segmentation(
        segmentation_stream,
        [control_coverage_profiles, experiment_coverage_profiles]
    )
    for (transcript_id, (segmentation, control_coverage, experiment_coverage)) in aligned_stream:
        yield comparison_infos(transcript_id,
                               control_coverage.coverage,
                               experiment_coverage.coverage,
                               segmentation)
