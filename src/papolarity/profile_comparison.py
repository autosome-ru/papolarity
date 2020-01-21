from sklearn.linear_model import LinearRegression
import numpy as np
from math import log2
from collections import namedtuple
from .polarity_score import polarity_score
from .utils import common_subsequence

WeigthedPoints = namedtuple('WeigthedPoints', ['xs', 'ys', 'weights'])

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
        return {'center': None}

    center_info = WeigthedPoints([], [], [])
    weighted_center_info = WeigthedPoints([], [], [])
    every_point_info = WeigthedPoints([], [], [])

    if log_mode:
        rate_transform = log2
    else:
        rate_transform = lambda x: x
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
        center_info.xs.append(rel_coord)
        center_info.ys.append(value)
        center_info.weights.append(1)

    return {
        'center': slope_by_points(center_info.xs, center_info.ys, center_info.weights),
    }

def slope_by_points(xs, ys, weights):
    if len(xs) < 2:
        return None
    xs = np.array(xs)
    model = LinearRegression()
    model.fit(xs.reshape(-1, 1), ys, sample_weight=weights)
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
