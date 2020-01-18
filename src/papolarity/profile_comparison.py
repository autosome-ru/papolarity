from sklearn.linear_model import LinearRegression
import numpy as np
from math import log

from collections import namedtuple
WeigthedPoints = namedtuple('WeigthedPoints', ['xs', 'ys', 'weights'])

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
    total_coverage_control = np.sum(control_sums)
    total_coverage_experiment = np.sum(experiment_sums)

    assert len(segmentation.segments) == len(control_sums) == len(experiment_sums)
    profile_len = segmentation.segmentation_length

    center_info = WeigthedPoints([], [], [])
    weighted_center_info = WeigthedPoints([], [], [])
    every_point_info = WeigthedPoints([], [], [])

    num_segments = segmentation.num_segments
    if log_mode:
        rate_transform = log
    else:
        rate_transform = lambda x: x
    for (segment, control_segment_sum, experiment_segment_sum) in zip(segmentation.segments, control_sums, experiment_sums):
        mean_control = control_segment_sum / segment.length
        mean_experiment = experiment_segment_sum / segment.length
        normalized_mean_control = (mean_control + 1) / (total_coverage_control + num_segments)
        normalized_mean_experiment = (mean_experiment + 1) / (total_coverage_experiment + num_segments)
        detrended_profile = normalized_mean_experiment / normalized_mean_control
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
    xs = np.array(xs)
    model = LinearRegression()
    model.fit(xs.reshape(-1, 1), ys, sample_weight=weights)
    return model.coef_[0]

def discrepancy(control_profile, experiment_profile, segmentation):
    control_sums = segmentwise_sums(segmentation, control_profile)
    experiment_sums = segmentwise_sums(segmentation, experiment_profile)
    assert len(control_sums) == len(experiment_sums) == len(segmentation.segments)
    return discrepancy_by_segment_counts(control_sums, experiment_sums)

def discrepancy_by_segment_counts(control_sums, experiment_sums):
    assert len(control_sums) == len(experiment_sums)

    control_profile_sum = np.sum(control_sums)
    experiment_profile_sum = np.sum(experiment_sums)

    control_normed = control_sums / control_profile_sum  if control_profile_sum != 0  else control_sums
    experiment_normed = experiment_sums / experiment_profile_sum  if experiment_profile_sum != 0  else experiment_sums
    difference = 0
    for (control_segment_sum, experiment_segment_sum) in zip(control_normed, experiment_normed):
        difference += abs(control_segment_sum - experiment_segment_sum)
    return difference
