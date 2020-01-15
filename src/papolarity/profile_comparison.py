from sklearn.linear_model import LinearRegression
import numpy as np
from math import log

from collections import namedtuple
WeigthedPoints = namedtuple('WeigthedPoints', ['xs', 'ys', 'weights'])

def slope_by_profiles(control_profile, experiment_profile, segmentation):
    cumsum_coverage_control = np.hstack([0, np.cumsum(control_profile)])
    cumsum_coverage_experiment = np.hstack([0, np.cumsum(experiment_profile)])
    total_coverage_control = cumsum_coverage_control[-1]
    total_coverage_experiment = cumsum_coverage_experiment[-1]
    assert len(control_profile) == len(experiment_profile)
    assert len(control_profile) == segmentation.segmentation_length
    profile_len = len(control_profile)

    center_info = WeigthedPoints([], [], [])
    weighted_center_info = WeigthedPoints([], [], [])
    every_point_info = WeigthedPoints([], [], [])

    num_segments = segmentation.num_segments
    for segment in segmentation.segments:
        start = segment.start
        stop = segment.stop
        mean_control = (cumsum_coverage_control[stop] - cumsum_coverage_control[start]) / (stop - start)
        mean_experiment = (cumsum_coverage_experiment[stop] - cumsum_coverage_experiment[start]) / (stop - start)
        normalized_mean_control = (mean_control + 1) / (total_coverage_control + num_segments)
        normalized_mean_experiment = (mean_experiment + 1) / (total_coverage_experiment + num_segments)
        detrended_profile = normalized_mean_experiment / normalized_mean_control
        value = log(detrended_profile)

        coord = (start + stop - 1) / 2
        rel_coord = coord / profile_len

        # mode: 'center'
        center_info.xs.append(rel_coord)
        center_info.ys.append(value)
        center_info.weights.append(1)

        # mode: 'weighted_center'
        weight = stop - start
        weighted_center_info.xs.append(rel_coord)
        weighted_center_info.ys.append(value)
        weighted_center_info.weights.append(weight)

        # mode: 'every_point'
        for pos in range(start, stop):
            rel_pos = pos / profile_len
            every_point_info.xs.append(rel_pos)
            every_point_info.ys.append(value)
            every_point_info.weights.append(1)

    return {
        'center': slope_by_points(center_info.xs, center_info.ys, center_info.weights),
        'weighted_center': slope_by_points(weighted_center_info.xs, weighted_center_info.ys, weighted_center_info.weights),
        'every_point': slope_by_points(every_point_info.xs, every_point_info.ys, every_point_info.weights),
    }

def slope_by_points(xs, ys, weights):
    xs = np.array(xs)
    model = LinearRegression()
    model.fit(xs.reshape(-1, 1), ys, sample_weight=weights)
    return model.coef_[0]

def profile_difference(control_profile, experiment_profile, segmentation):
    assert len(control_profile) == len(experiment_profile)

    cumsum_coverage_control = np.hstack([0, np.cumsum(control_profile)])
    cumsum_coverage_experiment = np.hstack([0, np.cumsum(experiment_profile)])
    control_profile_sum = cumsum_coverage_control[-1]
    experiment_profile_sum = cumsum_coverage_experiment[-1]

    cumsum_coverage_control_normed = cumsum_coverage_control / control_profile_sum  if control_profile_sum != 0  else cumsum_coverage_control
    cumsum_coverage_experiment_normed = cumsum_coverage_experiment / experiment_profile_sum  if experiment_profile_sum != 0  else cumsum_coverage_experiment
    difference = 0
    for segment in segmentation.segments:
        start = segment.start
        stop = segment.stop
        control_segment = cumsum_coverage_control_normed[stop] - cumsum_coverage_control_normed[start]
        experiment_segment = cumsum_coverage_experiment_normed[stop] - cumsum_coverage_experiment_normed[start]
        difference += abs(control_segment - experiment_segment)
    return difference / len(control_profile)
