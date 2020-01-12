from sklearn.linear_model import LinearRegression
import numpy as np
from math import log

def slope_by_profiles(control_profile, experiment_profile, segmentation, mode='center'):
    total_coverage_control = sum(control_profile)
    total_coverage_experiment = sum(experiment_profile)
    assert len(control_profile) == len(experiment_profile)
    assert len(control_profile) == segmentation.segmentation_length
    profile_len = len(control_profile)
    model = LinearRegression()
    xs = []
    ys = []
    sample_weights = []
    num_segments = segmentation.num_segments
    for segment in segmentation.segments:
        start = segment.start
        stop = segment.stop
        mean_control = np.mean(control_profile[start:stop])
        mean_experiment = np.mean(experiment_profile[start:stop])
        normalized_mean_control = (mean_control + 1) / (total_coverage_control + num_segments)
        normalized_mean_experiment = (mean_experiment + 1) / (total_coverage_experiment + num_segments)
        detrended_profile = normalized_mean_experiment / normalized_mean_control
        log_detrended_profile = log(detrended_profile)
        if mode == 'center':
            coord = (start + stop - 1) / 2
            rel_coord = coord / profile_len
            xs.append(rel_coord)
            ys.append(log_detrended_profile)
            sample_weights.append(1)
        elif mode == 'weighted_center':
            coord = (start + stop - 1) / 2
            rel_coord = coord / profile_len
            weight = stop - start
            xs.append(rel_coord)
            ys.append(log_detrended_profile)
            sample_weights.append(weight)
        elif mode == 'every_point':
            for pos in range(start, stop):
                rel_coord = pos / profile_len
                xs.append(rel_coord)
                ys.append(log_detrended_profile)
                sample_weights.append(1)
        else:
            raise NotImplementedError()
        # sample_weights.append(stop - start)
    xs = np.array(xs)
    model.fit(xs.reshape(-1, 1), ys, sample_weight=sample_weights)
    slope = model.coef_[0]
    return slope

def profile_difference(control_profile, experiment_profile, segmentation):
    assert len(control_profile) == len(experiment_profile)
    control_profile_sum = np.sum(control_profile)
    experiment_profile_sum = np.sum(experiment_profile)
    normed_control_profile = control_profile / control_profile_sum  if control_profile_sum != 0  else control_profile
    normed_experiment_profile = experiment_profile / experiment_profile_sum  if experiment_profile_sum != 0  else experiment_profile
    difference = 0
    for segment in segmentation.segments:
        start = segment.start
        stop = segment.stop
        mean_control = np.mean(normed_control_profile[start:stop])
        mean_experiment = np.mean(normed_experiment_profile[start:stop])
        difference += (stop - start) * abs(mean_control - mean_experiment)
    return difference / len(control_profile)
