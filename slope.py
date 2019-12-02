from coding_transcript_info import CodingTranscriptInfo
import itertools
import sys
from collections import namedtuple
from sklearn.linear_model import LinearRegression
from polarity_score import polarity_score
import pasio
from coverage_profile import transcript_coverages_from_alignment
from pooling import starjoin_sorted, pooling
import numpy as np
from math import log
from dto.coverage_comparison_stats import CoverageComparisonStats

import logging
logger = logging.getLogger('pasio')
logger.setLevel(logging.WARNING)

class TranscriptComparator:
    # start and stop codon can have piles of reads, so we usually want to drop them
    # flank lengths to drop are of 15nt ~= half-ribosome (~half of riboseq footprint length)
    def __init__(self, cds_info_by_transcript, splitter, drop_start_flank=15, drop_stop_flank=15):
        self.cds_info_by_transcript = cds_info_by_transcript
        self.splitter = splitter
        self.drop_start_flank = drop_start_flank
        self.drop_stop_flank = drop_stop_flank

    def compare_multiple_alignments(self, alignment_control, alignment_experiment):
        coverages_control_iter = transcript_coverages_from_alignment(alignment_control, sort_transcripts=True)
        coverages_experiment_iter = transcript_coverages_from_alignment(alignment_experiment, sort_transcripts=True)

        for (transcript_id, (_, coverage_control), (_, coverage_experiment)) in starjoin_sorted(coverages_control_iter, coverages_experiment_iter, key=lambda txid, coverage: txid):
            if transcript_id not in self.cds_info_by_transcript:
                continue
            transcript_info = self.cds_info_by_transcript[transcript_id]
            info = self.compare_profiles(transcript_info, coverage_control, coverage_experiment)
            if info:
                yield info

    def compare_profiles(self, transcript_info, coverage_control, coverage_experiment):
        cds_profile_control = transcript_info.cds_profile(coverage_control)
        cds_profile_experiment = transcript_info.cds_profile(coverage_experiment)

        # start and stop codon can have piles of reads, so we usually want to drop them
        cds_profile_control = cds_profile_control[self.drop_start_flank : -self.drop_stop_flank]
        cds_profile_experiment = cds_profile_experiment[self.drop_start_flank : -self.drop_stop_flank]
        pooled_cds_coverage = pooling([cds_profile_control, cds_profile_experiment])
        if len(pooled_cds_coverage) == 0:
            return None

        segments = list(pasio.segments_with_scores(pooled_cds_coverage, self.splitter)) # [(start, stop, lambda), ...]
        # stable_cds_profile_control = stabilize_profile(cds_profile_control, segments)
        # stable_cds_profile_experiment = stabilize_profile(cds_profile_experiment, segments)
        return CoverageComparisonStats.make_from_profiles(transcript_info, cds_profile_control, cds_profile_experiment, segments)


def slope_by_profiles(control_profile, experiment_profile, segments, mode='center'):
    total_coverage_control = sum(control_profile)
    total_coverage_experiment = sum(experiment_profile)
    assert len(control_profile) == len(experiment_profile)
    profile_len = len(control_profile)
    model = LinearRegression()
    xs = []
    ys = []
    sample_weights = []
    num_segments = len(segments)
    for segment in segments:
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


def profile_difference(control_profile, experiment_profile, segments):
    assert len(control_profile) == len(experiment_profile)
    control_profile_sum = np.sum(control_profile)
    experiment_profile_sum = np.sum(experiment_profile)
    normed_control_profile = control_profile / control_profile_sum  if control_profile_sum != 0  else control_profile
    normed_experiment_profile = experiment_profile / experiment_profile_sum  if experiment_profile_sum != 0  else experiment_profile
    difference = 0
    for segment in segments:
        start = segment.start
        stop = segment.stop
        mean_control = np.mean(normed_control_profile[start:stop])
        mean_experiment = np.mean(normed_experiment_profile[start:stop])
        difference += (stop - start) * abs(mean_control - mean_experiment)
    return difference / len(control_profile)

# by default takes max by geometric mean of cds-mean between experiments
def choose_best_transcript(slope_data, key=lambda info: info.geom_mean_coverage()):
    by_gene = lambda info: info.gene_id
    slope_data = sorted(slope_data, key=by_gene)
    for (gene_id, infos) in itertools.groupby(slope_data, by_gene):
        best_transcript_info = max(infos, key=key)
        yield best_transcript_info

def stabilize_profile(profile, segments):
    stable_profile = np.zeros_like(profile)
    for segment in segments:
        start = segment.start
        stop = segment.stop
        stable_profile[start:stop] = np.mean(profile[start:stop])
    return stable_profile
