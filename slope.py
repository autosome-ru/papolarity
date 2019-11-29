from coding_transcript_info import CodingTranscriptInfo
import itertools
import sys
from collections import namedtuple
from sklearn.linear_model import LinearRegression
from polarity_score import polarity_score
import pasio
from coverage_profile import transcript_coverages_in_file
from pooling import starjoin_sorted, pooling
import numpy as np
from math import log

import logging
logger = logging.getLogger('pasio')
logger.setLevel(logging.WARNING)


coverage_stats_fields = [
    'transcript_info', 'slope',
    'control_mean_coverage', 'experiment_mean_coverage',
    'control_total_coverage', 'experiment_total_coverage',
    'control_polarity_score', 'experiment_polarity_score',
    'num_segments',
    'multipoint_slope', 'profile_difference',
]
class CoverageComparisonStats(namedtuple('CoverageComparisonStats', coverage_stats_fields)):
    @property
    def gene_id(self):
        return self.transcript_info.gene_id

    @property
    def transcript_id(self):
        return self.transcript_info.transcript_id

    @property
    def transcript_length(self):
        return self.transcript_info.transcript_length
    
    @property
    def cds_start(self):
        return self.transcript_info.cds_start

    @property
    def cds_stop(self):
        return self.transcript_info.cds_stop

    def __repr__(self):
        fields = [
            self.gene_id, self.transcript_id, self.transcript_length, self.cds_start, self.cds_stop,
            self.slope,
            self.control_mean_coverage, self.experiment_mean_coverage,
            self.control_total_coverage, self.experiment_total_coverage,
            self.control_polarity_score, self.experiment_polarity_score,
            self.num_segments,
            self.multipoint_slope, self.profile_difference,
        ]
        return '\t'.join(map(str, fields))

    @classmethod
    def header(cls):
        fields = [
            'gene_id', 'transcript_id', 'transcript_length', 'cds_start', 'cds_stop',
            'slope',
            'control_mean_coverage', 'experiment_mean_coverage',
            'control_total_coverage', 'experiment_total_coverage',
            'control_polarity_score', 'experiment_polarity_score',
            'num_segments',
            'multipoint_slope', 'profile_difference',
        ]
        return '\t'.join(fields)

    @classmethod
    def from_string(cls, line):
        row = line.rstrip('\n').split('\t')
        gene_id, transcript_id, transcript_length, cds_start, cds_stop, \
            slope, \
            control_mean_coverage, experiment_mean_coverage, \
            control_total_coverage, experiment_total_coverage, \
            control_polarity_score, experiment_polarity_score, \
            num_segments, \
            multipoint_slope, profile_difference, \
            *rest = row
        info = {
            'transcript_info': CodingTranscriptInfo(gene_id, transcript_id, int(transcript_length), int(cds_start), int(cds_stop)),
            'slope': float(slope),
            'control_mean_coverage': float(control_mean_coverage), 'experiment_mean_coverage': float(experiment_mean_coverage),
            'control_total_coverage': int(control_total_coverage), 'experiment_total_coverage': int(experiment_total_coverage),
            'control_polarity_score': float(control_polarity_score), 'experiment_polarity_score': float(experiment_polarity_score),
            'num_segments': int(num_segments),
            'multipoint_slope': float(multipoint_slope),
            'profile_difference': float(profile_difference),
        }
        return cls(**info)

    @classmethod
    def make_from_profiles(cls, transcript_info, cds_profile_control, cds_profile_experiment, segments):
        info = {
            'transcript_info': transcript_info,
            'slope': slope_by_profiles(cds_profile_control, cds_profile_experiment, segments, mode='center'),
            'control_mean_coverage': np.mean(cds_profile_control), 'experiment_mean_coverage': np.mean(cds_profile_experiment),
            'control_total_coverage': np.sum(cds_profile_control), 'experiment_total_coverage': np.sum(cds_profile_experiment),
            'control_polarity_score': polarity_score(cds_profile_control), 'experiment_polarity_score': polarity_score(cds_profile_experiment),
            'num_segments': len(segments),
            'multipoint_slope': slope_by_profiles(cds_profile_control, cds_profile_experiment, segments, mode='every_point'),
            'profile_difference': profile_difference(cds_profile_control, cds_profile_experiment, segments),
        }
        return cls(**info)

    @classmethod
    def each_in_file(cls, filename):
        with open(filename) as f:
            header = f.readline()
            for line in f:
                yield cls.from_string(line)

    @classmethod
    def print(cls, infos, file=sys.stdout, extended=False):
        if extended:
            infos = list(infos)
            best_transcript_ids = set(info.transcript_id for info in choose_best_transcript(infos))

            print(cls.header() + '\t' + '\t'.join(['geom_mean_coverage', 'polarity_difference', 'is_best_transcript']), file=file)
            for info in infos:
                if info.transcript_id in best_transcript_ids:
                    is_best_transcript = '+'
                else:
                    is_best_transcript = '-'
                print(str(info) + '\t' + '\t'.join(map(str, [info.geom_mean_coverage(), info.polarity_difference(), is_best_transcript])), file=file)
        else:
            print(cls.header(), file=file)
            for info in infos:
                print(info, file=file)

    def geom_mean_coverage(self):
        return (self.control_mean_coverage * self.experiment_mean_coverage) ** 0.5

    def polarity_difference(self):
        return self.experiment_polarity_score - self.control_polarity_score

class TranscriptComparator:
    # start and stop codon can have piles of reads, so we usually want to drop them
    # flank lengths to drop are of 15nt ~= half-ribosome (~half of riboseq footprint length)
    def __init__(self, cds_info_by_transcript, splitter, drop_start_flank=15, drop_stop_flank=15):
        self.cds_info_by_transcript = cds_info_by_transcript
        self.splitter = splitter
        self.drop_start_flank = drop_start_flank
        self.drop_stop_flank = drop_stop_flank

    def compare_multiple_alignments(self, alignment_control, alignment_experiment):
        coverages_control_iter = transcript_coverages_in_file(alignment_control, sort_transcripts=True)
        coverages_experiment_iter = transcript_coverages_in_file(alignment_experiment, sort_transcripts=True)

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
