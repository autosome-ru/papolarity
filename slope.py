import sys
import sklearn.linear_model
from annotation import load_transcript_cds_info
from coverage_profile import transcript_coverages_in_file
from pooling import starjoin_sorted, pooling
import pasio_wrapper
from polarity_score import polarity_score
import numpy as np
import math
from pybedtools import BedTool

import traceback

def stabilize_profile(profile, segments):
    stable_profile = np.zeros_like(profile)
    for (start, stop, *_) in segments:
        stable_profile[start:stop] = np.mean(profile[start:stop])
    return stable_profile

def slope_by_profiles(control_profile, experiment_profile, segments):
    total_coverage_control = sum(control_profile)
    total_coverage_experiment = sum(experiment_profile)
    assert len(control_profile) == len(experiment_profile)
    profile_len = len(control_profile)
    model = sklearn.linear_model.LinearRegression()
    xs = []
    ys = []
    sample_weights = []
    num_segments = len(segments)
    for (start, stop, *_) in segments:
        mean_control = np.mean(control_profile[start:stop])
        mean_experiment = np.mean(experiment_profile[start:stop])
        normalized_mean_control = (mean_control + 1) / (total_coverage_control + num_segments)
        normalized_mean_experiment = (mean_experiment + 1) / (total_coverage_experiment + num_segments)
        detrended_profile = normalized_mean_experiment / normalized_mean_control
        log_detrended_profile = math.log(detrended_profile)
        coord = (start + stop - 1) / 2
        rel_coord = coord / profile_len
        xs.append(rel_coord)
        ys.append(log_detrended_profile)
        sample_weights.append(1)
        # sample_weights.append(stop - start)
    xs = np.array(xs)
    model.fit(xs.reshape(-1, 1), ys, sample_weight=sample_weights)
    slope = model.coef_[0]
    return slope

def compare_profiles(transcript_info, coverage_control, coverage_experiment, drop_start_stop=True):
    cds_profile_control = transcript_info.cds_profile(coverage_control)
    cds_profile_experiment = transcript_info.cds_profile(coverage_experiment)

    # start and stop codon can have piles of reads, so we usually want to drop them
    if drop_start_stop:
        flank = 15  # half-ribosome (~half of riboseq footprint length)
        cds_profile_control = cds_profile_control[flank:-flank]
        cds_profile_experiment = cds_profile_experiment[flank:-flank]
    pooled_cds_coverage = pooling([cds_profile_control, cds_profile_experiment])
    if len(pooled_cds_coverage) == 0:
        return None
    segments = list( pasio_wrapper.stable_segments(pooled_cds_coverage, splitter) ) # [(start, stop, lambda), ...]
    stable_cds_profile_control = stabilize_profile(cds_profile_control, segments)
    stable_cds_profile_experiment = stabilize_profile(cds_profile_experiment, segments)

    slope = slope_by_profiles(cds_profile_control, cds_profile_experiment, segments)

    info = {
        'transcript_info': transcript_info,
        'slope': slope,
        'control_cds_mean': np.mean(cds_profile_control), 'experiment_cds_mean': np.mean(cds_profile_experiment),
        'polarity_score_control': polarity_score(cds_profile_control), 'polarity_score_experiment': polarity_score(cds_profile_experiment),
    }
    return info

def compare_multiple_alignments(cds_info_by_transcript, alignment_control, alignment_experiment):
    coverages_control_iter = transcript_coverages_in_file(alignment_control, sort_transcripts=True)
    coverages_experiment_iter = transcript_coverages_in_file(alignment_experiment, sort_transcripts=True)

    for (transcript_id, (_, coverage_control), (_, coverage_experiment)) in starjoin_sorted(coverages_control_iter, coverages_experiment_iter, key=lambda txid, coverage: txid):
        if transcript_id not in cds_info_by_transcript:
            continue
        transcript_info = cds_info_by_transcript[transcript_id]
        info = compare_profiles(transcript_info, coverage_control, coverage_experiment, drop_start_stop=True)
        if info:
            yield info

splitter = pasio_wrapper.pasio_splitter()

cds_annotation_fn = sys.argv[1] # 'gencode.vM22.cds_features.tsv'
cds_info_by_transcript = load_transcript_cds_info(cds_annotation_fn)

alignment_control_fn = sys.argv[2]
alignment_experiment_fn = sys.argv[3]

alignment_control = BedTool(alignment_control_fn)
alignment_experiment = BedTool(alignment_experiment_fn)

header = [
    'gene_id', 'transcript_id', 'transcript_length', 'cds_start', 'cds_stop',
    'slope',
    'control_cds_mean', 'experiment_cds_mean',
    'polarity_score_control', 'polarity_score_experiment',
]

print('\t'.join(header))

for info in compare_multiple_alignments(cds_info_by_transcript, alignment_control, alignment_experiment):
    field_values = [
        info['transcript_info'],
        info['slope'],
        info['control_cds_mean'], info['experiment_cds_mean'],
        info['polarity_score_control'], info['polarity_score_experiment'],
    ]
    print(*field_values, sep='\t')
