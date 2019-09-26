import sys
import sklearn.linear_model
from annotation import load_transcript_cds_info
from coverage_profile import transcript_coverages_in_file
from pooling import starjoin_sorted, pooling
import pasio_wrapper
import numpy as np
import math
from pybedtools import BedTool

def stabilize_profile(profile, segments):
    stable_profile = np.zeros_like(profile)
    for (start, stop, *_) in segments:
        stable_profile[start:stop] = np.mean(profile[start:stop])
    return stable_profile

def slope_by_profiles(control_profile, experiment_profile, segments):
    total_coverage_control = sum(control_profile)
    total_coverage_experiment = sum(experiment_profile)
    model = sklearn.linear_model.LinearRegression()
    xs = []
    ys = []
    sample_weights = []
    num_segments = len(segments)
    for (start, stop, *_) in segments:
        mean_control = np.mean(control_profile[start:stop])
        mean_experiment = np.mean(experiment_profile[start:stop])
        normalized_mean_control = (mean_control + 1) / (total_coverage_control + num_segments)
        normalized_mean_experiments = (mean_experiments + 1) / (total_coverage_experiments + num_segments)
        detrended_profile = normalized_mean_experiments / normalized_mean_control
        log_detrended_profile = math.log(detrended_profile)
        xs.append([(start + stop) / 2])
        ys.append(log_detrended_profile)
        sample_weights.append(1)
        # sample_weights.append(stop - start)
    model.fit(xs, ys, sample_weight=sample_weights)
    slope = model.coef_[0]
    return slope

splitter = pasio_wrapper.pasio_splitter()

cds_annotation_fn = sys.argv[1] # 'gencode.vM22.cds_features.tsv'
cds_info_by_transcript = load_transcript_cds_info(cds_annotation_fn)

alignment_1_fn = sys.argv[2]
alignment_2_fn = sys.argv[3]

alignment_1 = BedTool(alignment_1_fn)
alignment_2 = BedTool(alignment_2_fn)

coverages_1_iter = transcript_coverages_in_file(alignment_1, cds_info_by_transcript, sort_transcripts=True)
coverages_2_iter = transcript_coverages_in_file(alignment_2, cds_info_by_transcript, sort_transcripts=True)

for (transcript_info, (_, coverage_1), (_, coverage_2)) in starjoin_sorted(coverages_1_iter, coverages_2_iter, key=lambda txinfo, coverage: txinfo):
    cds_profile_1 = transcript_info.cds_profile(coverage_1)
    cds_profile_2 = transcript_info.cds_profile(coverage_2)

    pooled_cds_coverage = pooling([cds_profile_1, cds_profile_2])
    segments = list( pasio_wrapper.stable_segments(pooled_cds_coverage, splitter) ) # [(start, stop, lambda), ...]
    stable_cds_profile_1 = stabilize_profile(cds_profile_1, segments)
    stable_cds_profile_2 = stabilize_profile(cds_profile_2, segments)

    try:
        slope = slope_by_profiles(cds_profile_1, cds_profile_2, segments)
    except:
        # num_errors += 1
        print('Error', file = sys.stderr)
        traceback.print_exc(file = sys.stderr)
        continue

    
    field_values = [
        # transcript_info_1.gene_id, transcript_id,
        # transcript_info_1.transcript_length, transcript_info_1.cds_start, transcript_info_1.cds_stop,
        # transcript_info_1.read_number, transcript_info_1.mean_cds_coverage, transcript_info_1.polarity_score,
        # transcript_info_2.read_number, transcript_info_2.mean_cds_coverage, transcript_info_2.polarity_score,
        # slope,
        transcript_info, slope,
    ]
    print(*field_values, sep='\t')
    print(*stable_cds_profile_1, sep='\t')
    print(*cds_profile_1, sep='\t')
    print(*cds_profile_2, sep='\t')
    print(*stable_cds_profile_2, sep='\t')
    # coverage_diffs.append( CoverageDifference(*field_values) )
    # print(f'Num errors: {num_errors}', file = sys.stderr)
    # return coverage_diffs





# for (transcript_info, profile) in :
#     print(transcript_info, ','.join(map(str, cds_profile)), sep='\t')
