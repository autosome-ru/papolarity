import sys
import sklearn.linear_model
from annotation import load_transcript_cds_info
from coverage_profile import transcript_coverages_in_file
from pooling import starjoin_sorted, pooling
import pasio_wrapper
import numpy as np
import math

def slope_by_profiles(cds_profile_1, cds_profile_2, pooled_cds_coverage, splitter):
    total_coverage_1 = sum(cds_profile_1)
    total_coverage_2 = sum(cds_profile_2)
    model = sklearn.linear_model.LinearRegression()
    xs = []
    ys = []
    sample_weights = []
    for (start, stop, _mean_count) in pasio_wrapper.stabile_segments(pooled_cds_coverage, splitter):
        mean_1 = np.mean(cds_profile_1[start:stop])
        mean_2 = np.mean(cds_profile_2[start:stop])
        coeff = ((mean_2 + 1) / (mean_1 + 1)) * ((total_coverage_1 + 1) / (total_coverage_2 + 1))
        coeff = math.log(coeff)
        xs.append([(start + stop) / 2])
        ys.append(coeff)
        sample_weights.append(1)
        # sample_weights.append(stop - start + 1)
    model.fit(xs, ys, sample_weight=sample_weights)
    slope = model.coef_[0]
    return slope

splitter = pasio_wrapper.pasio_splitter()

cds_annotation_fn = sys.argv[1] # 'gencode.vM22.cds_features.tsv'
cds_info_by_transcript = load_transcript_cds_info(cds_annotation_fn)

alignment_1_fn = sys.argv[2]
alignment_2_fn = sys.argv[3]
coverages_1_iter = transcript_coverages_in_file(alignment_1_fn, cds_info_by_transcript, sort_transcripts=True)
coverages_2_iter = transcript_coverages_in_file(alignment_2_fn, cds_info_by_transcript, sort_transcripts=True)

for (transcript_info, (_, coverage_1), (_, coverage_2)) in starjoin_sorted(coverages_1_iter, coverages_2_iter, key=lambda txinfo, coverage: txinfo):
    # print(transcript_info)
    cds_profile_1 = transcript_info.cds_profile(coverage_1)
    cds_profile_2 = transcript_info.cds_profile(coverage_2)

    pooled_cds_coverage = pooling([cds_profile_1, cds_profile_2])

    try:
        slope = slope_by_profiles(cds_profile_1, cds_profile_2, pooled_cds_coverage, splitter)
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
    # coverage_diffs.append( CoverageDifference(*field_values) )
    # print(f'Num errors: {num_errors}', file = sys.stderr)
    # return coverage_diffs





# for (transcript_info, profile) in :
#     print(transcript_info, ','.join(map(str, cds_profile)), sep='\t')
