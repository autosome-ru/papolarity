from coding_transcript_info import CodingTranscriptInfo
import itertools
import sys
import seaborn as sns
import matplotlib.pyplot as plt

def read_slope_data(filename):
    with open(filename) as f:
        header = f.readline()
        for line in f:
            row = line.rstrip('\n').split('\t')
            gene_id, transcript_id, transcript_length, cds_start, cds_stop, \
                slope, \
                control_mean_coverage, experiment_mean_coverage, \
                control_polarity_score, experiment_polarity_score = row

            info = {
                'transcript_info': CodingTranscriptInfo(gene_id, transcript_id, int(transcript_length), int(cds_start), int(cds_stop)),
                'slope': float(slope),
                'control_mean_coverage': float(control_mean_coverage), 'experiment_mean_coverage': float(experiment_mean_coverage),
                'control_polarity_score': float(control_polarity_score), 'experiment_polarity_score': float(experiment_polarity_score),
            }
            yield info


# by default takes max by geometric mean of cds-mean between experiments
def choose_best_transcript(slope_data, key=lambda info: info['control_mean_coverage'] * info['experiment_mean_coverage']):
    slope_data = sorted(slope_data, key=lambda info: info['transcript_info'].gene_id)
    for (gene_id, infos) in itertools.groupby(slope_data, lambda info: info['transcript_info'].gene_id):
        best_transcript_info = max(infos, key=key)
        yield best_transcript_info

filename = sys.argv[1]
slope_data = list(read_slope_data(filename))
slopes_data = [info for info in choose_best_transcript(slope_data)]
slope_data = [info for info in slope_data  if info['control_mean_coverage'] * info['experiment_mean_coverage'] >= 100]

slope_data = list(choose_best_transcript(slope_data))
print(len(slope_data))

slopes = [info['slope'] for info in slope_data]
polarities_control = [info['control_polarity_score'] for info in slope_data]
polarities_experiment = [info['experiment_polarity_score'] for info in slope_data]
polarity_deltas = [info['experiment_polarity_score'] - info['control_polarity_score'] for info in slope_data]

plt.figure()
sns.kdeplot(slopes)
# plt.xlim(-0.02,0.02)
plt.axvline(x=0, color='green')
plt.savefig('slope.png')

plt.figure()
sns.kdeplot(polarities_control, color='b')
sns.kdeplot(polarities_experiment, color='r')
# plt.xlim(-0.02,0.02)
plt.axvline(x=0, color='green')
plt.savefig('polarities.png')

plt.figure()
sns.kdeplot(polarity_deltas)
# plt.xlim(-0.02,0.02)
plt.axvline(x=0, color='green')
plt.savefig('polarity_deltas.png')

plt.figure()
sns.scatterplot(polarity_deltas, slopes)
plt.savefig('polarity_delta_vs_slope.png')
