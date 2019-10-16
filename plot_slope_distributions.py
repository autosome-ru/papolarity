from slope import CoverageComparisonStats, choose_best_transcript
import sys
import seaborn as sns
import matplotlib.pyplot as plt

filename = sys.argv[1]
slope_data = list(CoverageComparisonStats.each_in_file(filename))
slope_data = [info for info in slope_data  if info.geom_mean_coverage() >= 10]
slope_data = list(choose_best_transcript(slope_data))
print(len(slope_data))

slopes = [info.slope for info in slope_data]
polarities_control = [info.control_polarity_score for info in slope_data]
polarities_experiment = [info.experiment_polarity_score for info in slope_data]
polarity_deltas = [info.experiment_polarity_score - info.control_polarity_score for info in slope_data]

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
