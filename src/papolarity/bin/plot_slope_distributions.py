import argparse
import os.path
import seaborn as sns
import matplotlib.pyplot as plt
from ..dto.coverage_comparison_stats import CoverageComparisonStats

def configure_argparser(argparser=None):
    if not argparser:
        argparser = argparse.ArgumentParser(prog="plot_slope_distributions", description="Plot distributions of polarities and slopes over genes")
    argparser.add_argument('coverages_properties', metavar='coverage_properties.tsv', help='Coverage comparison statistics')
    return argparser

def main():
    argparser = configure_argparser()
    args = argparser.parse_args()
    invoke(args)

def invoke(args):
    filename = args.coverages_properties
    filename_wo_ext = os.path.splitext(filename)[0]
    dirname = os.path.dirname(filename_wo_ext)
    basename = os.path.basename(filename_wo_ext)

    slope_data = list(CoverageComparisonStats.each_in_file(filename))
    slope_data = [info for info in slope_data  if info.geom_mean_coverage() >= 10]
    slope_data = list(CoverageComparisonStats.choose_best_transcript(slope_data))
    print(len(slope_data))

    slopes = [info.slope for info in slope_data]
    slopes_multipoint = [info.multipoint_slope for info in slope_data]
    profile_difference = [info.profile_difference for info in slope_data]
    polarities_control = [info.control_polarity_score for info in slope_data]
    polarities_experiment = [info.experiment_polarity_score for info in slope_data]
    polarity_deltas = [info.experiment_polarity_score - info.control_polarity_score for info in slope_data]

    plt.figure()
    plt.title(basename)
    sns.kdeplot(slopes)
    # plt.xlim(-0.02,0.02)
    plt.axvline(x=0, color='green')
    plt.savefig(f'{dirname}/slope_{basename}.png')

    plt.figure()
    plt.title(basename)
    sns.kdeplot(slopes_multipoint)
    # plt.xlim(-0.02,0.02)
    plt.axvline(x=0, color='green')
    plt.savefig(f'{dirname}/slope_multipoint_{basename}.png')

    plt.figure()
    plt.title(basename)
    sns.kdeplot(profile_difference)
    # plt.xlim(-0.02,0.02)
    plt.axvline(x=0, color='green')
    plt.savefig(f'{dirname}/profile_diffs_{basename}.png')

    plt.figure()
    plt.title(basename)
    sns.kdeplot(polarities_control, color='b', label='control')
    sns.kdeplot(polarities_experiment, color='r', label='experiment')
    # plt.xlim(-0.02,0.02)
    plt.axvline(x=0, color='green')
    plt.savefig(f'{dirname}/polarities_{basename}.png')

    plt.figure()
    plt.title(basename)
    sns.kdeplot(polarity_deltas)
    # plt.xlim(-0.02,0.02)
    plt.axvline(x=0, color='green')
    plt.savefig(f'{dirname}/polarity_deltas_{basename}.png')
