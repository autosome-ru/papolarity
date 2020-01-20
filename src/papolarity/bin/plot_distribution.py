import argparse
import seaborn as sns
import matplotlib.pyplot as plt
from ..tsv_reader import each_in_tsv
from ..utils import drop_none

def configure_argparser(argparser=None):
    if not argparser:
        argparser = argparse.ArgumentParser(prog="plot_distribution", description="Plot distributions of a feature or several features")
    argparser.add_argument('table', metavar='table.tsv', help='Table in tab-separated format')
    argparser.add_argument('--fields', nargs='*', required=True, help='Fields to plot on the same figure')
    argparser.add_argument('--labels', nargs='*', help='Set legend labels (by default field names are used)')
    argparser.add_argument('--title', help="Add plot title")
    argparser.add_argument('--xlim', nargs=2, type=float, help="Add limits for X-axis (in form `--xlim min max`)")
    argparser.add_argument('--ylim', nargs=2, type=float, help="Add limits for Y-axis (in form `--ylim min max`)")
    argparser.add_argument('--zero-line', help="Add vertical line at zero of specified color")

    has_legend_group = argparser.add_mutually_exclusive_group(required=True)
    has_legend_group.add_argument('--legend', action='store_true', dest='has_legend', help="Table has legend")
    has_legend_group.add_argument('--no-legend', action='store_false', dest='has_legend', help="Tables doesn't have legend")
    
    argparser.add_argument('--display', action='store_true', help="Display resulting figure")
    argparser.add_argument('--output-file', '-o', dest='output_file', help="Store results at this path")
    return argparser

def main():
    argparser = configure_argparser()
    args = argparser.parse_args()
    invoke(args)

def invoke(args):
    dtype = lambda x: float(x) if x != '' else None
    data = list(each_in_tsv(args.table))

    fields = args.fields

    if len(fields) == 0:
        raise ValueError(f'Specify at least one field to plot')

    if args.labels:
        if len(fields) != len(args.labels):
            raise ValueError(f'Number of labels ({len(args.labels)}) should be equal to number of fields ({len(fields)})')
        labels = args.labels
    else:
        labels = fields

    for field in fields:
        for row in data:
            row[field] = dtype(row[field])

    plt.figure()

    if args.title:
        plt.title(args.title)
    
    for (field, label) in zip(fields, labels):
        values = [row[field] for row in data]
        sns.kdeplot(drop_none(values), label = label, gridsize=10000, clip=args.xlim,)

    if args.xlim:
        plt.xlim(*args.xlim)

    if args.ylim:
        plt.ylim(*args.ylim)

    if args.zero_line:
        plt.axvline(x=0, color=args.zero_line)

    if args.output_file:
        plt.savefig(args.output_file)

    if args.display:
        plt.show()