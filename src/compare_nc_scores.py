import argparse
import matplotlib.pyplot as plt
import seaborn as sns

def argparser():
    '''
    Set up simple program arguments so we can work on the command line.
    '''
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description='Compare NC scores computed in different ways.')
    ap.add_argument('nc1', type=argparse.FileType('r'),
                    help='Path to JSON file containing clusters.')
    ap.add_argument('nc2', type=argparse.FileType('r'),
                    help='Another path to a JSON file containing clusters.')
    ap.add_argument('outfile', type=str,
                    help='Where to put the output PDF')
    return ap


def read_nc_scores(scorefile):
    scores = dict()
    for line in scorefile:
        cols = line.split()
        if len(cols) != 3:
            raise ValueError(f'Expected three columns, but got "{line}".')
        left, right, val = cols
        if left < right:
            scores[left, right] = float(val)
        elif right < left:
            scores[right, left] = float(val)
    return scores


def plot_scores(nc1, nc2):
    s1 = set(nc1.keys())
    s2 = set(nc2.keys())
    
    joint = s1 & s2
    unique_s1 = s1 - s2
    unique_s2 = s2 - s1

    x = list()
    y = list()
    for pair in joint:
        x.append(nc1[pair])
        y.append(nc2[pair])

    for pair in unique_s1:
        x.append(nc1[pair])
        y.append(0.0)

    for pair in unique_s2:
        x.append(0.0)
        y.append(nc2[pair])

    return sns.jointplot(x=x, y=y, kind='scatter', alpha=0.3, color='blue')


def main():
    ap = argparser()
    args = ap.parse_args()

    nc1 = read_nc_scores(args.nc1)
    nc2 = read_nc_scores(args.nc2)

    fig = plot_scores(nc1, nc2)
    fig.ax_joint.set_xlabel(args.nc1)
    fig.ax_joint.set_ylabel(args.nc2)
    plt.savefig(args.outfile)
    
    
    


if __name__ == "__main__":
    main()

