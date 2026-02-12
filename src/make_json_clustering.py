import argparse
from enum import IntEnum
import json
import sys

description_text = """
Some clustering tools output clusterings in a line-based format like
  acc   cluster-id
or 
  cluster id   acc. 

This tool converts this plain-text format to a JSON format enabling
straightforward parsing in other tools.
"""

def argparser():
    '''
    Set up simple program arguments so we can work on the command line.
    '''
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description=description_text)
    ap.add_argument('-left', type=argparse.FileType('r'),
                    help='Path to text file with two white-space delimited columns, and the cluster-id in the left column.')
    ap.add_argument('-right', type=argparse.FileType('r'),
                    help='Same as -left, but with the cluster-id in the right column.')
    return ap

class Column(IntEnum):
    LEFT = 0
    RIGHT = 1
    

def convert(handle, cluster_column, seq_column):
    """
    Handle is an open readable file.
    clustercolumn and seq_columns are integers, 0 for left column and 1 for right column.
    """
    d = dict()
    for line in handle:
        line = line.strip()
        columns = line.split()
        family = columns[cluster_column]
        prot = columns[seq_column]
        
        if family in d:
            d[family].append(prot)
        else:
            d[family] = [prot]
    return d



        
def main():
    ap = argparser()
    args = ap.parse_args()
    d = None
    
    if args.left:
        d = convert(args.left, Column.LEFT, Column.RIGHT)
    elif args.right:
        d = convert(args.right, Column.RIGHT, Column.LEFT)
    else:
        print('Wrong usage. Look at the options.', file=sys.stderr)
        sys.exit(1)

    print(json.dumps(d))    


if __name__ == "__main__":
    main()
    
