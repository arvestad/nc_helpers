import argparse
import json
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
import sys

description_text = '''
Compare two clusters stored in a JSON format.

Author: Jon Bryder
'''

def argparser():
    '''
    Set up simple program arguments so we can work on the command line.
    '''
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description=description_text)
    ap.add_argument('jsonfile1', type=argparse.FileType('r'),
                    help='Path to JSON file containing clusters.')
    ap.add_argument('jsonfile2', type=argparse.FileType('r'),
                    help='Another path to a JSON file containing clusters.')
    return ap


def load_json(handle) -> dict | None:
    """Load and return a JSON file, or None on failure."""
    try:
        jsondata = json.load(handle)
        if not jsondata:        # In case empty file
            sys.exit(f'No data in {handle.name}', file=sys.stderr)
        return jsondata
    except:
        sys.exit(f'Could not load JSON data from {handle.name}')


def strip_species_prefix(seq_id: str) -> str:
    """
    Strip species prefix from sequence IDs if present.
    e.g. 'Acamar.WP_012160728.1' -> 'WP_012160728.1'

    Detects the Species.ID format by checking if the part before
    the first dot looks like a species code (letters only, no digits).
    """
    parts = seq_id.split(".", 1)
    if len(parts) == 2 and parts[0].isalpha():
        return parts[1]
    return seq_id


def invert_clusters(data: dict) -> dict[str, str]:
    """
    Convert cluster-centric format to sequence-centric format.

    Handles two formats automatically:
      - Daniel/UniRef format: {'components': [[seq1, seq2], [seq3]], ...metadata...}
      - Orthogroups format:   {'cluster_id': [seq1, seq2], ...}

    Also auto-detects and strips species prefixes (e.g. 'Acamar.WP_xxx')
    so files with different ID formats can be compared directly.

    Returns: {sequence_id: cluster_label}
    """
    # Daniel/UniRef format: clusters live under a 'components' key
    if "components" in data and isinstance(data["components"], list):
        clusters = data["components"]
        return {
            strip_species_prefix(seq): str(i)
            for i, cluster in enumerate(clusters)
            for seq in cluster
        }

    # Orthogroups format: direct {cluster_id: [seq_ids]}
    return {
        strip_species_prefix(seq): cluster_id
        for cluster_id, members in data.items()
        if isinstance(members, list)
        for seq in members
    }



def compare_clusterings(inv1: dict, inv2: dict) -> None:
    """
    Compare two sequence->cluster mappings and print ARS and NMI scores.
    Only sequences present in both files are compared.
    """
    keys1, keys2 = set(inv1.keys()), set(inv2.keys())
    common = keys1 & keys2
    only_in_1 = keys1 - keys2
    only_in_2 = keys2 - keys1

    if only_in_1:
        print(f"  Note: {len(only_in_1)} sequences only in file 1 (excluded)", file=sys.stderr)
    if only_in_2:
        print(f"  Note: {len(only_in_2)} sequences only in file 2 (excluded)", file=sys.stderr)

    if not common:
        sys.exit("Error: No sequences in common between the two files.", file=sys.stderr)
    else:
        print(f"  Sequences in common: {len(common)}", file=sys.stderr)

    # Build parallel label lists sorted for reproducibility
    sorted_seqs = sorted(common)
    labels1 = [inv1[s] for s in sorted_seqs]
    labels2 = [inv2[s] for s in sorted_seqs]

    ars = adjusted_rand_score(labels1, labels2)
    nmi = normalized_mutual_info_score(labels1, labels2)

    print(f"\nAdjusted Rand Score   : {ars:.4f}")
    print(f"Normalized Mutual Info: {nmi:.4f}")


def main():
    ap = argparser()
    args = ap.parse_args()
    
    json1 = load_json(args.jsonfile1)
    json2 = load_json(args.jsonfile2)

    print("\nProcessing file 1...", file=sys.stderr)
    inv1 = invert_clusters(json1)
    print(f"  {len(inv1)} sequences found", file=sys.stderr)

    print("Processing file 2...", file=sys.stderr)
    inv2 = invert_clusters(json2)
    print(f"  {len(inv2)} sequences found", file=sys.stderr)

    print("\nComparing clusterings...", file=sys.stderr)
    compare_clusterings(inv1, inv2)


if __name__ == "__main__":
    main()
