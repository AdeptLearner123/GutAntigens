from collections import defaultdict
from tqdm import tqdm
from constants import *
import psutil

KMER_SIZE = 3

def build_kmer_table(seqs):
    table = defaultdict(lambda: list())
    for seq in tqdm(seqs):
        kmers = [ seq[i:i+KMER_SIZE] for i in range(len(seq) - KMER_SIZE + 1) ]
        for kmer in kmers:
            table[kmer].append(seq)
        memory = psutil.virtual_memory()
        available_memory = memory.available
        print(f"Available memory: {available_memory / (1024 ** 3):.2f} GB", end="\r")

    return table


def lookup_kmer_table(seq, table):
    kmers = [ seq[i:i+KMER_SIZE] for i in range(len(seq) - KMER_SIZE + 1) ]
    matches = [ set(table[kmer]) for kmer in kmers ]
    return set.intersection(*matches)


def get_contained_seqs():
    with open(UNIQUE_COMPLETE_SEQUENCES) as f:
        complete_seqs = f.read().splitlines()
    
    kmer_table = build_kmer_table(complete_seqs)
    contained_seqs = []
    uncontained_seqs = []

    with open(UNIQUE_PARTIAL_SEQUENCES) as f:
        partial_seqs = f.read().splitlines()
    
    for seq in tqdm(partial_seqs):
        if (len(lookup_kmer_table(seq, kmer_table)) > 0):
            contained_seqs.append(seq)
        else:
            uncontained_seqs.append(seq)
    
    print(f"Contained sequences: {len(contained_seqs)}")
    print(f"Uncontained sequences: {len(uncontained_seqs)}")

    with open(PARTIAL_CONTAINED, "w") as f:
        f.write("\n".join(contained_seqs))


def get_collapsed_seqs():
    with open(UNIQUE_PARTIAL_SEQUENCES) as f:
        partial_seqs = f.read().splitlines()
    with open(PARTIAL_CONTAINED) as f:
        contained_seqs = f.read().splitlines()
    uncontained_seqs = list(set(partial_seqs) - set(contained_seqs))
    kmer_table = build_kmer_table(uncontained_seqs)
    uncollapsed_seqs = []
    collapsed_seqs = []
    
    for seq in tqdm(uncontained_seqs):
        if (len(lookup_kmer_table(seq, kmer_table)) > 1):
            collapsed_seqs.append(seq)
        else:
            uncollapsed_seqs.append(seq)
    
    print(f"Collapsed sequences: {len(collapsed_seqs)}")
    print(f"Uncollapsed sequences: {len(uncollapsed_seqs)}")

    with open(PARTIAL_COLLAPSED, "w") as f:
        f.write("\n".join(collapsed_seqs))
    with open(PARTIAL_UNCOLLAPSED, "w") as f:
        f.write("\n".join(uncollapsed_seqs))


def main():
    get_contained_seqs()
    get_collapsed_seqs()


if __name__ == "__main__":
    main()