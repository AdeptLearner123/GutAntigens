from collections import Counter
from tqdm import tqdm
from constants import *
from Bio import SeqIO
import math
import matplotlib.pyplot as plt
import numpy as np

def main():
    seq_samples = {}
    frequencies = Counter()

    with open(COMPLETE_SEQUENCES) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            samples = line.split(" # ")[1].split(";")
            seqid = line.split(">")[1].split(" # ")[0]
            seq_samples[seqid] = samples
            idx = math.ceil(math.log10(len(samples)))
            frequencies[idx] += 1
    
    sorted_freqs = sorted(frequencies.items(), key=lambda item: item[0], reverse=True)
    print("Sorted freqs")
    print(sorted_freqs)

    cluster_freqs = Counter()
    with open("genomes/clusters/cluster80.clstr") as f:
        cluster_samples = set()
        for line in tqdm(f):
            if line.startswith(">"):
                if len(cluster_samples) > 0:
                    idx = math.ceil(math.log10(len(cluster_samples)))
                    cluster_freqs[idx] += 1
                cluster_samples = set()
            else:
                seqid = line.split(">")[1].split("...")[0]
                cluster_samples.update(seq_samples[seqid])
    
    print("Sorted cluster freqs")
    sorted_cluster_freqs = sorted(cluster_freqs.items(), key=lambda item: item[0], reverse=True)
    print(sorted_cluster_freqs)


if __name__ == "__main__":
    main()