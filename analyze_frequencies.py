from collections import Counter
from tqdm import tqdm
from constants import *
from Bio import SeqIO

def main():
    seq_samples = {}
    frequencies = Counter()

    for record in tqdm(SeqIO.parse(COMPLETE_SEQUENCES, "fasta")):
        seq = str(record.seq)
        samples = record.description.split(" # ")[1].split(";")
        seq_samples[seq] = samples
        frequencies[len(samples)] += 1
    
    sorted_freqs = sorted(frequencies.items(), key=lambda item: item[1], reverse=True)
    print("Sorted freqs")
    print(sorted_freqs)

    cluster_freqs = Counter()
    with open("genomes/clusters/cluster99.clstr") as f:
        cluster_samples = set()
        for line in tqdm(f):
            if line.startswith(">"):
                if len(cluster_samples) > 0:
                    cluster_freqs[len(cluster_samples)] += 1
                cluster_samples = set()
            else:
                seq = line.split("\t")[1]
                cluster_samples.update(seq_samples[seq])
    
    print("Sorted cluster freqs")
    sorted_cluster_freqs = sorted(frequencies.items(), key=lambda item: item[1], reverse=True)
    print(sorted_cluster_freqs)


if __name__ == "__main__":
    main()