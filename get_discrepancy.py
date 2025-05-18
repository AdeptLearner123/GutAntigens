from constants import *
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():
    my_sequences = set()
    bio_sequences = set()
    
    for record in tqdm(SeqIO.parse(COMPLETE_SEQUENCES, "fasta")):
        bio_sequences.add(str(record.seq))
    
    with open("genomes/results/unique_complete_sequences") as f:
        my_sequences.update(f.read().splitlines())
    
    difference = list(my_sequences - bio_sequences)
    print(f"Diff: {len(difference)}")
    for seq in difference:
        print(seq)


if __name__ == "__main__":
    main()