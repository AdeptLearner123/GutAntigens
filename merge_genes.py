import os
from constants import *
from tqdm import tqdm
from Bio import SeqIO
from collections import Counter

def generate_nt_sequences():
    os.makedirs(NT_SEQUENCES_DIR, exist_ok=True)
    fasta_files = [ file for file in os.listdir(GENE_CALLS_DIR) if file.endswith(".fasta") ]
    
    for file in tqdm(fasta_files):
        out_path = os.path.join(NT_SEQUENCES_DIR, file.replace(".fasta", ".nt"))
        if not os.path.exists(out_path):
            fasta_path = os.path.join(GENE_CALLS_DIR, file)
            lines = []
            
            for record in SeqIO.parse(fasta_path, "fasta"):
                parts = record.description.split(" # ")
                name = parts[0]
                partial = parts[4].split("partial=")[1].split(";")[0]
                direction = parts[3]
                partial = partial if direction == "1" else partial[::-1]
                lines.append(f"{name}|{partial}")
                lines.append(str(record.seq))
            
            with open(out_path, "w") as f:
                f.write("\n".join(lines))


def count_types():
    nt_files = [ file for file in os.listdir(NT_SEQUENCES_DIR) if file.endswith(".nt") ]
    counter = Counter()

    for file in tqdm(nt_files):
        with open(os.path.join(NT_SEQUENCES_DIR, file)) as f:
            sequence_headers = [ line for i, line in enumerate(f.read().splitlines()) if i % 2 == 0 ]
            partial_types = []
            partial_types = [ header.split("|")[1] for header in sequence_headers ]
            counter.update(partial_types)
    
    print(counter)


def count_collapsed_incomplete():
    nt_files = [ file for file in os.listdir(NT_SEQUENCES_DIR) if file.endswith(".nt") ]
    complete = set()

    for file in tqdm(nt_files):
        with open(os.path.join(NT_SEQUENCES_DIR, file)) as f:
            lines = f.read().splitlines()
            for i in range(len(lines)):
                if i % 2 == 0:
                    header = lines[i]
                    sequence = lines[i + 1]
                    partial_type = header.split("|")[1]
                    if partial_type == "00":
                        complete.add(sequence)

    print(f"Unique complete: {len(complete)}")


def main():
    #generate_nt_sequences()
    #count_types()
    count_collapsed_incomplete()

if __name__ == "__main__":
    main()