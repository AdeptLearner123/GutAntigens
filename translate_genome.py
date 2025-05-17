import os
from constants import *
from tqdm import tqdm
from Bio import SeqIO

def translate(file_name):
    fasta_path = os.path.join(GENE_CALLS_DIR, file_name)
    output_file_name = file_name.replace(".fasta", ".aa")
    output_path = os.path.join(AA_SEQUENCES_DIR, output_file_name)

    if not os.path.exists(output_path):
        lines = []
        for record in SeqIO.parse(fasta_path, "fasta"):
                parts = record.description.split(" # ")
                name = parts[0]
                partial = parts[4].split("partial=")[1].split(";")[0]
                direction = parts[3]
                partial = partial if direction == "1" else partial[::-1]
                lines.append(f"{name}|{partial}")
                lines.append(str(record.translate(stop_symbol="").seq))

        with open(output_path, "w") as f:
            print("Writing to " + output_path)
            f.write("\n".join(lines))


def main():
    os.makedirs(AA_SEQUENCES_DIR, exist_ok=True)
    fasta_files = [ file for file in os.listdir(GENE_CALLS_DIR) if file.endswith(".fasta") ]
    for file in tqdm(fasta_files):
        translate(file)


if __name__ == "__main__":
    main()