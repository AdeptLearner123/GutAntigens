import os
import tempfile
from constants import *
from collections import defaultdict
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def translate_genomes():
    os.makedirs(AA_SEQUENCES_DIR, exist_ok=True)
    fasta_files = [ file for file in os.listdir(GENE_CALLS_DIR) if file.endswith(".fasta") ]
    for file in tqdm(fasta_files):
        translate(file)


def translate(file_name):
    input_path = os.path.join(GENE_CALLS_DIR, file_name)
    output_path = os.path.join(AA_SEQUENCES_DIR, file_name)

    # IDEA: to simplify this, remove os.path.exists and tempfile writing. This keeps the code clean, and if the user wants to optimize this, they can introduce the logic again later.
    if not os.path.exists(output_path):
        with tempfile.NamedTemporaryFile(
            mode='w',
            dir=AA_SEQUENCES_DIR,
            prefix=os.path.basename(output_path) + '.tmp.',
            delete=False
        ) as tmp_handle:
            tmp_filename = tmp_handle.name

            for record in SeqIO.parse(input_path, "fasta"):
                print("Processing " + record.id)
                 # Translate nucleotide to protein (assumes standard table)
                protein_seq = record.seq.translate(stop_symbol="")
                
                # Create a new SeqRecord for the protein sequence
                protein_record = SeqRecord(
                    protein_seq,
                    id=record.id,
                    description=record.description
                )
                
                # Write to output file
                SeqIO.write(protein_record, tmp_handle, "fasta")
            
            tmp_handle.close()
            os.replace(tmp_filename, output_path)


def collect_sequences():
    complete_sequences = defaultdict(lambda: [])
    partial_sequences = defaultdict(lambda: [])
    os.makedirs(RESULTS_DIR, exist_ok=True)

    for file in tqdm(os.listdir(AA_SEQUENCES_DIR)):
        sample_name = file.replace(".fasta", "")
        for record in SeqIO.parse(os.path.join(AA_SEQUENCES_DIR, file), "fasta"):
            sequence = str(record.seq)
            is_complete = record.description.split("partial=")[1].split(";")[0] == "00"
            if is_complete:
                complete_sequences[sequence].append(sample_name)
            else:
                partial_sequences[sequence].append(sample_name)
    
    with open(COMPLETE_SEQUENCES, "w") as output_handle:
        for i, sequence in tqdm(list(enumerate(complete_sequences))):
            samples_description = ";".join(complete_sequences[sequence])
            record = SeqRecord(
                Seq(sequence),
                id=f"{i}",
                description=f"{i} # {samples_description}"
            )
            SeqIO.write(record, output_handle, "fasta")

    with open(PARTIAL_SEQUENCES, "w") as output_handle:
        for i, sequence in tqdm(list(enumerate(partial_sequences))):
            samples_description = ";".join(partial_sequences[sequence])
            record = SeqRecord(
                Seq(sequence),
                id=f"{i}",
                description=f"{i} # {samples_description}"
            )
            SeqIO.write(record, output_handle, "fasta")


def main():
    translate_genomes()
    collect_sequences()


if __name__ == "__main__":
    main()