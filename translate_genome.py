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
    sequence_samples = defaultdict(lambda: [])
    sequence_completions = dict()
    os.makedirs(RESULTS_DIR, exist_ok=True)

    for file in tqdm(os.listdir(AA_SEQUENCES_DIR)):
        sample_name = file.replace(".fasta", "")
        for record in SeqIO.parse(os.path.join(AA_SEQUENCES_DIR, file), "fasta"):
            sequence = str(record.seq)
            is_complete = record.description.split("partial=")[1].split(";")[0] == "00"
            sequence_samples[sequence].append(sample_name)

            completion_str = "Complete" if is_complete else "Incomplete"
            if sequence not in sequence_completions:
                sequence_completions[sequence] = completion_str
            elif sequence_completions[sequence] != completion_str:
                sequence_completions[sequence] = "Mixed"
    
    with open(SEQUENCES, "w") as output_handle:
        for i, sequence in tqdm(list(enumerate(sequence_samples))):
            samples_description = ";".join(sequence_samples[sequence])
            record = SeqRecord(
                Seq(sequence),
                id=f"{i}",
                description=f"{i} # {sequence_completions[sequence]}|{samples_description}"
            )
            SeqIO.write(record, output_handle, "fasta")


def main():
    translate_genomes()
    collect_sequences()


if __name__ == "__main__":
    main()