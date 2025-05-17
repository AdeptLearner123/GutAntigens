import os
from tqdm import tqdm
from collections import Counter
from constants import *


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    unique_complete_sequences = set()
    unique_partial_sequences = set()
    total_sequences = 0
    frequencies = Counter()
    files = os.listdir(AA_SEQUENCES_DIR)
    
    for file in tqdm(files):
        path = os.path.join(AA_SEQUENCES_DIR, file)
        with open(path) as f:
            lines = f.read().splitlines()
            header_sequences = list(zip(lines[::2], lines[1::2]))
            partial_type_sequences = [ (header.split("|")[1], sequence) for header, sequence in header_sequences ]
            complete_sequences = [ sequence for partial_type, sequence in partial_type_sequences if partial_type == "00" ]
            partial_sequences = [ sequence for partial_type, sequence in partial_type_sequences if partial_type != "00" ]
            
            unique_complete_sequences.update(complete_sequences)
            unique_partial_sequences.update(partial_sequences)
            total_sequences += len(header_sequences)
            frequencies.update(complete_sequences)

    print(f"Total sequences: {total_sequences}")
    print(f"Unique complete sequences: {len(unique_complete_sequences)}")
    print(f"Unque partial sequences {len(unique_partial_sequences)}")

    with open(UNIQUE_COMPLETE_SEQUENCES, "w") as f:
        f.write("\n".join(unique_complete_sequences))
    
    with open(UNIQUE_PARTIAL_SEQUENCES, "w") as f:
        f.write("\n".join(unique_partial_sequences))

    with open(SEQUENCE_FREQUENCIES, "w") as f:
        lines = [ f"{seq}:{freq}" for seq, freq in frequencies.items() ]
        f.write("\n".join(lines))


if __name__ == "__main__":
    main()