import os
from tqdm import tqdm
from collections import defaultdict
from constants import *

SUBSTRING_LEN = 10

def build_substring_dict(sequences):
    substring_dict = defaultdict(lambda: set())
    for i, sequence in tqdm(list(enumerate(sequences))):
        substrings = [sequence[i:i+SUBSTRING_LEN] for i in range(len(sequence) - SUBSTRING_LEN + 1)]
        for substring in substrings:
            substring_dict[substring].add(i)
    return substring_dict


def is_contained(sequence, substring_dict, sequences):
    substrings = [sequence[i:i+SUBSTRING_LEN] for i in range(0, len(sequence) - len(sequence) % SUBSTRING_LEN, SUBSTRING_LEN)]
    containing_sets = [ substring_dict[substring] for substring in substrings ]
    containing_indices = set.intersection(*containing_sets)

    for i in containing_indices:
        if sequence in sequences[i]:
            return True
    return False


def main():
    with open(UNIQUE_COMPLETE_SEQUENCES) as f:
        complete_sequences = f.read().splitlines()
    with open(UNIQUE_PARTIAL_SEQUENCES) as f:
        partial_sequences = f.read().splitlines()
    
    substring_dict = build_substring_dict(complete_sequences)
    contained_partial_sequences = [ sequence for sequence in tqdm(partial_sequences) if is_contained(sequence, substring_dict, complete_sequences) ]
    print(f"Partial sequences contained in complete sequences {len(contained_partial_sequences)}")


if __name__ == "__main__":
    main()