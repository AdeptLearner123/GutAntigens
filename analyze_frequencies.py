from collections import Counter
from constants import *

def main():
    with open(SEQUENCE_FREQUENCIES) as f:
        lines = f.read().splitlines()
    frequencies = [ line.split(":")[1] for line in lines ]
    counter = Counter(frequencies)
    sorted_counter = sorted(counter.items(), key=lambda item: item[1], reverse=True)
    print(list(sorted_counter))


if __name__ == "__main__":
    main()