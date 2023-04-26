import sys

with open(sys.argv[1]) as handle:
    for line in handle:
        splits = line.strip().split()
        print(f'{splits[1]}\t{splits[0]}')
