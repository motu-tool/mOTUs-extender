import sys
import Bio.SeqIO.FastaIO as FastaIO

input_f = sys.argv[1]
output_f = sys.argv[2]

with open(input_f) as inhandle, open(output_f, 'w') as outhandle:
    for (header, sequence) in FastaIO.SimpleFastaParser(inhandle):
        header = header.split()[0]
        outhandle.write(f'>{header}\n{sequence}\n')
