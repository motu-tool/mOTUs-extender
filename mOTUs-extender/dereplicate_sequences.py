# Taken from Lucas/Hans specI update


import Bio.SeqIO.FastaIO as FastaIO
import sys, collections

infasta = sys.argv[1]
outfasta = open(sys.argv[2], 'w')
outcluster = open(sys.argv[3], 'w')

sequence_2_header = collections.defaultdict(list)
seq2count = collections.Counter()

with open(infasta) as handle:
    for (header, sequence) in FastaIO.SimpleFastaParser(handle):
        header = header.split('\t')[0]
        sequence_2_header[sequence].append(header)
        seq2count[sequence] += 1

for cnt, (sequence, ab) in enumerate(seq2count.most_common(), 1):
    headers = sequence_2_header[sequence]
    outfasta.write('>cluster_{};length={};size={}\n{}\n'.format(cnt, len(sequence), len(headers), sequence))
    outcluster.write('>cluster_{};length={};size={}\n'.format(cnt, len(sequence), len(headers)))
    for header in headers:
        outcluster.write('\t{}\n'.format(header))

outfasta.close()
outcluster.close()