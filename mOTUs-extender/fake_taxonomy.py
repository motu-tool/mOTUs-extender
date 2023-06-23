import sys


novel = set()
with open(sys.argv[1]) as handle:
    for line in handle:
        if 'Novel' in line:
            genome = line.strip().split('\t')[0]
            novel.add(genome)
#with open('test/genomes.tax') as handle:
#    for line in handle:
#        line = line.strip().split('\t')
#        print(line)

for n in sorted(list(novel)):
    tmp = [n] + ['0 NA'] * 7
    tmp = '\t'.join(tmp)
    print(tmp)
