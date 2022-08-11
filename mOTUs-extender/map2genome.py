import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file) as inhandle, open(output_file, 'w') as outhandle:
    for line in inhandle:
        
        splits = line.strip().split()
        gene = splits[0]
        genome = splits[1]
        cog = gene.split('.')[-1]
        outhandle.write(f'{gene}\t{genome}\t{cog}\n')
