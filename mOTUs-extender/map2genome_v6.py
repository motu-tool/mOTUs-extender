import sys

input_file = sys.argv[1]
output_file = sys.argv[2]
mgs = ['COG0172','COG0012','COG0215','COG0525','COG0495','COG0541','COG0533','COG0552','COG0016','COG0018']
with open(input_file) as inhandle, open(output_file, 'w') as outhandle:
    for line in inhandle:
        
        splits = line.strip().split()
        gene = splits[0]
        genome = splits[1]
        thismg = []
        for mg in mgs:
            if mg in gene:
                thismg.append(mg)
        if len(thismg) != 1:
            print('Either no MG or multiple MGs are found. SHould find exactly one. Quitting')
            exit(1)
        cog = thismg[0]
        #cog = 'COG' + gene.split('.COG')[1].split('.')[0].split('-')[0]
        #print(gene, genome, cog)
        outhandle.write(f'{gene}\t{genome}\t{cog}\n')
