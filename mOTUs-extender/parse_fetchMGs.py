import glob
import sys
import Bio.SeqIO.FastaIO as FastaIO

motu_mgs = set(['COG0012',
                   'COG0016',
                   'COG0018',
                   'COG0172',
                   'COG0215',
                   'COG0495',
                   'COG0525',
                   'COG0533',
                   'COG0541',
                   'COG0552'])


fmp = sys.argv[1]
samplename = sys.argv[3]
outfolder = sys.argv[2]
fna_files = glob.glob(fmp + '/*fna')
faa_files = glob.glob(fmp + '/*faa')

motus_cogs = []

for fna_file in fna_files:
    cog = fna_file.split('/')[-1].rsplit('.', 1)[0]
    with open(outfolder + '/'  + cog + '.fna', 'w') as outhandle:
        with open(fna_file) as inhandle:
            for (header, sequence) in FastaIO.SimpleFastaParser(inhandle):
                outhandle.write(f'>{samplename}.{cog}\n{sequence}\n')
                if cog in motu_mgs:
                    motus_cogs.append(cog)

with open(outfolder + '/motus.mgs.count', 'w') as outhandle:
    for cog in sorted(motus_cogs):
        outhandle.write(f'{cog}\n')
    
from pathlib import Path
if len(motus_cogs) > 5:
    Path(outfolder + '/motus.mgs.count.ok').touch()

for faa_file in faa_files:
    cog = faa_file.split('/')[-1].rsplit('.', 1)[0]
    with open(outfolder + '/'  + cog + '.faa', 'w') as outhandle:
        with open(faa_file) as inhandle:
            for (header, sequence) in FastaIO.SimpleFastaParser(inhandle):
                outhandle.write(f'>{samplename}.{cog}\n{sequence}\n')
