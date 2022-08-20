import sys
import glob
import Bio.SeqIO.FastaIO as FastaIO
from Bio.Seq import Seq
import pathlib
import subprocess
import collections
import pprint


AAHEADER = 'aaheader'
NTHEADER = 'ntheader'
NTUNPADDED = 'ntunpadded'
AAUNPADDED = 'aaunpadded'
NTPADDED = 'ntpadded'
NTPADDEDSTART = 'ntpaddedstart'
NTPADDEDSTOP = 'ntpaddedstop'
NTPADDEDLENGTH = 'ntpaddedlength'
NTUNPADDEDLENGTH = 'ntunpaddedlength'
NTSCAFFOLD = 'ntscaffold'

mOTUs_MGs = set(['COG0012',
            'COG0016',
            'COG0018',
            'COG0172',
            'COG0215',
            'COG0495',
            'COG0525',
            'COG0533',
            'COG0541',
            'COG0552'])

def get_sequences(seq_file):
    sequences = []
    with open(seq_file) as handle:
        for (header, sequence) in FastaIO.SimpleFastaParser(handle):
            sequences.append((header, sequence.upper()))
    return sequences

def revcomp(sequence):
    return (sequence[0] + '-rev', str(Seq(sequence[1]).reverse_complement()))





def main():
    #$sequence_file $new_database_folder/dbs/$genome_ID $fetchMGs_folder
    genome_file = sys.argv[1]
    output_folder = sys.argv[2]
    fetchMGs_folder = sys.argv[3]
    genome_name = sys.argv[4]
    vsearch_database = sys.argv[5]

    #load all MGs

    mg_files_fna = glob.glob(f'{fetchMGs_folder}/*fna')
    cog_2_data = {}
    for cog in mOTUs_MGs:
        cog_2_data[cog] = {}
    for mg_file_fna in mg_files_fna:
        mg_file_faa = mg_file_fna.replace('.fna', '.faa')
        if mg_file_fna.split('/')[-1].replace('.fna', '') in mOTUs_MGs:
            cogid = mg_file_fna.split('/')[-1].replace('.fna', '')
            fna_sequences = get_sequences(mg_file_fna)
            faa_sequences = get_sequences(mg_file_faa)
            if len(fna_sequences) != len(faa_sequences) or (len(fna_sequences) + len(faa_sequences)) > 2:
                print(f'Genome {genome_file} and associated fetchMGs folder have wrong number of called genes in file {mg_file_fna} or {mg_file_faa}')
                raise  Exception
            if len(fna_sequences) != 0:
                cog_2_data[cogid][NTHEADER] = fna_sequences[0][0]
                cog_2_data[cogid][NTUNPADDED] = fna_sequences[0][1]
                cog_2_data[cogid][AAHEADER] = faa_sequences[0][0]
                cog_2_data[cogid][AAUNPADDED] = faa_sequences[0][1]
                cog_2_data[cogid][NTUNPADDEDLENGTH] = len(cog_2_data[cogid][NTUNPADDED])


    #pprint.pprint(cog_2_data)


    # find the sequences in the genome and generate padded fna sequences

    genome_sequences = [sequence for sequence in get_sequences(genome_file)]
    rev_genome_sequences = [revcomp(sequence) for sequence in genome_sequences]
    genome_sequences = genome_sequences + rev_genome_sequences


    for cog, data in cog_2_data.items():
        if len(data) == 0:
            continue
        ntunpadded = data[NTUNPADDED]
        found = False
        for scaffold, genome_sequence in genome_sequences:
            if ntunpadded in genome_sequence:
                splits = genome_sequence.split(ntunpadded)
                prefix = splits[0]
                suffix = splits[1]

                if len(prefix)>100:
                    prefix = prefix[len(prefix)-100:]

                if len(suffix)>100:
                    suffix = suffix[:100]
                ntpadded = prefix + ntunpadded + suffix
                startpos = len(prefix) + 1
                endpos = len(prefix + ntunpadded)
                ntpaddedlength = len(ntpadded)
                cog_2_data[cog][NTPADDED] = ntpadded
                cog_2_data[cog][NTPADDEDLENGTH] = ntpaddedlength
                cog_2_data[cog][NTPADDEDSTART] = startpos
                cog_2_data[cog][NTPADDEDSTOP] = endpos
                cog_2_data[cog][NTSCAFFOLD] = scaffold
                found = True
                break
        if not found:
            print(f'Didnt find the cog {cog} sequence in genome {genome_file}')
            raise Exception


    pathlib.Path(output_folder).mkdir(exist_ok=True, parents=True)
    # Output files:
    # GENOME_NAME.padded.coords --> used in generateDB
    # GENOME_NAME.{COG}     SCAFFOLD    NTPADDEDSTART   NTPADDEDSTOP --> Only for existing cogs
    padded_coords_file = f'{output_folder}/{genome_name}.padded.coords'
    print(f'Writing {padded_coords_file}')
    with open(padded_coords_file, 'w') as outhandle:
        for cog, data in cog_2_data.items():
            if len(data) == 0:
                continue
            outhandle.write(f'{genome_name}.{cog}\t{data[NTSCAFFOLD]}\t{data[NTPADDEDSTART]}\t{data[NTPADDEDSTOP]}\n')



    # GENOME_NAME.genes.len --> used in addMarkerGenes
    # GENOME_NAME.{COG}         NTUNPADDEDLENGTH

    genes_len_file = f'{output_folder}/{genome_name}.genes.len'
    print(f'Writing {genes_len_file}')
    with open(genes_len_file, 'w') as outhandle:
        for cog, data in cog_2_data.items():
            if len(data) == 0:
                continue
            outhandle.write(f'{genome_name}.{cog}\t{data[NTUNPADDEDLENGTH]}\n')

    # GENOME_NAME.padded.fasta --> used in generateDB
    # GENOME_NAME.{COG}
    # NTPADDED

    padded_fasta_file = f'{output_folder}/{genome_name}.padded.fasta'
    print(f'Writing {padded_fasta_file}')
    with open(padded_fasta_file, 'w') as outhandle:
        for cog, data in cog_2_data.items():
            if len(data) == 0:
                continue
            outhandle.write(f'>{genome_name}.{cog}\n{data[NTPADDED]}\n')

    # GENOME_NAME.map2genome --> needed everywhere
    # GENOME_NAME.{COG}     GENOME_NAME
    map2genome_file = f'{output_folder}/{genome_name}.map2genome'
    print(f'Writing {map2genome_file}')
    with open(map2genome_file, 'w') as outhandle:
        for cog, data in cog_2_data.items():
            if len(data) == 0:
                continue
            outhandle.write(f'{genome_name}.{cog}\t{genome_name}\n')


    # GENOME_NAME_markerGenes2
    markergenes2_folder = f'{output_folder}/{genome_name}_markerGenes2/'
    pathlib.Path(markergenes2_folder).mkdir(exist_ok=True, parents=True)
    # GENOME_NAME_markerGenes2/all.notab.fna
    # GENOME_NAME.{COG}
    # NTUNPADDED
    all_notab_fasta_file = f'{markergenes2_folder}all.notab.fna'
    print(f'Writing {all_notab_fasta_file}')
    with open(all_notab_fasta_file, 'w') as outhandle:
        for cog, data in cog_2_data.items():
            if len(data) == 0:
                continue
            outhandle.write(f'>{genome_name}.{cog}\n{data[NTUNPADDED]}\n')

    # create the individual notab.fna files
    for cog, data in cog_2_data.items():
        with open(f'{markergenes2_folder}{cog}.notab.fna', 'w') as outhandle:
            if len(data) == 0:
                continue
            else:
                outhandle.write(f'>{genome_name}.{cog}\n{data[NTUNPADDED]}\n')
    # create the SAMPLE.map file
    mapfile = f'{output_folder}/{genome_name}.map'
    with open(mapfile, 'w') as outhandle:
        for cog, data in cog_2_data.items():
            if len(data) == 0:
                continue
            upl = data['ntunpaddedlength']
            outhandle.write(f'{genome_name}.{cog}\t{upl}\t{cog}\t{genome_name}\n')
    # Empty vsearch files
    vsearch_folder = f'{output_folder}/vsearch/'
    pathlib.Path(vsearch_folder).mkdir(exist_ok=True, parents=True)
    for cog, data in cog_2_data.items():
        empty_m8_file = f'{vsearch_folder}{cog}.distances_vs_db.m8'
        print(f'Create empty m8 file {empty_m8_file}')
        pathlib.Path(empty_m8_file).touch(exist_ok=True)

    all_notab_m8_file = f'{markergenes2_folder}all.notab.m8'
    command = f'vsearch --threads 1 --usearch_global {all_notab_fasta_file} --strand both --db {vsearch_database} --id 0.0  --maxaccepts 1000 --maxaccepts 1000 --mincols 20 --blast6out {all_notab_m8_file}'
    print(f'Creating alignment file: {all_notab_m8_file}')
    try:
        print(f'Executing command: {command}')
        returncode = subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print('Command {} failed with message:\t{}'.format(e.cmd, e.stderr))
        exit(1)

    cog_2_alignments = collections.defaultdict(list)
    with open(all_notab_m8_file) as handle:
        for line in handle:
            cog = line.strip().split()[0].rsplit('.', 1)[1]
            cog_2_alignments[cog].append(line.strip())
    for cog, alignments in cog_2_alignments.items():
        empty_m8_file = f'{vsearch_folder}{cog}.distances_vs_db.m8'
        print(f'Filling empty m8 file {empty_m8_file}')
        with open(empty_m8_file, 'w') as outhandle:
            for alignment in alignments:
                outhandle.write(f'{alignment}\n')












if __name__ == '__main__':
    main()
