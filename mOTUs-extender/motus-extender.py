import argparse
import collections
import sys
import logging
import pathlib
import subprocess
import shutil



MOTUS_EXTENDER_VERSION = '3.0.1'
SCRIPT_DIR = pathlib.Path(__file__).parent.resolve()

########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################


class CapitalisedHelpFormatter(argparse.HelpFormatter):
    def add_usage(self, usage, actions, groups, prefix=None):
        if prefix is None:
            prefix = ''
        return super(CapitalisedHelpFormatter, self).add_usage(usage, actions, groups, prefix)



########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################


def startup():
    """Generic startup method
    """
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=logging.INFO)
    logging.info('Starting mOTUs-extender')



########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################

def shutdown(status=0):
    """Generic shutdown method
    """
    if status == 0:
        logging.info(f'mOTUs-extender finished successfully')
    else:
        logging.info(f'mOTUs-extender finished with error code {status}')
    sys.exit(status)



########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################

def check_call(command: str):
    """
    Simple wrapper to execute check_call and catch exceptions
    :param command:
    :return:
    """

    returncode = -1
    try:
        logging.info(f'Executing: {command}')
        returncode = subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        logging.error('Command {} failed with message:\t{}'.format(e.cmd, e.stderr))
        shutdown(1)


########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################

def execute_motus_prepare(args):
    parser = argparse.ArgumentParser(
        description='mOTUs-extender - extend the mOTUs database using your own genomes', usage=f'''
Program:    mOTUs-extender - extend the mOTUs database using your own genomes
Version:    {MOTUS_EXTENDER_VERSION}
Reference:  Milanese et al. Microbial abundance, activity and population 
            genomic profiling with mOTUs2; Nature Communications 10, 
            Article number: 1014 (2019). PMID: 30833550; 
            doi: 10.1038/s41467-019-08844-4
                
Usage: motus-extender prepare [options]

Input options:
    -w  DIR   Work folder for the extension
    ''', formatter_class=CapitalisedHelpFormatter, add_help=False)

    parser.add_argument('-w', action='store', dest='w', help='Work folder for the extension', required=True)

    if len(args) == 0:
        parser.print_usage()
        shutdown(1)

    results = parser.parse_args(args)

    #######
    #Input#
    #######
    input_workfolder = pathlib.Path(results.w).resolve()

    if input_workfolder.is_file():
        logging.error('Folder exists and is file. Select folder.')
        shutdown(1)
    if input_workfolder.is_dir():
        logging.error('Folder exists. Delete folder or select a non-existing folder.')
        shutdown(1)

    logging.info(f'Creating work folder {input_workfolder}')
    input_workfolder.mkdir(parents=True, exist_ok=True)
    logging.info(f'Downloading mOTUs {MOTUS_EXTENDER_VERSION} database')
    orig_db_folder = input_workfolder.joinpath('orig_db')
    orig_db_folder.mkdir(parents=True, exist_ok=True)
    download_command = f'curl https://zenodo.org/record/5140350/files/db_mOTU_v3.0.1.tar.gz -o {orig_db_folder}/db_mOTU_v3.0.1.tar.gz'
    check_call(download_command)
    logging.info(f'Uncompressing mOTUs {MOTUS_EXTENDER_VERSION} database')
    untar_command = f'tar -xzvf {orig_db_folder}/db_mOTU_v3.0.1.tar.gz -C {orig_db_folder}'
    check_call(untar_command)
    logging.info(f'Preparing mOTUs {MOTUS_EXTENDER_VERSION} database for extension')
    prepare_db_command = f'python {SCRIPT_DIR}/prepare_mOTUs.py {orig_db_folder}/db_mOTU/ {input_workfolder}/temp_db_folder'
    check_call(prepare_db_command)
    marker_file = input_workfolder.joinpath('motus-extender-prepare.done')
    marker_file.touch()

########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################


def execute_motus_membership(args):
    parser = argparse.ArgumentParser(
        description='mOTUs-extender - extend the mOTUs database using your own genomes', usage=f'''
Program:    mOTUs-extender - extend the mOTUs database using your own genomes
Version:    {MOTUS_EXTENDER_VERSION}
Reference:  Milanese et al. Microbial abundance, activity and population 
            genomic profiling with mOTUs2; Nature Communications 10, 
            Article number: 1014 (2019). PMID: 30833550; 
            doi: 10.1038/s41467-019-08844-4

Usage: motus-extender membership [options]

Input options:
    -w  DIR   Work folder for the extension
    -g  DIR   Folder with genomes. File extension: .fa
Other options:
    -t  INT               Number of threads. [4]
    ''', formatter_class=CapitalisedHelpFormatter, add_help=False)

    parser.add_argument('-w', action='store', dest='w', help='Work folder for the extension', required=True)
    parser.add_argument('-g', action='store', dest='g', help='Folder with genomes. File extension: .fa', required=True)
    parser.add_argument('-t', action='store', dest='threads', help='Number of threads.',default=4, type=int)
    parser.add_argument('-n', action='store_true', dest='n', help='dryrun')

    if len(args) == 0:
        parser.print_usage()
        shutdown(1)

    results = parser.parse_args(args)

    #######
    #Input#
    #######
    input_workfolder = pathlib.Path(results.w).resolve()
    genomes_folder = pathlib.Path(results.g).resolve()
    threads = results.threads
    dryrun = results.n


    # check correctness of workfolder

    if not input_workfolder.is_dir():
        logging.error('Work folder does not exist.')
        shutdown(1)
    prepare_marker_file = input_workfolder.joinpath('motus-extender-prepare.done')
    membership_marker_file = input_workfolder.joinpath('motus-extender-membership.done')
    if not prepare_marker_file.exists():
        logging.error('motus-extender prepare command wasn\'t executed or didn\'t finish. Rerun the motus-extender prepare command')
        shutdown(1)
    if membership_marker_file.exists():
        logging.error('motus-extender membership finished successfully before. Delete the marker file if new genomes were added.')
        shutdown(0)

    # copy mOTU.v3.0.mOTU-LG.map.gz to workfolder
    motu_map_file = f'{input_workfolder}/mOTU.v3.0.mOTU-LG.map'
    unzip_command = f'gunzip -c {SCRIPT_DIR}/mOTU.v3.0.mOTU-LG.map.gz > {motu_map_file}'
    check_call(unzip_command)
    # check the number of genomes

    genomes = sorted(list(genomes_folder.glob('*fa')))
    logging.info(f'Found {len(genomes)} genomes in {genomes_folder} path')
    chunk_size = 300


    genome_batches = [genomes[i * chunk_size:(i + 1) * chunk_size] for i in range((len(genomes) + chunk_size - 1) // chunk_size )]



    temp_db_folder = input_workfolder.joinpath('temp_db_folder')
    for genome_batch in genome_batches:
        tmp_genomes = ','.join([str(genome) for genome in genome_batch])
        sm_command = f'snakemake -s {SCRIPT_DIR}/0prod_fetchMGs.py --config mapfile={motu_map_file} infolder={genomes_folder}/ outfolder={input_workfolder}/ temp_db_folder={temp_db_folder}/ scriptfolder={SCRIPT_DIR} genomes={tmp_genomes} -j {threads} -k'
        if dryrun:
            sm_command = sm_command + ' -n -q'
        snakemake_folder = pathlib.Path.cwd().joinpath('.snakemake')
        if snakemake_folder.exists:
            shutil.rmtree(snakemake_folder, ignore_errors=True)
        check_call(sm_command)
        if snakemake_folder.exists:
            shutil.rmtree(snakemake_folder, ignore_errors=True)

    genome_2_assignation = {}
    if not dryrun:
        logging.info('Processing genome assignation')
        for cnt, genome in enumerate(genomes, 1):

            if cnt % 100 == 0:
                logging.info(f'Processed {cnt} / {len(genomes)}')
            genome_name = str(genome.name).rsplit('.fa')[0]

            combined_distance_file = input_workfolder.joinpath(f'extension/dbs/{genome_name}/vsearch/combined.normalized.distances_vs_db.m8')
            mg_count_file = input_workfolder.joinpath(f'genes/{genome_name}/fetchMGs/{genome_name}-bestMGs-renamed/motus.mgs.count')

            if not mg_count_file.exists():
                logging.error('Not all genomes in the list were processed. E.g. {genome_name} is missing. Quitting')
                shutdown(1)

            mg_count = len(set(mg_count_file.read_text().splitlines()))
            if not combined_distance_file.exists():
                genome_2_assignation[genome_name]= ('NotEnoughMGs', mg_count)
            else:
                cutoff = 99.0
                pre_assigned_motu_distance_2_pre_assigned_motu = collections.defaultdict(list)
                for line in combined_distance_file.read_text().splitlines():
                    splits = line.strip().split('\t')
                    pre_assigned_motu = splits[1]
                    pre_assigned_motu_distance = float(splits[2])
                    if pre_assigned_motu_distance < cutoff:
                        continue
                    pre_assigned_motu_distance_2_pre_assigned_motu[pre_assigned_motu_distance].append(pre_assigned_motu)
                assigned_motu = 'Novel'
                if len(pre_assigned_motu_distance_2_pre_assigned_motu) != 0:
                    assigned_motu = ';'.join(pre_assigned_motu_distance_2_pre_assigned_motu[max(pre_assigned_motu_distance_2_pre_assigned_motu.keys())])
                genome_2_assignation[genome_name]= (assigned_motu, mg_count)
        logging.info(f'Processed {cnt} / {len(genomes)}')

        membership_file = input_workfolder.joinpath('mOTUs.membership.tsv')
        novel = 0
        refmotus = 0
        metamotus = 0
        extmotus = 0
        excluded = 0
        with open(membership_file, 'w') as handle:
            for genome_name, (assigned_motu, mg_count) in genome_2_assignation.items():
                if 'ref' in assigned_motu:
                    refmotus += 1
                elif 'meta' in assigned_motu:
                    metamotus += 1
                elif 'ext' in assigned_motu:
                    extmotus += 1
                elif 'Novel' in assigned_motu:
                    novel += 1
                else:
                    excluded += 1

                handle.write(f'{genome_name}\t{assigned_motu}\t{mg_count}\n')
        logging.info(f'Total genomes: {len(genomes)}')
        logging.info(f'\t Assigned to refmotus \t\t\t{refmotus}')
        logging.info(f'\t Assigned to metamotus \t\t\t{metamotus}')
        logging.info(f'\t Assigned to extmotus \t\t\t{extmotus}')
        logging.info(f'\t Genomes not assigned to any mOTU \t{novel}')
        logging.info(f'\t Filtered genomes (<6 MGs) \t\t{excluded}')
        membership_marker_file.touch()



########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################

def execute_motus_createdb(args):
    parser = argparse.ArgumentParser(
        description='mOTUs-extender - extend the mOTUs database using your own genomes', usage=f'''
Program:    mOTUs-extender - extend the mOTUs database using your own genomes
Version:    {MOTUS_EXTENDER_VERSION}
Reference:  Milanese et al. Microbial abundance, activity and population 
            genomic profiling with mOTUs2; Nature Communications 10, 
            Article number: 1014 (2019). PMID: 30833550; 
            doi: 10.1038/s41467-019-08844-4

Usage: motus-extender createdb [options]

Input options:
    -w  DIR   Work folder for the extension
    -a  FILE  7 rank taxonomy file with NCBI naming. 
              One line per genome
    -p STR    Database prefix, 1-5 letters 
Other options:
    -t  INT               Number of threads. [4]
    ''', formatter_class=CapitalisedHelpFormatter, add_help=False)

    parser.add_argument('-w', action='store', dest='w', help='Work folder for the extension', required=True)
    parser.add_argument('-t', action='store', dest='threads', help='Number of threads.', default=4, type=int)
    parser.add_argument('-a', action='store', dest='a', help='7 rank taxonomy file with NCBI naming.', required=True)
    parser.add_argument('-p', action='store', dest='p', help='database extension prefix, 1-5 letters', required=True)

    if len(args) == 0:
        parser.print_usage()
        shutdown(1)

    results = parser.parse_args(args)

    #######
    #Input#
    #######
    input_workfolder = pathlib.Path(results.w).resolve()
    threads = results.threads
    taxonomy_file = pathlib.Path(results.a).resolve()
    db_prefix = results.p

    # check correctness of workfolder

    if not input_workfolder.is_dir():
        logging.error('Work folder does not exist.')
        shutdown(1)
    prepare_marker_file = input_workfolder.joinpath('motus-extender-prepare.done')
    membership_marker_file = input_workfolder.joinpath('motus-extender-membership.done')
    createdb_marker_file = input_workfolder.joinpath('motus-extender-createdb.done')
    if not prepare_marker_file.exists():
        logging.error('motus-extender prepare command wasn\'t executed or didn\'t finish. Rerun the motus-extender prepare command')
        shutdown(1)
    if not membership_marker_file.exists():
        logging.error('motus-extender membership command wasn\'t executed or didn\'t finish. Rerun the motus-extender prepare command')
        shutdown(1)
    if createdb_marker_file.exists():
        logging.error('motus-extender createdb command already successfully finished before. Quitting')
        shutdown(0)

    if not taxonomy_file.exists():
        logging.error(f'Taxonomy file {taxonomy_file} does not exist. Please provide a valid file')
        shutdown(1)

    # copy mOTU.v3.0.mOTU-LG.map.gz to workfolder
    motu_map_file = f'{input_workfolder}/mOTU.v3.0.mOTU-LG.map'
    membership_file = input_workfolder.joinpath('mOTUs.membership.tsv')
    novel_genomes  = [line.strip().split('\t')[0] for line in membership_file.read_text().splitlines() if 'Novel' in line]
    all_taxonomy_entries = {line.strip().split('\t')[0]: line.strip().split('\t')[1:]  for line in taxonomy_file.read_text().splitlines() if not line.startswith('#')}
    genome_list_file = input_workfolder.joinpath('genomes.list')
    with open(genome_list_file, 'w') as handle:
        for novel_genome in novel_genomes:
            novel_genome_tax = all_taxonomy_entries.get(novel_genome, f'No taxonomy found for {novel_genome}')
            if len(novel_genome_tax) != 7:
                logging.error(f'Taxonomy line for {novel_genome} malformed. Should be 7 column rank taxonomy but looks like: {novel_genome_tax}')
                shutdown(1)
            handle.write(f'{novel_genome}\n')


    create_db_command = f'sh {SCRIPT_DIR}/extend_mOTUs_generateDB.sh {genome_list_file} {db_prefix} {taxonomy_file} {input_workfolder}/extension/ {SCRIPT_DIR}/ {input_workfolder}/temp_db_folder/ {threads}'
    check_call(create_db_command)

    versions_command = f"sed -i 's/manual_update/{MOTUS_EXTENDER_VERSION}_{db_prefix}/g'{input_workfolder}/extension/{db_prefix}/db_mOTU/db_mOTU_versions"
    check_call(versions_command)



    # create a mapping file for all genomes and their motus associations
    createdb_marker_file.touch()








########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################

def main():

    parser = argparse.ArgumentParser(
        description='mOTUs-extender - extend the mOTUs database using your own genomes', usage=f'''
Program:    mOTUs-extender - extend the mOTUs database using your own genomes  
Version:    {MOTUS_EXTENDER_VERSION}
Reference:  Milanese et al. Microbial abundance, activity and population 
	        genomic profiling with mOTUs2; Nature Communications 10, 
            Article number: 1014 (2019). PMID: 30833550; 
            doi: 10.1038/s41467-019-08844-4

Usage: motus-extender <command> [options]
Command:
-- General
    prepare      Prepare the extender
    membership   Identify novel genomes
    createdb     Add novel genomes to database
        
    ''', formatter_class=CapitalisedHelpFormatter, add_help=False)

    parser.add_argument('command', help='Subcommand to run: prepare|membership|createdb', choices=['prepare', 'membership', 'createdb'])

    args = parser.parse_args(sys.argv[1:2])
    startup()


    command = args.command
    if command == 'prepare':
        execute_motus_prepare(sys.argv[2:])
        shutdown(0)
    elif command == 'membership':
        execute_motus_membership(sys.argv[2:])
        shutdown(0)
    elif command == 'createdb':
        execute_motus_createdb(sys.argv[2:])
        shutdown(0)
    else:
        parser.print_usage()
        shutdown(1)
    shutdown(0)


if __name__ == '__main__':
    main()
