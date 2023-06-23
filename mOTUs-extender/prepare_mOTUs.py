
import argparse
import sys
import Bio.SeqIO.FastaIO as FastaIO
import pathlib
import subprocess
import shutil
import os


def main(argv=None):
    
    parser = argparse.ArgumentParser(description='Prepare the mOTUs database for the extension pipeline - Removes marker genes that are part of the unassigned mOTU.', add_help = True)
    parser.add_argument('input', action="store", help='Input: The orginal db_mOTU folder.')
    parser.add_argument('output', action="store", help='Output: The temporary mOTU database folder')
    args = parser.parse_args()

    inputfolder = args.input
    outputfolder = args.output

    nomin1_output_folder = outputfolder + '/db_mOTU/'
    pathlib.Path(nomin1_output_folder).mkdir(parents=True, exist_ok=True)
    min1_output_folder = outputfolder + '/min1_db_mOTU/'
    pathlib.Path(min1_output_folder).mkdir(parents=True, exist_ok=True)


    '''
    1. check of the input folder format is correct
    '''
    versions_file = inputfolder + '/db_mOTU_versions'

    mgc_2_motus_file = inputfolder + '/db_mOTU_MAP_MGCs_to_mOTUs.tsv'
    mgc_2_gene_file = inputfolder + '/db_mOTU_MAP_genes_to_MGCs.tsv'
    db_mOTU_MAP_MGCs_to_mOTUs_inline_file = inputfolder + '/db_mOTU_MAP_MGCs_to_mOTUs_in-line.tsv'
    db_mOTU_padding_coordinates_NR_file = inputfolder + '/db_mOTU_padding_coordinates_NR.tsv'
    db_mOTU_taxonomy_refmOTUs_short_names_file = inputfolder + '/db_mOTU_taxonomy_ref-mOTUs_short_names.tsv'
    db_mOTU_genes_length_NR_file = inputfolder + '/db_mOTU_genes_length_NR'
    db_mOTU_taxonomy_CAMI_file = inputfolder + '/db_mOTU_taxonomy_CAMI.tsv'
    db_mOTU_taxonomy_refmOTUs_file = inputfolder + '/db_mOTU_taxonomy_ref-mOTUs.tsv'
    README_file = inputfolder + '/README'
    db_mOTU_taxonomy_metamOTUs_file = inputfolder + '/db_mOTU_taxonomy_meta-mOTUs.tsv'
    db_mOTU_test_folder = inputfolder + '/db_mOTU_test'
    db_mOTU_padding_coordinates_CEN_file = inputfolder + '/db_mOTU_padding_coordinates_CEN.tsv'
    db_mOTU_DB_CEN_fastaannotations_file = inputfolder + '/db_mOTU_DB_CEN.fasta.annotations'
    db_mOTU_bam_header_NR_file = inputfolder + '/db_mOTU_bam_header_NR'
    db_mOTU_bam_header_CEN_file = inputfolder + '/db_mOTU_bam_header_CEN'
    nr_fasta_file = inputfolder + '/db_mOTU_DB_NR.fasta'
    db_mOTU_DB_CEN_file = inputfolder + '/db_mOTU_DB_CEN.fasta'









    print('####################################')
    print('#### Checking db_mOTU_versions #####')
    print('####################################')
    with open(versions_file) as handle:
        for line in handle:
            splits = line.strip().split()
            if splits[0] != '#':
                if splits[1] == '3.1.0':
                    print(f'Version of {splits[0]} is {splits[1]}. Correct version. Continue')
                else:
                    print(f'Version of {splits[0]} is {splits[1]}. WRONG VERSION. Stopping')
                    sys.exit(1)













    print('###################################################')
    print('##### Processing db_mOTU_MAP_MGCs_to_mOTUs.tsv ####')
    print('###################################################')

    number_of_mgcs = 0
    min1_mgcs = set()
    min1_mgcs_2_length = {}

    min1_db_mOTU_MAP_MGCs_to_mOTUs_file = min1_output_folder + '/db_mOTU_MAP_MGCs_to_mOTUs.tsv'
    nomin1_db_mOTU_MAP_MGCs_to_mOTUs_file = nomin1_output_folder + '/db_mOTU_MAP_MGCs_to_mOTUs.tsv'

    with open(mgc_2_motus_file) as inhandle, open(min1_db_mOTU_MAP_MGCs_to_mOTUs_file, 'w') as min1_handle, open(nomin1_db_mOTU_MAP_MGCs_to_mOTUs_file, 'w') as nomin1_handle:

        for line in inhandle:
            number_of_mgcs += 1
            splits = line.strip().split()
            if splits[1] == 'unassigned':
                min1_handle.write(f'{line.strip()}\n')
                min1_mgcs.add(splits[0])
            else:
                nomin1_handle.write(f'{line.strip()}\n')

    print(f'Found {number_of_mgcs} MGCs of which {len(min1_mgcs)} MGCs belong to unassigned')




















    print('###################################################')
    print('#### Processing db_mOTU_MAP_genes_to_MGCs.tsv #####')
    print('###################################################')

    number_of_marker_genes = 0
    min1_genes = set()
    min1_db_mOTU_MAP_genes_to_MGC_file = min1_output_folder + '/db_mOTU_MAP_genes_to_MGCs.tsv'
    nomin1_db_mOTU_MAP_genes_to_MGC_file = nomin1_output_folder + '/db_mOTU_MAP_genes_to_MGCs.tsv'
    with open(mgc_2_gene_file) as inhandle, open(min1_db_mOTU_MAP_genes_to_MGC_file, 'w') as min1_handle, open(nomin1_db_mOTU_MAP_genes_to_MGC_file, 'w') as nomin1_handle:
        for line in inhandle:
            number_of_marker_genes += 1

            splits = line.strip().split()
            marker_gene = splits[0]
            mgc = splits[3]
            if mgc in min1_mgcs:
                min1_genes.add(marker_gene)
                min1_handle.write(f'{line.strip()}\n')
            else:
                nomin1_handle.write(f'{line.strip()}\n')

    print(f'Found {number_of_marker_genes} MGs of which {len(min1_genes)} MGs belong to unassigned')














    print('####################################################')
    print('# Processing db_mOTU_MAP_MGCs_to_mOTUs_in-line.tsv #')
    print('####################################################')


    min1_db_mOTU_MAP_MGCs_to_mOTUs_inline_file = min1_output_folder + '/db_mOTU_MAP_MGCs_to_mOTUs_in-line.tsv'
    nomin1_db_mOTU_MAP_MGCs_to_mOTUs_inline_file = nomin1_output_folder + '/db_mOTU_MAP_MGCs_to_mOTUs_in-line.tsv'
    min = 0
    nomin = 0
    with open(db_mOTU_MAP_MGCs_to_mOTUs_inline_file) as inhandle, open(min1_db_mOTU_MAP_MGCs_to_mOTUs_inline_file, 'w') as min1_handle, open(nomin1_db_mOTU_MAP_MGCs_to_mOTUs_inline_file, 'w') as nomin1_handle:
        for line in inhandle:
            if line.strip().startswith('unassigned'):
                min1_handle.write(f'{line.strip()}\n')
                min += 1
            else:
                nomin1_handle.write(f'{line.strip()}\n')
                nomin += 1

    print(f'unassigned mOTU written: {min}')
    print(f'ref-mOTU/meta-mOTU mOTUs written: {nomin}')













    print('####################################################')
    print('## Processing db_mOTU_padding_coordinates_NR.tsv ###')
    print('####################################################')

    min1_db_mOTU_padding_coordinates_NR_file = min1_output_folder + '/db_mOTU_padding_coordinates_NR.tsv'
    nomin1_db_mOTU_padding_coordinates_NR_file = nomin1_output_folder + '/db_mOTU_padding_coordinates_NR.tsv'

    min = 0
    nomin = 0
    with open(db_mOTU_padding_coordinates_NR_file) as inhandle, open(min1_db_mOTU_padding_coordinates_NR_file, 'w') as min1_handle, open(nomin1_db_mOTU_padding_coordinates_NR_file, 'w') as nomin1_handle:
        for line in inhandle:
            splits = line.strip().split()
            if splits[0] in min1_genes:
                min += 1
                min1_mgcs_2_length[splits[0]] = ( int(splits[2]), int(splits[3]) )

                min1_handle.write(f'{line.strip()}\n')
            else:
                nomin += 1
                nomin1_handle.write(f'{line.strip()}\n')

    print(f'-1 MGs written: {min}')
    print(f'ref-mOTU/meta-mOTU MGs written: {nomin}')





    print('####################################################')
    print('## Processing db_mOTU_padding_coordinates_CEN.tsv ##')
    print('####################################################')

    min1_db_mOTU_padding_coordinates_CEN_file = min1_output_folder + '/db_mOTU_padding_coordinates_CEN.tsv'
    nomin1_db_mOTU_padding_coordinates_CEN_file = nomin1_output_folder + '/db_mOTU_padding_coordinates_CEN.tsv'

    min = 0
    nomin = 0
    with open(db_mOTU_padding_coordinates_CEN_file) as inhandle, open(min1_db_mOTU_padding_coordinates_CEN_file, 'w') as min1_handle, open(nomin1_db_mOTU_padding_coordinates_CEN_file, 'w') as nomin1_handle:
        for line in inhandle:
            splits = line.strip().split()
            if splits[0].startswith('unassigned.'):
                min += 1
                min1_handle.write(f'{line.strip()}\n')
            else:
                nomin += 1
                nomin1_handle.write(f'{line.strip()}\n')

    print(f'unassigned MGs written: {min}')
    print(f'ref-mOTU/meta-mOTU MGs written: {nomin}')












    print('#####################################################')
    print('Copying db_mOTU_taxonomy_ref-mOTUs_short_names.tsv')
    print('#####################################################')

    nomin_db_mOTU_taxonomy_refmOTUs_short_names_file = nomin1_output_folder + '/db_mOTU_taxonomy_ref-mOTUs_short_names.tsv'
    shutil.copyfile(db_mOTU_taxonomy_refmOTUs_short_names_file, nomin_db_mOTU_taxonomy_refmOTUs_short_names_file)

    print('#####################################################')
    print('Copying db_mOTU_taxonomy_CAMI.tsv')
    print('#####################################################')

    nomin_db_mOTU_taxonomy_CAMI_file = nomin1_output_folder + '/db_mOTU_taxonomy_CAMI.tsv'
    shutil.copyfile(db_mOTU_taxonomy_CAMI_file, nomin_db_mOTU_taxonomy_CAMI_file)

    print('#####################################################')
    print('Copying db_mOTU_taxonomy_ref-mOTUs.tsv')
    print('#####################################################')

    nomin_db_mOTU_taxonomy_refmOTUs_file = nomin1_output_folder + '/db_mOTU_taxonomy_ref-mOTUs.tsv'
    shutil.copyfile(db_mOTU_taxonomy_refmOTUs_file, nomin_db_mOTU_taxonomy_refmOTUs_file)


    print('#####################################################')
    print('Copying README')
    print('#####################################################')

    nomin_README_file = nomin1_output_folder + '/README'
    shutil.copyfile(README_file, nomin_README_file)


    print('#####################################################')
    print('Copying db_mOTU_taxonomy_meta-mOTUs.tsv')
    print('#####################################################')

    nomin_db_mOTU_taxonomy_metamOTUs_file = nomin1_output_folder + '/db_mOTU_taxonomy_meta-mOTUs.tsv'
    shutil.copyfile(db_mOTU_taxonomy_metamOTUs_file, nomin_db_mOTU_taxonomy_metamOTUs_file)


    print('#####################################################')
    print('Copying db_mOTU_test')
    print('#####################################################')

    nomin_db_mOTU_test_folder = nomin1_output_folder + '/db_mOTU_test'

    if os.path.exists(nomin_db_mOTU_test_folder):
        shutil.rmtree(nomin_db_mOTU_test_folder)
    shutil.copytree(db_mOTU_test_folder, nomin_db_mOTU_test_folder)






    print('#####################################################')
    print('Processing db_mOTU_genes_length_NR')
    print('#####################################################')


    min1_db_mOTU_genes_length_NR_file = min1_output_folder + '/db_mOTU_genes_length_NR'
    nomin1_db_mOTU_genes_length_NR_file = nomin1_output_folder + '/db_mOTU_genes_length_NR'

    min = 0
    nomin = 0
    with open(db_mOTU_genes_length_NR_file) as inhandle, open(min1_db_mOTU_genes_length_NR_file, 'w') as min1_handle, open(nomin1_db_mOTU_genes_length_NR_file, 'w') as nomin1_handle:
        for line in inhandle:
            splits = line.strip().split()
            if splits[0] in min1_genes:
                min += 1
                min1_handle.write(f'{line.strip()}\n')
            else:
                nomin += 1
                nomin1_handle.write(f'{line.strip()}\n')

    print(f'unassigned MGs written: {min}')
    print(f'ref-mOTU/meta-mOTU MGs written: {nomin}')


    print('#####################################################')
    print('Processing db_mOTU_DB_CEN.fasta.annotations')
    print('#####################################################')



    min1_db_mOTU_DB_CEN_fastaannotations_file = min1_output_folder + '/db_mOTU_DB_CEN.fasta.annotations'
    nomin1_db_mOTU_DB_CEN_fastaannotations_file = nomin1_output_folder + '/db_mOTU_DB_CEN.fasta.annotations'

    min = 0
    nomin = 0
    with open(db_mOTU_DB_CEN_fastaannotations_file) as inhandle, open(min1_db_mOTU_DB_CEN_fastaannotations_file, 'w') as min1_handle, open(nomin1_db_mOTU_DB_CEN_fastaannotations_file, 'w') as nomin1_handle:
        for line in inhandle:
            splits = line.strip().split()
            if splits[1].startswith('unassigned.'):
                min += 1
                min1_handle.write(f'{line.strip()}\n')
            else:
                nomin += 1
                nomin1_handle.write(f'{line.strip()}\n')

    print(f'unassigned MGs written: {min}')
    print(f'ref-mOTU/meta-mOTU MGs written: {nomin}')








    print('#####################################################')
    print('Processing db_mOTU_bam_header_CEN')
    print('#####################################################')

    min1_db_mOTU_bam_header_CEN_file = min1_output_folder + '/db_mOTU_bam_header_CEN'
    nomin1_db_mOTU_bam_header_CEN_file = nomin1_output_folder + '/db_mOTU_bam_header_CEN'

    min = 0
    nomin = 0
    with open(db_mOTU_bam_header_CEN_file) as inhandle, open(min1_db_mOTU_bam_header_CEN_file,'w') as min1_handle, open(nomin1_db_mOTU_bam_header_CEN_file, 'w') as nomin1_handle:
        for line in inhandle:
            splits = line.strip().split()
            if splits[1].startswith('SN:unassigned.'):
                min += 1
                min1_handle.write(f'{line.strip()}\n')
            else:
                nomin += 1
                nomin1_handle.write(f'{line.strip()}\n')

    print(f'unassigned MGs written: {min}')
    print(f'ref-mOTU/meta-mOTU MGs written: {nomin}')



    print('#####################################################')
    print('Processing db_mOTU_bam_header_NR')
    print('#####################################################')

    min1_db_mOTU_bam_header_NR_file = min1_output_folder + '/db_mOTU_bam_header_NR'
    nomin1_db_mOTU_bam_header_NR_file = nomin1_output_folder + '/db_mOTU_bam_header_NR'

    min = 0
    nomin = 0
    with open(db_mOTU_bam_header_NR_file) as inhandle, open(min1_db_mOTU_bam_header_NR_file,'w') as min1_handle, open(nomin1_db_mOTU_bam_header_NR_file, 'w') as nomin1_handle:
        for line in inhandle:
            splits = line.strip().split()
            if splits[1].replace('SN:', '') in min1_genes:
                min += 1
                min1_handle.write(f'{line.strip()}\n')
            else:
                nomin += 1
                nomin1_handle.write(f'{line.strip()}\n')

    print(f'unassigned MGs written: {min}')
    print(f'ref-mOTU/meta-mOTU MGs written: {nomin}')






    print('####################################')
    print('### Processing db_mOTU_DB_NR.fasta ####')
    print('####################################')

    min1_nr_fasta_file = min1_output_folder + '/db_mOTU_DB_NR.fasta'
    min1_nr_fasta_unpadded_file = min1_output_folder + '/db_mOTU_DB_NR_unpadded.fasta'
    nomin1_nr_fasta_file = nomin1_output_folder + '/db_mOTU_DB_NR.fasta'

    min = 0
    nomin = 0

    with open(nr_fasta_file) as inhandle, open(min1_nr_fasta_file, 'w') as min1_handle, open(nomin1_nr_fasta_file, 'w') as nomin1_handle, open(min1_nr_fasta_unpadded_file, 'w') as min1_unpadded_handle:
        for cnt, (header, sequence) in enumerate(FastaIO.SimpleFastaParser(inhandle), 1):
            if cnt % 50000 == 0:

                print(f'{cnt}/{number_of_marker_genes} MGs read')
            if header in min1_genes:
                min += 1
                min1_handle.write(f'>{header}\n{sequence}\n')
                subseq = sequence[min1_mgcs_2_length[header][0]:min1_mgcs_2_length[header][1] + 1]
                min1_unpadded_handle.write(f'>{header}\n{subseq}\n')
            else:
                nomin += 1
                nomin1_handle.write(f'>{header}\n{sequence}\n')
    print(f'unassigned MGs written: {min}')
    print(f'ref-mOTU/meta-mOTU MGs written: {nomin}')



    print('####################################')
    print('### Processing db_mOTU_DB_CEN.fasta ####')
    print('####################################')

    min1_cen_fasta_file = min1_output_folder + '/db_mOTU_DB_CEN.fasta'
    nomin1_cen_fasta_file = nomin1_output_folder + '/db_mOTU_DB_CEN.fasta'

    min = 0
    nomin = 0

    with open(db_mOTU_DB_CEN_file) as inhandle, open(min1_cen_fasta_file, 'w') as min1_handle, open(nomin1_cen_fasta_file, 'w') as nomin1_handle:
        for cnt, (header, sequence) in enumerate(FastaIO.SimpleFastaParser(inhandle), 1):
            if cnt % 50000 == 0:
                print(f'{cnt}/163789 MGs read')
            if header.startswith('unassigned.'):
                min += 1
                min1_handle.write(f'>{header}\n{sequence}\n')
            else:
                nomin += 1
                nomin1_handle.write(f'>{header}\n{sequence}\n')
    print(f'unassigned MGs written: {min}')
    print(f'ref-mOTU/meta-mOTU MGs written: {nomin}')



    print(f'Start building vsearch database on {nomin1_nr_fasta_file}')
    subprocess.check_call(f'vsearch --makeudb_usearch {nomin1_nr_fasta_file} --output {nomin1_nr_fasta_file}.udb', shell=True)
    print(f'Finished building vsearch database on {nomin1_nr_fasta_file}')


    print('####################################')
    print('### Finished mOTUs DB preparation ##')
    print('####################################')


    
if __name__ == '__main__':
    status = main()
    sys.exit(status)    
