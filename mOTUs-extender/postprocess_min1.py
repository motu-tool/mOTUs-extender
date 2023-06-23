

import sys
import argparse
import subprocess
import glob
import Bio.SeqIO.FastaIO as FastaIO
import pathlib
import collections
#Example
#


def main():
    parser = argparse.ArgumentParser(description='Postprocess the new motus database by adding -1 MGs that are distant enough to new ext-MGs', add_help = True)
    parser.add_argument('db', action="store", help='The updated database folder (same as in prepare_mOTUs.py)')
    parser.add_argument('reffasta', action="store", help='The padded fasta file with sequences that will be extended.')
    parser.add_argument('newdbfolder', action="store", help='The folder with the file with the new ext-MGs')
    parser.add_argument('cutoffs', action="store", help='10MGs and their cutoff as a file')
    parser.add_argument('threads', action="store", help='threads')
    args = parser.parse_args()

    db_folder = args.db
    newdbfolder = args.newdbfolder
    threads = int(args.threads)
    cutoffs_file = args.cutoffs
    maxaccepts = 1000
    maxrejects = 1000
    min_coverage = 0.8

    ext_mgs_file = args.reffasta

    min1_vs_ext_mgs_file = newdbfolder + '/min1_vs_padded.m8'

    min1_unpadded_fasta_file = db_folder + '/min1_db_mOTU/db_mOTU_DB_NR_unpadded.fasta'

    #if not pathlib.Path(min1_vs_ext_mgs_file).exists():
    subprocess.check_call(f'vsearch --threads {threads} --usearch_global {min1_unpadded_fasta_file} --db {ext_mgs_file} --id 0.8  --maxaccepts {maxaccepts} --maxrejects {maxrejects} --mincols 20 --alnout {min1_vs_ext_mgs_file}.aln --blast6out {min1_vs_ext_mgs_file}',shell=True)


    mg_2_cutoff = {}
    with open(cutoffs_file) as handle:
        for line in handle:
            splits = line.strip().split(',')
            mg_2_cutoff[splits[0]] = float(splits[1])


    min1_genes_to_remove = set()



    with open(min1_vs_ext_mgs_file) as handle:
        for line in handle:
            splits = line.strip().split()

            query = splits[0]
            cog = query.split('.')[1].split('-')[0]#query.rsplit('.', 1)[-1]
            cutoff = mg_2_cutoff[cog]
            ref = splits[1]
            percid = float(splits[2])
            alength = float(splits[3])
            qlength = float(splits[7])
            reflength = float(splits[9])
            qcov =  alength / qlength
            rcov =  alength / reflength
            acov = qcov
            if rcov > acov:
                acov = rcov
            if acov < min_coverage:
                continue

            if percid < cutoff:
                continue

            min1_genes_to_remove.add(query)





    '''
    Update the files
    '''
    inputfolder = db_folder + '/min1_db_mOTU/'
    outputfolder = newdbfolder + '/updated_min1_db_mOTU/'

    mgc_2_gene_file = inputfolder + '/db_mOTU_MAP_genes_to_MGCs.tsv'
    mgc_2_motus_file = inputfolder + '/db_mOTU_MAP_MGCs_to_mOTUs.tsv'
    db_mOTU_MAP_MGCs_to_mOTUs_inline_file = inputfolder + '/db_mOTU_MAP_MGCs_to_mOTUs_in-line.tsv'
    db_mOTU_padding_coordinates_NR_file = inputfolder + '/db_mOTU_padding_coordinates_NR.tsv'
    db_mOTU_genes_length_NR_file = inputfolder + '/db_mOTU_genes_length_NR'

    nr_fasta_file = inputfolder + '/db_mOTU_DB_NR.fasta'
    db_mOTU_padding_coordinates_CEN_file = inputfolder + '/db_mOTU_padding_coordinates_CEN.tsv'
    db_mOTU_DB_CEN_fastaannotations_file = inputfolder + '/db_mOTU_DB_CEN.fasta.annotations'



    db_mOTU_DB_CEN_file = inputfolder + '/db_mOTU_DB_CEN.fasta'

    '''
    Not needed motus files for -1
    '''
    #db_mOTU_taxonomy_refmOTUs_short_names_file = inputfolder + '/db_mOTU_taxonomy_ref-mOTUs_short_names.tsv'
    #db_mOTU_taxonomy_CAMI_file = inputfolder + '/db_mOTU_taxonomy_CAMI.tsv'
    #db_mOTU_taxonomy_refmOTUs_file = inputfolder + '/db_mOTU_taxonomy_ref-mOTUs.tsv'
    #README_file = inputfolder + '/README'
    #db_mOTU_bam_header_NR_file = inputfolder + '/db_mOTU_bam_header_NR'
    #db_mOTU_bam_header_CEN_file = inputfolder + '/db_mOTU_bam_header_CEN'
    #db_mOTU_taxonomy_metamOTUs_file = inputfolder + '/db_mOTU_taxonomy_meta-mOTUs.tsv'
    #db_mOTU_test_folder = inputfolder + '/db_mOTU_test'



    pathlib.Path(outputfolder).mkdir(parents=True, exist_ok=True)




    mgc_2_mgs = collections.defaultdict(set)
    mgcs_to_remove = set()
    mgs_to_remove = set()

    # collecting mgcs and mgs to remove

    with open(mgc_2_gene_file) as inhandle:
        for line in inhandle:
            splits = line.strip().split()
            mg = splits[0]
            mgc = splits[3]
            mgc_2_mgs[mgc].add(mg)
    for mgc, mgs in mgc_2_mgs.items():
        if len(min1_genes_to_remove.intersection(mgs)) != 0:
            #print(f'{mgc}\t{len(min1_genes_to_remove.intersection(mgs))}\t{len(mgs)}')
            mgcs_to_remove.add(mgc)
            for mg in mgs:
                mgs_to_remove.add(mg)

    print(f'{len(mgcs_to_remove)} MGCs will be removed from the -1 mOTU')
    print(f'{len(mgs_to_remove)} MGs will be removed from the -1 mOTU')

    min1_genes_to_remove = None # just to find errors



    updated_mgc_2_gene_file = outputfolder + 'db_mOTU_MAP_genes_to_MGCs.tsv'
    print(f'Writing {updated_mgc_2_gene_file}')

    all_mgs = 0
    rm_mgs = 0
    with open(mgc_2_gene_file) as inhandle, open(updated_mgc_2_gene_file,'w') as outhandle:
        for line in inhandle:
            all_mgs += 1
            splits = line.strip().split()
            mg = splits[0]
            mgc = splits[3]

            if splits[0] not in mgs_to_remove:
                outhandle.write(f'{line}')
                continue
            rm_mgs += 1
    print(f'Wrote {all_mgs - rm_mgs} MGs. Removed {rm_mgs} MGs')
    print('################################################################')



    updated_mgc_2_motus_file = outputfolder + 'db_mOTU_MAP_MGCs_to_mOTUs.tsv'
    print(f'Writing {updated_mgc_2_motus_file}')
    written_mgcs = 0
    with open(mgc_2_motus_file) as inhandle, open(updated_mgc_2_motus_file, 'w') as outhandle:
        for line in inhandle:
            splits = line.strip().split()
            mgc = splits[0]
            if mgc in mgcs_to_remove:
                continue
            else:
                written_mgcs += 1
                outhandle.write(line)

    print(f'Wrote {written_mgcs} MGCs')
    print('################################################################')






    updated_mOTU_MAP_MGCs_to_mOTUs_inline_file = outputfolder + 'db_mOTU_MAP_MGCs_to_mOTUs_in-line.tsv'
    print(f'Writing {updated_mOTU_MAP_MGCs_to_mOTUs_inline_file}')
    with open(db_mOTU_MAP_MGCs_to_mOTUs_inline_file) as inhandle, open(updated_mOTU_MAP_MGCs_to_mOTUs_inline_file, 'w') as outhandle:
        mgcs = []
        for line in inhandle:
            allmgcs = line.strip().split()[1].split(';')
            mgcs = ';'.join([mgc for mgc in allmgcs if mgc not in mgcs_to_remove])
        outhandle.write(f'unassigned\t{mgcs}')

    print('################################################################')














    updated_mOTU_padding_coordinates_NR_file = outputfolder + 'db_mOTU_padding_coordinates_NR.tsv'
    written_mgs = 0
    removed_mgs = 0
    print(f'Writing {updated_mOTU_padding_coordinates_NR_file}')
    with open(db_mOTU_padding_coordinates_NR_file) as inhandle, open(updated_mOTU_padding_coordinates_NR_file, 'w') as outhandle:
        for line in inhandle:
            splits = line.strip().split()
            mg = splits[0]
            if mg not in mgs_to_remove:
                written_mgs += 1
                outhandle.write(line)
            else:
                removed_mgs += 1
    print(f'Written MGs {written_mgs}\nRemoved MGs {removed_mgs}')
    print('################################################################')









    updated_mOTU_genes_length_NR_file = outputfolder + 'db_mOTU_genes_length_NR'
    written_mgs = 0
    removed_mgs = 0
    print(f'Writing {updated_mOTU_genes_length_NR_file}')
    with open(db_mOTU_genes_length_NR_file) as inhandle, open(updated_mOTU_genes_length_NR_file, 'w') as outhandle:
        for line in inhandle:
            splits = line.strip().split()
            mg = splits[0]
            if mg not in mgs_to_remove:
                written_mgs += 1
                outhandle.write(line)
            else:
                removed_mgs += 1
    print(f'Written MGs {written_mgs}\nRemoved MGs {removed_mgs}')
    print('################################################################')




    updated_nr_fasta_file = outputfolder + 'db_mOTU_DB_NR.fasta'
    print(f'Writing {updated_nr_fasta_file}')
    written_mgs = 0
    removed_mgs = 0


    with open(nr_fasta_file) as inhandle, open(updated_nr_fasta_file, 'w') as outhandle:
        for cnt, (header, sequence) in enumerate(FastaIO.SimpleFastaParser(inhandle), 1):

            if header not in mgs_to_remove:
                written_mgs += 1
                outhandle.write(f'>{header}\n{sequence}\n')
            else:
                removed_mgs += 1
    print(f'Written MGs {written_mgs}\nRemoved MGs {removed_mgs}')
    print('################################################################')









    updated_db_mOTU_padding_coordinates_CEN_file = outputfolder + '/db_mOTU_padding_coordinates_CEN.tsv'
    print(f'Writing {updated_db_mOTU_padding_coordinates_CEN_file}')
    written_mgs = 0
    removed_mgs = 0
    with open(db_mOTU_padding_coordinates_CEN_file) as inhandle, open(updated_db_mOTU_padding_coordinates_CEN_file, 'w') as outhandle:
        for line in inhandle:
            splits = line.strip().split()
            mg = splits[0]#.replace('NA.', '')
            mg_orig_name = mg[::-1].replace('_', '.', 1)[::-1]
            if mg_orig_name not in mgs_to_remove:
                outhandle.write(line)
                written_mgs += 1
            else:
                removed_mgs += 1

    print(f'Written MGs {written_mgs}\nRemoved MGs {removed_mgs} --> Numbers are different to the previous ones as these are centroid sequences')
    print('################################################################')






    updated_db_mOTU_DB_CEN_fastaannotations_file = outputfolder + '/db_mOTU_DB_CEN.fasta.annotations'

    print(f'Writing {updated_db_mOTU_DB_CEN_fastaannotations_file}')
    written_mgs = 0
    removed_mgs = 0
    with open(db_mOTU_DB_CEN_fastaannotations_file) as inhandle, open(updated_db_mOTU_DB_CEN_fastaannotations_file, 'w') as outhandle:
        for line in inhandle:
            splits = line.strip().split()
            mg = splits[1]#.replace('NA.', '')
            mg_orig_name = mg[::-1].replace('_', '.', 1)[::-1]#mg.replace('_C', '.C')
            if mg_orig_name not in mgs_to_remove:
                outhandle.write(line)
                written_mgs += 1
            else:
                removed_mgs += 1
    print(f'Written MGs {written_mgs}\nRemoved MGs {removed_mgs} --> Numbers are different to the previous ones as these are centroid sequences')
    print('################################################################')






    updated_db_mOTU_DB_CEN_file = outputfolder + '/db_mOTU_DB_CEN.fasta'
    print(f'Reading {db_mOTU_DB_CEN_file}')
    print(f'Writing {updated_db_mOTU_DB_CEN_file}')
    written_mgs = 0
    removed_mgs = 0
    with open(db_mOTU_DB_CEN_file) as inhandle, open(updated_db_mOTU_DB_CEN_file, 'w') as outhandle:
        for cnt, (header, sequence) in enumerate(FastaIO.SimpleFastaParser(inhandle), 1):
            mg = header#.replace('NA.', '')
            mg_orig_name = mg[::-1].replace('_', '.', 1)[::-1]#mg.replace('_C', '.C')
            if mg_orig_name not in mgs_to_remove:
                written_mgs += 1
                outhandle.write(f'>{header}\n{sequence}\n')
            else:

                removed_mgs += 1

    print(f'Written MGs {written_mgs}\nRemoved MGs {removed_mgs}')
    print('################################################################')


    print('Succesfully finished')


if __name__ == '__main__':
    status = main()
    sys.exit(status)

