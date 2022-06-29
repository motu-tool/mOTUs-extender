#!/usr/bin/env python

import os
import sys
import argparse

# example call:
# python create_cami_tax.py motusfolder/db_mOTU/CAMI_Sep_2015_mOTUs_2.0.0.tsv  tempfolder/dbname/db_mOTU/mOTULG.taxonomy  tempfolder/dbname/db_mOTU/CAMI_Sep_2015_mOTUs_2.0.0.tsv

parser = argparse.ArgumentParser(description='This program creates the CAMI taxonomy for augmenting the mOTUs db', add_help=True)
parser.add_argument('old_file', action="store", help='old cami file in the original mOTUs directory')
parser.add_argument('augmented_meta_mOTU_tax', action="store", help='file mOTULG.taxonomy with the dbname_* at the end')
parser.add_argument('new_file', action="store", help='new cami file with the information of the dbname_*')
parser.add_argument('db_name', action="store", help='name of the new entries in the augmented_meta_mOTU_tax',default="newdbname")
args = parser.parse_args()

# open files
original_f = open(args.old_file, "r")
new_f = open(args.new_file, "w")

# write old file into new file
for i in original_f:
    if len(i) > 5: # this is to avoid the last "\n"
        new_f.write(i)

original_f.close()

# find the new lines that have to be added
tax_f = open(args.augmented_meta_mOTU_tax, "r")
lines_to_add = list()

for i in tax_f:
    if i.startswith(args.db_name):
        lines_to_add.append(i.rstrip())

tax_f.close()

# parse the lines that have to be added and write them in the new file
for i in lines_to_add:
    # parsed string
    final_s = ""

    # start parsing
    all_v = i.split("\t")

    # first value is the name of the motus
    final_s = final_s + all_v[0]

    # start to put info of the taxonomy levels - superkingdom
    ncbi_id = all_v[1].split(" ")[0]
    name = " ".join(all_v[1].split(" ")[1:])
    final_s = final_s + "\t" + ncbi_id + "\tsuperkingdom\t" + ncbi_id + "\t" + name

    # phylum
    k = 2
    ncbi_id_alone = all_v[k].split(" ")[0]
    ncbi_id = ncbi_id + "|" + ncbi_id_alone
    name = name + "|" + (" ".join(all_v[k].split(" ")[1:]) )
    final_s = final_s + "\t" + ncbi_id_alone + "\tphylum\t" + ncbi_id + "\t" + name

    # class
    k = 3
    ncbi_id_alone = all_v[k].split(" ")[0]
    ncbi_id = ncbi_id + "|" + ncbi_id_alone
    name = name + "|" + (" ".join(all_v[k].split(" ")[1:]) )
    final_s = final_s + "\t" + ncbi_id_alone + "\tclass\t" + ncbi_id + "\t" + name

    # order
    k = 4
    ncbi_id_alone = all_v[k].split(" ")[0]
    ncbi_id = ncbi_id + "|" + ncbi_id_alone
    name = name + "|" + (" ".join(all_v[k].split(" ")[1:]) )
    final_s = final_s + "\t" + ncbi_id_alone + "\torder\t" + ncbi_id + "\t" + name

    # family
    k = 5
    ncbi_id_alone = all_v[k].split(" ")[0]
    ncbi_id = ncbi_id + "|" + ncbi_id_alone
    name = name + "|" + (" ".join(all_v[k].split(" ")[1:]) )
    final_s = final_s + "\t" + ncbi_id_alone + "\tfamily\t" + ncbi_id + "\t" + name

    # genus
    k = 6
    ncbi_id_alone = all_v[k].split(" ")[0]
    ncbi_id = ncbi_id + "|" + ncbi_id_alone
    name = name + "|" + (" ".join(all_v[k].split(" ")[1:]) )
    final_s = final_s + "\t" + ncbi_id_alone + "\tgenus\t" + ncbi_id + "\t" + name

    # species
    k = 7
    ncbi_id_alone = all_v[k].split(" ")[0]
    ncbi_id = ncbi_id + "|" + ncbi_id_alone
    name = name + "|" + (" ".join(all_v[k].split(" ")[1:]) )
    final_s = final_s + "\t" + ncbi_id_alone + "\tspecies\t" + ncbi_id + "\t" + name

    # write the string
    new_f.write(final_s+"\n")

new_f.close()
