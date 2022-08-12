![alt text](https://raw.githubusercontent.com/motu-tool/mOTUs/master/pics/motu_logo.png)




mOTU extender
========

The mOTU profiler is a computational tool that estimates relative taxonomic abundance of known and currently unknown microbial community members using metagenomic shotgun sequencing data. The **mOTU extender** is a companion tool used to extend the current mOTUs database with user provided genomes.



If you are using mOTUs, please cite:

> Hans-Joachim Ruscheweyh, Alessio Milanese, Lucas Paoli, Nicolai Karcher, Quentin Clayssen, Marisa Isabell Metzger, Jakob Wirbel, Peer Bork, Daniel R. Mende, Georg Zeller & Shinichi Sunagawa.
> **Reference genome-independent taxonomic profiling of microbiomes with mOTUs3**; 
> BioRxiv, doi: [0.1101/2021.04.20.440600](https://doi.org/10.1101/2021.04.20.440600)


Pre-requisites
--------------
The mOTU extender was tested with:

* `python  3.6.7`
	* `biopython 1.79`
	* `numpy 1.18.2`
	* `scipy 1.1.0`
* `bwa 0.7.17-r1188`
* `samtools 1.9`
* `snakemake 6.15.5`
* `vsearch 2.14.1`
* `prodigal 2.6.3`
* `fetchMGs 1.2`


Installation
--------------

The mOTUs extender can be installed via GitHub:

```
git clone https://github.com/motu-tool/mOTUs-extender.git
```

Basic example
-------------

To use the mOTU extender you will need to provide 2 input data:
- A folder with genomes, one file per genome with the `.fa` file suffix.
- A 7 rank taxonomy (NCBI) file which contains the taxonomic assignation of each genome (see here for an example: `mOTUs-extender/mOTUs-extender/test/genomes.tax`)

We provide a test dataset with 50 genomes taken from ([Buck et. al, 2021](https://doi.org/10.1038/s41597-021-00910-1)).

The mOTUs extender is a 3 step pipeline which (1) prepares the current database, (2) checks what user genomes are already part of the current database and (3) adds genomes with novel taxonomic information to the database.

### 0. Test Data

``` bash
mkdir -p test_extension/genomes/
cp mOTUs-extender/mOTUs-extender/test/genomes/*gz test_extension/genomes/
gunzip test_extension/genomes/*gz
```

### 1. mOTUs extender - prepare

```
python mOTUs-extender/mOTUs-extender/motus-extender.py prepare -w test_extension/extension
```

**Note** This command will download the current mOTUs database (~3.5GB).

<details>
<summary>Commandline stdout</summary>
<br>

```bash

python mOTUs-extender/mOTUs-extender/motus-extender.py prepare -w test_extension/extension
2022-08-12 10:17:50,392 INFO: Starting mOTUs-extender
2022-08-12 10:17:50,393 INFO: Creating work folder test_extension/extension
2022-08-12 10:17:50,395 INFO: Downloading mOTUs 3.0.1 database
2022-08-12 10:17:50,396 INFO: Executing: curl https://zenodo.org/record/5140350/files/db_mOTU_v3.0.1.tar.gz -o test_extension/extension/orig_db/db_mOTU_v3.0.1.tar.gz
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 3342M  100 3342M    0     0  29.7M      0  0:01:52  0:01:52 --:--:-- 35.8M
2022-08-12 10:19:42,904 INFO: Uncompressing mOTUs 3.0.1 database
2022-08-12 10:19:42,904 INFO: Executing: tar -xzvf test_extension/extension/orig_db/db_mOTU_v3.0.1.tar.gz -C test_extension/extension/orig_db
db_mOTU/
db_mOTU/db_mOTU_bam_header_NR
db_mOTU/db_mOTU_padding_coordinates_NR.tsv
db_mOTU/db_mOTU_MAP_MGCs_to_mOTUs.tsv
db_mOTU/db_mOTU_taxonomy_ref-mOTUs.tsv
db_mOTU/db_mOTU_taxonomy_CAMI.tsv
db_mOTU/db_mOTU_DB_CEN.fasta.ann
db_mOTU/db_mOTU_DB_CEN.fasta.bwt
db_mOTU/db_mOTU_versions
db_mOTU/public_profiles/
db_mOTU/public_profiles/mOTUs.profiles_environments.gz
db_mOTU/public_profiles/mOTUs.profiles.gz
db_mOTU/db_mOTU_padding_coordinates_CEN.tsv
db_mOTU/db_mOTU_DB_NR.fasta.bwt
db_mOTU/db_mOTU_DB_NR.fasta.ann
db_mOTU/db_mOTU_bam_header_CEN
db_mOTU/db_mOTU_taxonomy_ref-mOTUs_short_names.tsv
db_mOTU/db_mOTU_MAP_genes_to_MGCs.tsv
db_mOTU/README
db_mOTU/db_mOTU_taxonomy_meta-mOTUs.tsv
db_mOTU/db_mOTU_DB_CEN.fasta
db_mOTU/db_mOTU_test/
db_mOTU/db_mOTU_test/test1.motus
db_mOTU/db_mOTU_test/test1_single.fastq
db_mOTU/db_mOTU_DB_NR.fasta.amb
db_mOTU/db_mOTU_DB_NR.fasta.sa
db_mOTU/db_mOTU_DB_NR.fasta.pac
db_mOTU/db_mOTU_DB_CEN.fasta.annotations
db_mOTU/db_mOTU_DB_CEN.fasta.pac
db_mOTU/db_mOTU_DB_CEN.fasta.sa
db_mOTU/db_mOTU_DB_NR.fasta
db_mOTU/db_mOTU_genes_length_NR
db_mOTU/db_mOTU_DB_CEN.fasta.amb
db_mOTU/db_mOTU_MAP_MGCs_to_mOTUs_in-line.tsv
2022-08-12 10:21:07,307 INFO: Preparing mOTUs 3.0.1 database for extension
2022-08-12 10:21:07,308 INFO: Executing: python mOTUs-extender/mOTUs-extender/prepare_mOTUs.py test_extension/extension/orig_db/db_mOTU/ test_extension/extension/temp_db_folder
####################################
#### Checking db_mOTU_versions #####
####################################
Version of motus is 2.6.0. Correct version. Continue
Version of nr is 2.6.0. Correct version. Continue
Version of cen is 2.6.0. Correct version. Continue
Version of append is 2.6.0. Correct version. Continue
Version of map_genes_to_mOTUs is 2.6.0. Correct version. Continue
Version of map_mOTUs_to_LGs is 2.6.0. Correct version. Continue
Version of runBWA is 2.6.0. Correct version. Continue
Version of specI_tax is 2.6.0. Correct version. Continue
Version of mOTULG_tax is 2.6.0. Correct version. Continue
###################################################
##### Processing db_mOTU_MAP_MGCs_to_mOTUs.tsv ####
###################################################
Found 521780 MGCs of which 213652 MGCs belong to unassigned
###################################################
#### Processing db_mOTU_MAP_genes_to_MGCs.tsv #####
###################################################
Found 1184175 MGs of which 222782 MGs belong to unassigned
####################################################
# Processing db_mOTU_MAP_MGCs_to_mOTUs_in-line.tsv #
####################################################
unassigned mOTU written: 1
ref-mOTU/meta-mOTU mOTUs written: 33570
####################################################
## Processing db_mOTU_padding_coordinates_NR.tsv ###
####################################################
-1 MGs written: 222782
ref-mOTU/meta-mOTU MGs written: 961393
####################################################
## Processing db_mOTU_padding_coordinates_CEN.tsv ##
####################################################
unassigned MGs written: 213652
ref-mOTU/meta-mOTU MGs written: 308128
#####################################################
Copying db_mOTU_taxonomy_ref-mOTUs_short_names.tsv
#####################################################
#####################################################
Copying db_mOTU_taxonomy_CAMI.tsv
#####################################################
#####################################################
Copying db_mOTU_taxonomy_ref-mOTUs.tsv
#####################################################
#####################################################
Copying README
#####################################################
#####################################################
Copying db_mOTU_taxonomy_meta-mOTUs.tsv
#####################################################
#####################################################
Copying db_mOTU_test
#####################################################
#####################################################
Processing db_mOTU_genes_length_NR
#####################################################
unassigned MGs written: 222782
ref-mOTU/meta-mOTU MGs written: 961393
#####################################################
Processing db_mOTU_DB_CEN.fasta.annotations
#####################################################
unassigned MGs written: 213652
ref-mOTU/meta-mOTU MGs written: 308128
#####################################################
Processing db_mOTU_bam_header_CEN
#####################################################
unassigned MGs written: 213652
ref-mOTU/meta-mOTU MGs written: 308128
#####################################################
Processing db_mOTU_bam_header_NR
#####################################################
unassigned MGs written: 222782
ref-mOTU/meta-mOTU MGs written: 961393
####################################
### Processing db_mOTU_DB_NR.fasta ####
####################################
50000/1184175 MGs read
100000/1184175 MGs read
150000/1184175 MGs read
200000/1184175 MGs read
250000/1184175 MGs read
300000/1184175 MGs read
350000/1184175 MGs read
400000/1184175 MGs read
450000/1184175 MGs read
500000/1184175 MGs read
550000/1184175 MGs read
600000/1184175 MGs read
650000/1184175 MGs read
700000/1184175 MGs read
750000/1184175 MGs read
800000/1184175 MGs read
850000/1184175 MGs read
900000/1184175 MGs read
950000/1184175 MGs read
1000000/1184175 MGs read
1050000/1184175 MGs read
1100000/1184175 MGs read
1150000/1184175 MGs read
unassigned MGs written: 222782
ref-mOTU/meta-mOTU MGs written: 961393
####################################
### Processing db_mOTU_DB_CEN.fasta ####
####################################
50000/163789 MGs read
100000/163789 MGs read
150000/163789 MGs read
200000/163789 MGs read
250000/163789 MGs read
300000/163789 MGs read
350000/163789 MGs read
400000/163789 MGs read
450000/163789 MGs read
500000/163789 MGs read
unassigned MGs written: 213652
ref-mOTU/meta-mOTU MGs written: 308128
Start building vsearch database on test_extension/extension/temp_db_folder/db_mOTU//db_mOTU_DB_NR.fasta
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/temp_db_folder/db_mOTU//db_mOTU_DB_NR.fasta 100%
1647961589 nt in 961393 seqs, min 403, max 9725, avg 1714
Masking 100%
Counting k-mers 100%
Creating k-mer index 100%
Writing UDB file 100%
Finished building vsearch database on test_extension/extension/temp_db_folder/db_mOTU//db_mOTU_DB_NR.fasta
####################################
### Finished mOTUs DB preparation ##
####################################
2022-08-12 10:30:51,695 INFO: mOTUs-extender finished successfully

```
</details>

### 2. mOTUs extender - membership

```
python mOTUs-extender/mOTUs-extender/motus-extender.py membership -w test_extension/extension/ -g test_extension/genomes/ -t 32
```

<details>
<summary>Commandline stdout</summary>
<br>

```
python mOTUs-extender/mOTUs-extender/motus-extender.py membership -w test_extension/extension/ -g test_extension/genomes/ -t 32
2022-08-12 11:20:36,781 INFO: Starting mOTUs-extender
2022-08-12 11:20:36,783 INFO: Executing: gunzip -c mOTUs-extender/mOTUs-extender/mOTU.v3.0.mOTU-LG.map.gz > test_extension/extension/mOTU.v3.0.mOTU-LG.map
2022-08-12 11:20:37,188 INFO: Found 50 genomes in test_extension/genomes path
2022-08-12 11:20:37,250 INFO: Executing: snakemake -s mOTUs-extender/mOTUs-extender/0prod_fetchMGs.py --config mapfile=test_extension/extension/mOTU.v3.0.mOTU-LG.map infolder=test_extension/genomes/ outfolder=test_extension/extension/ temp_db_folder=test_extension/extension/temp_db_folder/ scriptfolder=mOTUs-extender/mOTUs-extender genomes=test_extension/genomes/GCA_903824045.1.fa,test_extension/genomes/GCA_903824635.1.fa,test_extension/genomes/GCA_903824795.1.fa,test_extension/genomes/GCA_903826635.1.fa,test_extension/genomes/GCA_903835685.1.fa,test_extension/genomes/GCA_903839445.1.fa,test_extension/genomes/GCA_903841135.1.fa,test_extension/genomes/GCA_903842315.1.fa,test_extension/genomes/GCA_903843765.1.fa,test_extension/genomes/GCA_903845665.1.fa,test_extension/genomes/GCA_903853455.1.fa,test_extension/genomes/GCA_903853495.1.fa,test_extension/genomes/GCA_903854225.1.fa,test_extension/genomes/GCA_903854255.1.fa,test_extension/genomes/GCA_903857725.1.fa,test_extension/genomes/GCA_903869175.1.fa,test_extension/genomes/GCA_903869725.1.fa,test_extension/genomes/GCA_903871585.1.fa,test_extension/genomes/GCA_903875175.1.fa,test_extension/genomes/GCA_903879295.1.fa,test_extension/genomes/GCA_903883715.1.fa,test_extension/genomes/GCA_903884585.1.fa,test_extension/genomes/GCA_903888055.1.fa,test_extension/genomes/GCA_903894295.1.fa,test_extension/genomes/GCA_903895255.1.fa,test_extension/genomes/GCA_903899515.1.fa,test_extension/genomes/GCA_903900705.1.fa,test_extension/genomes/GCA_903902845.1.fa,test_extension/genomes/GCA_903902895.1.fa,test_extension/genomes/GCA_903907865.1.fa,test_extension/genomes/GCA_903907955.1.fa,test_extension/genomes/GCA_903913555.1.fa,test_extension/genomes/GCA_903914905.1.fa,test_extension/genomes/GCA_903915725.1.fa,test_extension/genomes/GCA_903915855.1.fa,test_extension/genomes/GCA_903918315.1.fa,test_extension/genomes/GCA_903920605.1.fa,test_extension/genomes/GCA_903920995.1.fa,test_extension/genomes/GCA_903924625.1.fa,test_extension/genomes/GCA_903928895.1.fa,test_extension/genomes/GCA_903934285.1.fa,test_extension/genomes/GCA_903935095.1.fa,test_extension/genomes/GCA_903936775.1.fa,test_extension/genomes/GCA_903937725.1.fa,test_extension/genomes/GCA_903939965.1.fa,test_extension/genomes/GCA_903945225.1.fa,test_extension/genomes/GCA_903945495.1.fa,test_extension/genomes/GCA_903953565.1.fa,test_extension/genomes/GCA_903959825.1.fa,test_extension/genomes/GCA_903961455.1.fa -j 32 -k
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 32
Rules claiming more threads will be scaled down.
Job stats:
job                 count    min threads    max threads
----------------  -------  -------------  -------------
all                     1              1              1
prodigalfetchMgs       50              1              1
total                  51              1              1

Select jobs to execute...

[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903920605.1.fa
    output: test_extension/extension/genes/GCA_903920605.1/GCA_903920605.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903920605.1/GCA_903920605.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903920605.1/GCA_903920605.1.prodigal_fetchMGs.command
    jobid: 37
    benchmark: test_extension/extension/genes/GCA_903920605.1/GCA_903920605.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903920605.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903857725.1.fa
    output: test_extension/extension/genes/GCA_903857725.1/GCA_903857725.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903857725.1/GCA_903857725.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903857725.1/GCA_903857725.1.prodigal_fetchMGs.command
    jobid: 15
    benchmark: test_extension/extension/genes/GCA_903857725.1/GCA_903857725.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903857725.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903939965.1.fa
    output: test_extension/extension/genes/GCA_903939965.1/GCA_903939965.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903939965.1/GCA_903939965.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903939965.1/GCA_903939965.1.prodigal_fetchMGs.command
    jobid: 45
    benchmark: test_extension/extension/genes/GCA_903939965.1/GCA_903939965.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903939965.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903842315.1.fa
    output: test_extension/extension/genes/GCA_903842315.1/GCA_903842315.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903842315.1/GCA_903842315.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903842315.1/GCA_903842315.1.prodigal_fetchMGs.command
    jobid: 8
    benchmark: test_extension/extension/genes/GCA_903842315.1/GCA_903842315.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903842315.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903888055.1.fa
    output: test_extension/extension/genes/GCA_903888055.1/GCA_903888055.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903888055.1/GCA_903888055.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903888055.1/GCA_903888055.1.prodigal_fetchMGs.command
    jobid: 23
    benchmark: test_extension/extension/genes/GCA_903888055.1/GCA_903888055.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903888055.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903920995.1.fa
    output: test_extension/extension/genes/GCA_903920995.1/GCA_903920995.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903920995.1/GCA_903920995.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903920995.1/GCA_903920995.1.prodigal_fetchMGs.command
    jobid: 38
    benchmark: test_extension/extension/genes/GCA_903920995.1/GCA_903920995.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903920995.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903869175.1.fa
    output: test_extension/extension/genes/GCA_903869175.1/GCA_903869175.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903869175.1/GCA_903869175.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903869175.1/GCA_903869175.1.prodigal_fetchMGs.command
    jobid: 16
    benchmark: test_extension/extension/genes/GCA_903869175.1/GCA_903869175.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903869175.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903907955.1.fa
    output: test_extension/extension/genes/GCA_903907955.1/GCA_903907955.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903907955.1/GCA_903907955.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903907955.1/GCA_903907955.1.prodigal_fetchMGs.command
    jobid: 31
    benchmark: test_extension/extension/genes/GCA_903907955.1/GCA_903907955.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903907955.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903824635.1.fa
    output: test_extension/extension/genes/GCA_903824635.1/GCA_903824635.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903824635.1/GCA_903824635.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903824635.1/GCA_903824635.1.prodigal_fetchMGs.command
    jobid: 2
    benchmark: test_extension/extension/genes/GCA_903824635.1/GCA_903824635.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903824635.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903945225.1.fa
    output: test_extension/extension/genes/GCA_903945225.1/GCA_903945225.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903945225.1/GCA_903945225.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903945225.1/GCA_903945225.1.prodigal_fetchMGs.command
    jobid: 46
    benchmark: test_extension/extension/genes/GCA_903945225.1/GCA_903945225.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903945225.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903843765.1.fa
    output: test_extension/extension/genes/GCA_903843765.1/GCA_903843765.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903843765.1/GCA_903843765.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903843765.1/GCA_903843765.1.prodigal_fetchMGs.command
    jobid: 9
    benchmark: test_extension/extension/genes/GCA_903843765.1/GCA_903843765.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903843765.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903945495.1.fa
    output: test_extension/extension/genes/GCA_903945495.1/GCA_903945495.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903945495.1/GCA_903945495.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903945495.1/GCA_903945495.1.prodigal_fetchMGs.command
    jobid: 47
    benchmark: test_extension/extension/genes/GCA_903945495.1/GCA_903945495.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903945495.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903902895.1.fa
    output: test_extension/extension/genes/GCA_903902895.1/GCA_903902895.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903902895.1/GCA_903902895.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903902895.1/GCA_903902895.1.prodigal_fetchMGs.command
    jobid: 29
    benchmark: test_extension/extension/genes/GCA_903902895.1/GCA_903902895.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903902895.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903826635.1.fa
    output: test_extension/extension/genes/GCA_903826635.1/GCA_903826635.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903826635.1/GCA_903826635.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903826635.1/GCA_903826635.1.prodigal_fetchMGs.command
    jobid: 4
    benchmark: test_extension/extension/genes/GCA_903826635.1/GCA_903826635.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903826635.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903875175.1.fa
    output: test_extension/extension/genes/GCA_903875175.1/GCA_903875175.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903875175.1/GCA_903875175.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903875175.1/GCA_903875175.1.prodigal_fetchMGs.command
    jobid: 19
    benchmark: test_extension/extension/genes/GCA_903875175.1/GCA_903875175.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903875175.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903915725.1.fa
    output: test_extension/extension/genes/GCA_903915725.1/GCA_903915725.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903915725.1/GCA_903915725.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903915725.1/GCA_903915725.1.prodigal_fetchMGs.command
    jobid: 34
    benchmark: test_extension/extension/genes/GCA_903915725.1/GCA_903915725.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903915725.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903959825.1.fa
    output: test_extension/extension/genes/GCA_903959825.1/GCA_903959825.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903959825.1/GCA_903959825.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903959825.1/GCA_903959825.1.prodigal_fetchMGs.command
    jobid: 49
    benchmark: test_extension/extension/genes/GCA_903959825.1/GCA_903959825.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903959825.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903853495.1.fa
    output: test_extension/extension/genes/GCA_903853495.1/GCA_903853495.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903853495.1/GCA_903853495.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903853495.1/GCA_903853495.1.prodigal_fetchMGs.command
    jobid: 12
    benchmark: test_extension/extension/genes/GCA_903853495.1/GCA_903853495.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903853495.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903900705.1.fa
    output: test_extension/extension/genes/GCA_903900705.1/GCA_903900705.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903900705.1/GCA_903900705.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903900705.1/GCA_903900705.1.prodigal_fetchMGs.command
    jobid: 27
    benchmark: test_extension/extension/genes/GCA_903900705.1/GCA_903900705.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903900705.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903935095.1.fa
    output: test_extension/extension/genes/GCA_903935095.1/GCA_903935095.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903935095.1/GCA_903935095.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903935095.1/GCA_903935095.1.prodigal_fetchMGs.command
    jobid: 42
    benchmark: test_extension/extension/genes/GCA_903935095.1/GCA_903935095.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903935095.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903879295.1.fa
    output: test_extension/extension/genes/GCA_903879295.1/GCA_903879295.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903879295.1/GCA_903879295.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903879295.1/GCA_903879295.1.prodigal_fetchMGs.command
    jobid: 20
    benchmark: test_extension/extension/genes/GCA_903879295.1/GCA_903879295.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903879295.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903835685.1.fa
    output: test_extension/extension/genes/GCA_903835685.1/GCA_903835685.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903835685.1/GCA_903835685.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903835685.1/GCA_903835685.1.prodigal_fetchMGs.command
    jobid: 5
    benchmark: test_extension/extension/genes/GCA_903835685.1/GCA_903835685.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903835685.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903915855.1.fa
    output: test_extension/extension/genes/GCA_903915855.1/GCA_903915855.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903915855.1/GCA_903915855.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903915855.1/GCA_903915855.1.prodigal_fetchMGs.command
    jobid: 35
    benchmark: test_extension/extension/genes/GCA_903915855.1/GCA_903915855.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903915855.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903961455.1.fa
    output: test_extension/extension/genes/GCA_903961455.1/GCA_903961455.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903961455.1/GCA_903961455.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903961455.1/GCA_903961455.1.prodigal_fetchMGs.command
    jobid: 50
    benchmark: test_extension/extension/genes/GCA_903961455.1/GCA_903961455.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903961455.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903854255.1.fa
    output: test_extension/extension/genes/GCA_903854255.1/GCA_903854255.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903854255.1/GCA_903854255.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903854255.1/GCA_903854255.1.prodigal_fetchMGs.command
    jobid: 14
    benchmark: test_extension/extension/genes/GCA_903854255.1/GCA_903854255.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903854255.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903854225.1.fa
    output: test_extension/extension/genes/GCA_903854225.1/GCA_903854225.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903854225.1/GCA_903854225.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903854225.1/GCA_903854225.1.prodigal_fetchMGs.command
    jobid: 13
    benchmark: test_extension/extension/genes/GCA_903854225.1/GCA_903854225.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903854225.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:38 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903902845.1.fa
    output: test_extension/extension/genes/GCA_903902845.1/GCA_903902845.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903902845.1/GCA_903902845.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903902845.1/GCA_903902845.1.prodigal_fetchMGs.command
    jobid: 28
    benchmark: test_extension/extension/genes/GCA_903902845.1/GCA_903902845.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903902845.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:39 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903936775.1.fa
    output: test_extension/extension/genes/GCA_903936775.1/GCA_903936775.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903936775.1/GCA_903936775.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903936775.1/GCA_903936775.1.prodigal_fetchMGs.command
    jobid: 43
    benchmark: test_extension/extension/genes/GCA_903936775.1/GCA_903936775.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903936775.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:39 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903839445.1.fa
    output: test_extension/extension/genes/GCA_903839445.1/GCA_903839445.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903839445.1/GCA_903839445.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903839445.1/GCA_903839445.1.prodigal_fetchMGs.command
    jobid: 6
    benchmark: test_extension/extension/genes/GCA_903839445.1/GCA_903839445.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903839445.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:39 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903883715.1.fa
    output: test_extension/extension/genes/GCA_903883715.1/GCA_903883715.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903883715.1/GCA_903883715.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903883715.1/GCA_903883715.1.prodigal_fetchMGs.command
    jobid: 21
    benchmark: test_extension/extension/genes/GCA_903883715.1/GCA_903883715.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903883715.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:39 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903918315.1.fa
    output: test_extension/extension/genes/GCA_903918315.1/GCA_903918315.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903918315.1/GCA_903918315.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903918315.1/GCA_903918315.1.prodigal_fetchMGs.command
    jobid: 36
    benchmark: test_extension/extension/genes/GCA_903918315.1/GCA_903918315.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903918315.1
    resources: tmpdir=/tmp


[Fri Aug 12 11:20:39 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903937725.1.fa
    output: test_extension/extension/genes/GCA_903937725.1/GCA_903937725.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903937725.1/GCA_903937725.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903937725.1/GCA_903937725.1.prodigal_fetchMGs.command
    jobid: 44
    benchmark: test_extension/extension/genes/GCA_903937725.1/GCA_903937725.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903937725.1
    resources: tmpdir=/tmp

Touching output file test_extension/extension/genes/GCA_903854225.1/GCA_903854225.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:30 2022]
Finished job 13.
1 of 51 steps (2%) done
Select jobs to execute...

[Fri Aug 12 11:21:30 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903934285.1.fa
    output: test_extension/extension/genes/GCA_903934285.1/GCA_903934285.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903934285.1/GCA_903934285.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903934285.1/GCA_903934285.1.prodigal_fetchMGs.command
    jobid: 41
    benchmark: test_extension/extension/genes/GCA_903934285.1/GCA_903934285.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903934285.1
    resources: tmpdir=/tmp

Touching output file test_extension/extension/genes/GCA_903853495.1/GCA_903853495.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:34 2022]
Finished job 12.
2 of 51 steps (4%) done
Select jobs to execute...

[Fri Aug 12 11:21:34 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903841135.1.fa
    output: test_extension/extension/genes/GCA_903841135.1/GCA_903841135.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903841135.1/GCA_903841135.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903841135.1/GCA_903841135.1.prodigal_fetchMGs.command
    jobid: 7
    benchmark: test_extension/extension/genes/GCA_903841135.1/GCA_903841135.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903841135.1
    resources: tmpdir=/tmp

Touching output file test_extension/extension/genes/GCA_903888055.1/GCA_903888055.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:45 2022]
Finished job 23.
3 of 51 steps (6%) done
Select jobs to execute...

[Fri Aug 12 11:21:45 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903924625.1.fa
    output: test_extension/extension/genes/GCA_903924625.1/GCA_903924625.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903924625.1/GCA_903924625.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903924625.1/GCA_903924625.1.prodigal_fetchMGs.command
    jobid: 39
    benchmark: test_extension/extension/genes/GCA_903924625.1/GCA_903924625.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903924625.1
    resources: tmpdir=/tmp

Touching output file test_extension/extension/genes/GCA_903826635.1/GCA_903826635.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:47 2022]
Finished job 4.
4 of 51 steps (8%) done
Select jobs to execute...

[Fri Aug 12 11:21:47 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903953565.1.fa
    output: test_extension/extension/genes/GCA_903953565.1/GCA_903953565.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903953565.1/GCA_903953565.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903953565.1/GCA_903953565.1.prodigal_fetchMGs.command
    jobid: 48
    benchmark: test_extension/extension/genes/GCA_903953565.1/GCA_903953565.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903953565.1
    resources: tmpdir=/tmp

Touching output file test_extension/extension/genes/GCA_903939965.1/GCA_903939965.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:50 2022]
Finished job 45.
5 of 51 steps (10%) done
Select jobs to execute...

[Fri Aug 12 11:21:50 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903824795.1.fa
    output: test_extension/extension/genes/GCA_903824795.1/GCA_903824795.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903824795.1/GCA_903824795.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903824795.1/GCA_903824795.1.prodigal_fetchMGs.command
    jobid: 3
    benchmark: test_extension/extension/genes/GCA_903824795.1/GCA_903824795.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903824795.1
    resources: tmpdir=/tmp

Touching output file test_extension/extension/genes/GCA_903835685.1/GCA_903835685.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:50 2022]
Finished job 5.
6 of 51 steps (12%) done
Select jobs to execute...

[Fri Aug 12 11:21:50 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903845665.1.fa
    output: test_extension/extension/genes/GCA_903845665.1/GCA_903845665.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903845665.1/GCA_903845665.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903845665.1/GCA_903845665.1.prodigal_fetchMGs.command
    jobid: 10
    benchmark: test_extension/extension/genes/GCA_903845665.1/GCA_903845665.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903845665.1
    resources: tmpdir=/tmp

Touching output file test_extension/extension/genes/GCA_903945225.1/GCA_903945225.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:50 2022]
Finished job 46.
7 of 51 steps (14%) done
Select jobs to execute...

[Fri Aug 12 11:21:50 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903895255.1.fa
    output: test_extension/extension/genes/GCA_903895255.1/GCA_903895255.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903895255.1/GCA_903895255.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903895255.1/GCA_903895255.1.prodigal_fetchMGs.command
    jobid: 25
    benchmark: test_extension/extension/genes/GCA_903895255.1/GCA_903895255.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903895255.1
    resources: tmpdir=/tmp

Touching output file test_extension/extension/genes/GCA_903935095.1/GCA_903935095.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:50 2022]
Finished job 42.
8 of 51 steps (16%) done
Select jobs to execute...

[Fri Aug 12 11:21:51 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903928895.1.fa
    output: test_extension/extension/genes/GCA_903928895.1/GCA_903928895.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903928895.1/GCA_903928895.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903928895.1/GCA_903928895.1.prodigal_fetchMGs.command
    jobid: 40
    benchmark: test_extension/extension/genes/GCA_903928895.1/GCA_903928895.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903928895.1
    resources: tmpdir=/tmp

Touching output file test_extension/extension/genes/GCA_903879295.1/GCA_903879295.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:51 2022]
Finished job 20.
9 of 51 steps (18%) done
Select jobs to execute...

[Fri Aug 12 11:21:51 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903871585.1.fa
    output: test_extension/extension/genes/GCA_903871585.1/GCA_903871585.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903871585.1/GCA_903871585.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903871585.1/GCA_903871585.1.prodigal_fetchMGs.command
    jobid: 18
    benchmark: test_extension/extension/genes/GCA_903871585.1/GCA_903871585.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903871585.1
    resources: tmpdir=/tmp

Touching output file test_extension/extension/genes/GCA_903959825.1/GCA_903959825.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:51 2022]
Finished job 49.
10 of 51 steps (20%) done
Select jobs to execute...

[Fri Aug 12 11:21:51 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903914905.1.fa
    output: test_extension/extension/genes/GCA_903914905.1/GCA_903914905.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903914905.1/GCA_903914905.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903914905.1/GCA_903914905.1.prodigal_fetchMGs.command
    jobid: 33
    benchmark: test_extension/extension/genes/GCA_903914905.1/GCA_903914905.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903914905.1
    resources: tmpdir=/tmp

Touching output file test_extension/extension/genes/GCA_903824635.1/GCA_903824635.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:52 2022]
Finished job 2.
11 of 51 steps (22%) done
Select jobs to execute...

[Fri Aug 12 11:21:52 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903913555.1.fa
    output: test_extension/extension/genes/GCA_903913555.1/GCA_903913555.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903913555.1/GCA_903913555.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903913555.1/GCA_903913555.1.prodigal_fetchMGs.command
    jobid: 32
    benchmark: test_extension/extension/genes/GCA_903913555.1/GCA_903913555.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903913555.1
    resources: tmpdir=/tmp

Touching output file test_extension/extension/genes/GCA_903883715.1/GCA_903883715.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:52 2022]
Finished job 21.
12 of 51 steps (24%) done
Select jobs to execute...

[Fri Aug 12 11:21:52 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903869725.1.fa
    output: test_extension/extension/genes/GCA_903869725.1/GCA_903869725.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903869725.1/GCA_903869725.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903869725.1/GCA_903869725.1.prodigal_fetchMGs.command
    jobid: 17
    benchmark: test_extension/extension/genes/GCA_903869725.1/GCA_903869725.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903869725.1
    resources: tmpdir=/tmp

Touching output file test_extension/extension/genes/GCA_903900705.1/GCA_903900705.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:53 2022]
Finished job 27.
13 of 51 steps (25%) done
Select jobs to execute...

[Fri Aug 12 11:21:53 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903824045.1.fa
    output: test_extension/extension/genes/GCA_903824045.1/GCA_903824045.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903824045.1/GCA_903824045.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903824045.1/GCA_903824045.1.prodigal_fetchMGs.command
    jobid: 1
    benchmark: test_extension/extension/genes/GCA_903824045.1/GCA_903824045.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903824045.1
    resources: tmpdir=/tmp

Touching output file test_extension/extension/genes/GCA_903842315.1/GCA_903842315.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:53 2022]
Finished job 8.
14 of 51 steps (27%) done
Select jobs to execute...

[Fri Aug 12 11:21:53 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903907865.1.fa
    output: test_extension/extension/genes/GCA_903907865.1/GCA_903907865.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903907865.1/GCA_903907865.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903907865.1/GCA_903907865.1.prodigal_fetchMGs.command
    jobid: 30
    benchmark: test_extension/extension/genes/GCA_903907865.1/GCA_903907865.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903907865.1
    resources: tmpdir=/tmp

Touching output file test_extension/extension/genes/GCA_903937725.1/GCA_903937725.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:54 2022]
Finished job 44.
15 of 51 steps (29%) done
Select jobs to execute...

[Fri Aug 12 11:21:54 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903853455.1.fa
    output: test_extension/extension/genes/GCA_903853455.1/GCA_903853455.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903853455.1/GCA_903853455.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903853455.1/GCA_903853455.1.prodigal_fetchMGs.command
    jobid: 11
    benchmark: test_extension/extension/genes/GCA_903853455.1/GCA_903853455.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903853455.1
    resources: tmpdir=/tmp

Touching output file test_extension/extension/genes/GCA_903857725.1/GCA_903857725.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:54 2022]
Finished job 15.
16 of 51 steps (31%) done
Select jobs to execute...

[Fri Aug 12 11:21:54 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903899515.1.fa
    output: test_extension/extension/genes/GCA_903899515.1/GCA_903899515.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903899515.1/GCA_903899515.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903899515.1/GCA_903899515.1.prodigal_fetchMGs.command
    jobid: 26
    benchmark: test_extension/extension/genes/GCA_903899515.1/GCA_903899515.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903899515.1
    resources: tmpdir=/tmp

Touching output file test_extension/extension/genes/GCA_903902845.1/GCA_903902845.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:54 2022]
Finished job 28.
17 of 51 steps (33%) done
Select jobs to execute...

[Fri Aug 12 11:21:54 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903884585.1.fa
    output: test_extension/extension/genes/GCA_903884585.1/GCA_903884585.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903884585.1/GCA_903884585.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903884585.1/GCA_903884585.1.prodigal_fetchMGs.command
    jobid: 22
    benchmark: test_extension/extension/genes/GCA_903884585.1/GCA_903884585.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903884585.1
    resources: tmpdir=/tmp

Touching output file test_extension/extension/genes/GCA_903918315.1/GCA_903918315.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:54 2022]
Finished job 36.
18 of 51 steps (35%) done
Select jobs to execute...

[Fri Aug 12 11:21:54 2022]
rule prodigalfetchMgs:
    input: test_extension/genomes/GCA_903894295.1.fa
    output: test_extension/extension/genes/GCA_903894295.1/GCA_903894295.1.prodigal_fetchMGs.done
    log: test_extension/extension/genes/GCA_903894295.1/GCA_903894295.1.prodigal_fetchMGs.log, test_extension/extension/genes/GCA_903894295.1/GCA_903894295.1.prodigal_fetchMGs.command
    jobid: 24
    benchmark: test_extension/extension/genes/GCA_903894295.1/GCA_903894295.1.prodigal_fetchMGs.benchmark
    wildcards: sample=GCA_903894295.1
    resources: tmpdir=/tmp

Touching output file test_extension/extension/genes/GCA_903907955.1/GCA_903907955.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:54 2022]
Finished job 31.
19 of 51 steps (37%) done
Touching output file test_extension/extension/genes/GCA_903869175.1/GCA_903869175.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:55 2022]
Finished job 16.
20 of 51 steps (39%) done
Touching output file test_extension/extension/genes/GCA_903843765.1/GCA_903843765.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:55 2022]
Finished job 9.
21 of 51 steps (41%) done
Touching output file test_extension/extension/genes/GCA_903920605.1/GCA_903920605.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:56 2022]
Finished job 37.
22 of 51 steps (43%) done
Touching output file test_extension/extension/genes/GCA_903875175.1/GCA_903875175.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:56 2022]
Finished job 19.
23 of 51 steps (45%) done
Touching output file test_extension/extension/genes/GCA_903915725.1/GCA_903915725.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:56 2022]
Finished job 34.
24 of 51 steps (47%) done
Touching output file test_extension/extension/genes/GCA_903945495.1/GCA_903945495.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:56 2022]
Finished job 47.
25 of 51 steps (49%) done
Touching output file test_extension/extension/genes/GCA_903902895.1/GCA_903902895.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:57 2022]
Finished job 29.
26 of 51 steps (51%) done
Touching output file test_extension/extension/genes/GCA_903839445.1/GCA_903839445.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:57 2022]
Finished job 6.
27 of 51 steps (53%) done
Touching output file test_extension/extension/genes/GCA_903920995.1/GCA_903920995.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:57 2022]
Finished job 38.
28 of 51 steps (55%) done
Touching output file test_extension/extension/genes/GCA_903961455.1/GCA_903961455.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:21:58 2022]
Finished job 50.
29 of 51 steps (57%) done
Touching output file test_extension/extension/genes/GCA_903915855.1/GCA_903915855.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:22:06 2022]
Finished job 35.
30 of 51 steps (59%) done
Touching output file test_extension/extension/genes/GCA_903936775.1/GCA_903936775.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:22:07 2022]
Finished job 43.
31 of 51 steps (61%) done
Touching output file test_extension/extension/genes/GCA_903854255.1/GCA_903854255.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:22:20 2022]
Finished job 14.
32 of 51 steps (63%) done
Touching output file test_extension/extension/genes/GCA_903841135.1/GCA_903841135.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:22:45 2022]
Finished job 7.
33 of 51 steps (65%) done
Touching output file test_extension/extension/genes/GCA_903869725.1/GCA_903869725.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:22:48 2022]
Finished job 17.
34 of 51 steps (67%) done
Touching output file test_extension/extension/genes/GCA_903894295.1/GCA_903894295.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:22:55 2022]
Finished job 24.
35 of 51 steps (69%) done
Touching output file test_extension/extension/genes/GCA_903914905.1/GCA_903914905.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:22:56 2022]
Finished job 33.
36 of 51 steps (71%) done
Touching output file test_extension/extension/genes/GCA_903824045.1/GCA_903824045.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:22:57 2022]
Finished job 1.
37 of 51 steps (73%) done
Touching output file test_extension/extension/genes/GCA_903884585.1/GCA_903884585.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:22:58 2022]
Finished job 22.
38 of 51 steps (75%) done
Touching output file test_extension/extension/genes/GCA_903953565.1/GCA_903953565.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:22:59 2022]
Finished job 48.
39 of 51 steps (76%) done
Touching output file test_extension/extension/genes/GCA_903913555.1/GCA_903913555.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:23:00 2022]
Finished job 32.
40 of 51 steps (78%) done
Touching output file test_extension/extension/genes/GCA_903934285.1/GCA_903934285.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:23:01 2022]
Finished job 41.
41 of 51 steps (80%) done
Touching output file test_extension/extension/genes/GCA_903845665.1/GCA_903845665.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:23:01 2022]
Finished job 10.
42 of 51 steps (82%) done
Touching output file test_extension/extension/genes/GCA_903928895.1/GCA_903928895.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:23:01 2022]
Finished job 40.
43 of 51 steps (84%) done
Touching output file test_extension/extension/genes/GCA_903871585.1/GCA_903871585.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:23:02 2022]
Finished job 18.
44 of 51 steps (86%) done
Touching output file test_extension/extension/genes/GCA_903824795.1/GCA_903824795.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:23:03 2022]
Finished job 3.
45 of 51 steps (88%) done
Touching output file test_extension/extension/genes/GCA_903853455.1/GCA_903853455.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:23:03 2022]
Finished job 11.
46 of 51 steps (90%) done
Touching output file test_extension/extension/genes/GCA_903899515.1/GCA_903899515.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:23:06 2022]
Finished job 26.
47 of 51 steps (92%) done
Touching output file test_extension/extension/genes/GCA_903907865.1/GCA_903907865.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:23:08 2022]
Finished job 30.
48 of 51 steps (94%) done
Touching output file test_extension/extension/genes/GCA_903895255.1/GCA_903895255.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:23:15 2022]
Finished job 25.
49 of 51 steps (96%) done
Touching output file test_extension/extension/genes/GCA_903924625.1/GCA_903924625.1.prodigal_fetchMGs.done.
[Fri Aug 12 11:23:18 2022]
Finished job 39.
50 of 51 steps (98%) done
Select jobs to execute...

[Fri Aug 12 11:23:18 2022]
localrule all:
    input: test_extension/extension/genes/GCA_903824045.1/GCA_903824045.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903824635.1/GCA_903824635.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903824795.1/GCA_903824795.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903826635.1/GCA_903826635.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903835685.1/GCA_903835685.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903839445.1/GCA_903839445.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903841135.1/GCA_903841135.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903842315.1/GCA_903842315.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903843765.1/GCA_903843765.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903845665.1/GCA_903845665.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903853455.1/GCA_903853455.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903853495.1/GCA_903853495.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903854225.1/GCA_903854225.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903854255.1/GCA_903854255.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903857725.1/GCA_903857725.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903869175.1/GCA_903869175.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903869725.1/GCA_903869725.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903871585.1/GCA_903871585.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903875175.1/GCA_903875175.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903879295.1/GCA_903879295.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903883715.1/GCA_903883715.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903884585.1/GCA_903884585.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903888055.1/GCA_903888055.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903894295.1/GCA_903894295.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903895255.1/GCA_903895255.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903899515.1/GCA_903899515.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903900705.1/GCA_903900705.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903902845.1/GCA_903902845.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903902895.1/GCA_903902895.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903907865.1/GCA_903907865.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903907955.1/GCA_903907955.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903913555.1/GCA_903913555.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903914905.1/GCA_903914905.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903915725.1/GCA_903915725.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903915855.1/GCA_903915855.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903918315.1/GCA_903918315.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903920605.1/GCA_903920605.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903920995.1/GCA_903920995.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903924625.1/GCA_903924625.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903928895.1/GCA_903928895.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903934285.1/GCA_903934285.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903935095.1/GCA_903935095.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903936775.1/GCA_903936775.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903937725.1/GCA_903937725.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903939965.1/GCA_903939965.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903945225.1/GCA_903945225.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903945495.1/GCA_903945495.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903953565.1/GCA_903953565.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903959825.1/GCA_903959825.1.prodigal_fetchMGs.done, test_extension/extension/genes/GCA_903961455.1/GCA_903961455.1.prodigal_fetchMGs.done
    jobid: 0
    resources: tmpdir=/tmp

[Fri Aug 12 11:23:18 2022]
Finished job 0.
51 of 51 steps (100%) done
Complete log: .snakemake/log/2022-08-12T112037.598847.snakemake.log
2022-08-12 11:23:18,293 INFO: Processing genome assignation
2022-08-12 11:23:18,571 INFO: Processed 50 / 50
2022-08-12 11:23:18,575 INFO: Total genomes: 50
2022-08-12 11:23:18,575 INFO: 	 Assigned to refmotus 			0
2022-08-12 11:23:18,575 INFO: 	 Assigned to metamotus 			0
2022-08-12 11:23:18,575 INFO: 	 Assigned to extmotus 			0
2022-08-12 11:23:18,575 INFO: 	 Genomes not assigned to any mOTU 	50
2022-08-12 11:23:18,575 INFO: 	 Filtered genomes (<6 MGs) 		0
2022-08-12 11:23:18,577 INFO: mOTUs-extender finished successfully

```
</details>

The membership file will show whether a genome is already represented by an existing mOTU of if it represents taxonomic novelty.

```
cat test_extension/extension/mOTUs.membership.tsv

GCA_903824045.1	Novel	10
GCA_903824635.1	Novel	10
GCA_903824795.1	Novel	10
GCA_903826635.1	Novel	10
GCA_903835685.1	Novel	10
GCA_903839445.1	Novel	10
GCA_903841135.1	Novel	10
GCA_903842315.1	Novel	10
GCA_903843765.1	Novel	10
GCA_903845665.1	Novel	10
GCA_903853455.1	Novel	10
GCA_903853495.1	Novel	8
GCA_903854225.1	Novel	9
GCA_903854255.1	Novel	10
...

```

**Note**: This step can, depending on the number of genomes, take from minutes to hours.


### 3. mOTUs extender - createdb

```
python mOTUs-extender/mOTUs-extender/motus-extender.py createdb  -w test_extension/extension/ -a mOTUs-extender/mOTUs-extender/test/genomes.tax -p NEWDB -t 32
```
<details>
<summary>Commandline stdout</summary>
<br>

```

python mOTUs-extender/mOTUs-extender/motus-extender.py createdb  -w test_extension/extension/ -a mOTUs-extender/mOTUs-extender/test/genomes.tax -p NEWDB -t 32
2022-08-12 14:44:05,495 INFO: Starting mOTUs-extender
2022-08-12 14:44:05,502 INFO: Executing: sh mOTUs-extender/mOTUs-extender/extend_mOTUs_generateDB.sh test_extension/extension/genomes.list NEWDB mOTUs-extender/mOTUs-extender/test/genomes.tax test_extension/extension/extension/ mOTUs-extender/mOTUs-extender/ test_extension/extension/temp_db_folder/ 32
+ fileListGenomes=test_extension/extension/genomes.list
+ newDBName=NEWDB
+ taxonomyFile=mOTUs-extender/mOTUs-extender/test/genomes.tax
+ new_database_folder=test_extension/extension/extension/
+ scriptDir=mOTUs-extender/mOTUs-extender/
+ mOTU_folder=test_extension/extension/temp_db_folder/
+ threads=32
+ cutoffsFile=mOTUs-extender/mOTUs-extender//cutoffs_fscore_specIAsRef.csv
+ mOTU_MG_file=mOTUs-extender/mOTUs-extender//10RefMGs.IDs
+ cutoff=99.0
+ fileListGenomes2=test_extension/extension/extension//NEWDB.newGenomes.filteredMissingMGs.txt
+ touch test_extension/extension/extension//NEWDB.newGenomes.filteredMissingMGs.txt
+ rm test_extension/extension/extension//NEWDB.newGenomes.filteredMissingMGs.txt
+ touch test_extension/extension/extension//NEWDB.newGenomes.filteredMissingMGs.txt
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903824045.1/GCA_903824045.1.map ']'
+ echo GCA_903824045.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903824635.1/GCA_903824635.1.map ']'
+ echo GCA_903824635.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903824795.1/GCA_903824795.1.map ']'
+ echo GCA_903824795.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903826635.1/GCA_903826635.1.map ']'
+ echo GCA_903826635.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903835685.1/GCA_903835685.1.map ']'
+ echo GCA_903835685.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903839445.1/GCA_903839445.1.map ']'
+ echo GCA_903839445.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903841135.1/GCA_903841135.1.map ']'
+ echo GCA_903841135.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903842315.1/GCA_903842315.1.map ']'
+ echo GCA_903842315.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903843765.1/GCA_903843765.1.map ']'
+ echo GCA_903843765.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903845665.1/GCA_903845665.1.map ']'
+ echo GCA_903845665.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903853455.1/GCA_903853455.1.map ']'
+ echo GCA_903853455.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903853495.1/GCA_903853495.1.map ']'
+ echo GCA_903853495.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903854225.1/GCA_903854225.1.map ']'
+ echo GCA_903854225.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903854255.1/GCA_903854255.1.map ']'
+ echo GCA_903854255.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903857725.1/GCA_903857725.1.map ']'
+ echo GCA_903857725.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903869175.1/GCA_903869175.1.map ']'
+ echo GCA_903869175.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903869725.1/GCA_903869725.1.map ']'
+ echo GCA_903869725.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903871585.1/GCA_903871585.1.map ']'
+ echo GCA_903871585.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903875175.1/GCA_903875175.1.map ']'
+ echo GCA_903875175.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903879295.1/GCA_903879295.1.map ']'
+ echo GCA_903879295.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903883715.1/GCA_903883715.1.map ']'
+ echo GCA_903883715.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903884585.1/GCA_903884585.1.map ']'
+ echo GCA_903884585.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903888055.1/GCA_903888055.1.map ']'
+ echo GCA_903888055.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903894295.1/GCA_903894295.1.map ']'
+ echo GCA_903894295.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903895255.1/GCA_903895255.1.map ']'
+ echo GCA_903895255.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903899515.1/GCA_903899515.1.map ']'
+ echo GCA_903899515.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903900705.1/GCA_903900705.1.map ']'
+ echo GCA_903900705.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903902845.1/GCA_903902845.1.map ']'
+ echo GCA_903902845.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903902895.1/GCA_903902895.1.map ']'
+ echo GCA_903902895.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903907865.1/GCA_903907865.1.map ']'
+ echo GCA_903907865.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903907955.1/GCA_903907955.1.map ']'
+ echo GCA_903907955.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903913555.1/GCA_903913555.1.map ']'
+ echo GCA_903913555.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903914905.1/GCA_903914905.1.map ']'
+ echo GCA_903914905.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903915725.1/GCA_903915725.1.map ']'
+ echo GCA_903915725.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903915855.1/GCA_903915855.1.map ']'
+ echo GCA_903915855.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903918315.1/GCA_903918315.1.map ']'
+ echo GCA_903918315.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903920605.1/GCA_903920605.1.map ']'
+ echo GCA_903920605.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903920995.1/GCA_903920995.1.map ']'
+ echo GCA_903920995.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903924625.1/GCA_903924625.1.map ']'
+ echo GCA_903924625.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903928895.1/GCA_903928895.1.map ']'
+ echo GCA_903928895.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903934285.1/GCA_903934285.1.map ']'
+ echo GCA_903934285.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903935095.1/GCA_903935095.1.map ']'
+ echo GCA_903935095.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903936775.1/GCA_903936775.1.map ']'
+ echo GCA_903936775.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903937725.1/GCA_903937725.1.map ']'
+ echo GCA_903937725.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903939965.1/GCA_903939965.1.map ']'
+ echo GCA_903939965.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903945225.1/GCA_903945225.1.map ']'
+ echo GCA_903945225.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903945495.1/GCA_903945495.1.map ']'
+ echo GCA_903945495.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903953565.1/GCA_903953565.1.map ']'
+ echo GCA_903953565.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903959825.1/GCA_903959825.1.map ']'
+ echo GCA_903959825.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//dbs/GCA_903961455.1/GCA_903961455.1.map ']'
+ echo GCA_903961455.1
+ read genome_ID
+ '[' -s test_extension/extension/extension//NEWDB.newGenomes.filteredMissingMGs.txt ']'
+ echo 'Adding genomes to database.'
Adding genomes to database.
+ mkdir -p test_extension/extension/extension//NEWDB/vsearch
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824045.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824635.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824795.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903826635.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903835685.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903839445.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903841135.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903842315.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903843765.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903845665.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853455.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853495.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854225.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854255.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903857725.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869175.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869725.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903871585.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903875175.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903879295.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903883715.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903884585.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903888055.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903894295.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903895255.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903899515.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903900705.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902845.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902895.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907865.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907955.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903913555.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903914905.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915725.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915855.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903918315.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920605.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920995.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903924625.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903928895.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903934285.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903935095.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903936775.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903937725.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903939965.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945225.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945495.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903953565.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903959825.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903961455.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ python mOTUs-extender/mOTUs-extender//extend_mOTUs_filterGenomesByDistance.py --cutoff 99.0 test_extension/extension/extension//NEWDB/vsearch/combined.normalized.distances_vs_db.m8 test_extension/extension/extension//NEWDB.newGenomes.filteredMissingMGs.txt
+ '[' -s test_extension/extension/extension//NEWDB/genomes.filtered.list ']'
+ echo 'Adding genomes listed in test_extension/extension/extension//NEWDB/genomes.filtered.list.'
Adding genomes listed in test_extension/extension/extension//NEWDB/genomes.filtered.list.
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903843765.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903913555.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903842315.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903835685.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903841135.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903884585.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853455.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903937725.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903883715.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903914905.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903894295.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903936775.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854225.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903953565.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915855.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903857725.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824635.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903875175.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902845.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853495.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854255.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903934285.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903839445.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945225.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903918315.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903826635.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915725.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903924625.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903879295.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920995.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903939965.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824045.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903845665.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907865.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903900705.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869725.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903895255.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903935095.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903888055.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903871585.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824795.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907955.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903959825.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920605.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903928895.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945495.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869175.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902895.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903899515.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903961455.1/vsearch/combined.normalized.distances_vs_db.m8
+ read genome_ID
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903843765.1/GCA_903843765.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903913555.1/GCA_903913555.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903842315.1/GCA_903842315.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903835685.1/GCA_903835685.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903841135.1/GCA_903841135.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903884585.1/GCA_903884585.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853455.1/GCA_903853455.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903937725.1/GCA_903937725.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903883715.1/GCA_903883715.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903914905.1/GCA_903914905.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903894295.1/GCA_903894295.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903936775.1/GCA_903936775.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854225.1/GCA_903854225.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903953565.1/GCA_903953565.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915855.1/GCA_903915855.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903857725.1/GCA_903857725.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824635.1/GCA_903824635.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903875175.1/GCA_903875175.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902845.1/GCA_903902845.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853495.1/GCA_903853495.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854255.1/GCA_903854255.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903934285.1/GCA_903934285.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903839445.1/GCA_903839445.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945225.1/GCA_903945225.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903918315.1/GCA_903918315.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903826635.1/GCA_903826635.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915725.1/GCA_903915725.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903924625.1/GCA_903924625.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903879295.1/GCA_903879295.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920995.1/GCA_903920995.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903939965.1/GCA_903939965.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824045.1/GCA_903824045.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903845665.1/GCA_903845665.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907865.1/GCA_903907865.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903900705.1/GCA_903900705.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869725.1/GCA_903869725.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903895255.1/GCA_903895255.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903935095.1/GCA_903935095.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903888055.1/GCA_903888055.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903871585.1/GCA_903871585.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824795.1/GCA_903824795.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907955.1/GCA_903907955.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903959825.1/GCA_903959825.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920605.1/GCA_903920605.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903928895.1/GCA_903928895.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945495.1/GCA_903945495.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869175.1/GCA_903869175.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902895.1/GCA_903902895.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903899515.1/GCA_903899515.1.map
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903961455.1/GCA_903961455.1.map
+ read genome_ID
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903843765.1/GCA_903843765.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903913555.1/GCA_903913555.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903842315.1/GCA_903842315.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903835685.1/GCA_903835685.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903841135.1/GCA_903841135.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903884585.1/GCA_903884585.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853455.1/GCA_903853455.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903937725.1/GCA_903937725.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903883715.1/GCA_903883715.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903914905.1/GCA_903914905.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903894295.1/GCA_903894295.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903936775.1/GCA_903936775.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854225.1/GCA_903854225.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903953565.1/GCA_903953565.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915855.1/GCA_903915855.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903857725.1/GCA_903857725.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824635.1/GCA_903824635.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903875175.1/GCA_903875175.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902845.1/GCA_903902845.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853495.1/GCA_903853495.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854255.1/GCA_903854255.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903934285.1/GCA_903934285.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903839445.1/GCA_903839445.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945225.1/GCA_903945225.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903918315.1/GCA_903918315.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903826635.1/GCA_903826635.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915725.1/GCA_903915725.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903924625.1/GCA_903924625.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903879295.1/GCA_903879295.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920995.1/GCA_903920995.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903939965.1/GCA_903939965.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824045.1/GCA_903824045.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903845665.1/GCA_903845665.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907865.1/GCA_903907865.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903900705.1/GCA_903900705.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869725.1/GCA_903869725.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903895255.1/GCA_903895255.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903935095.1/GCA_903935095.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903888055.1/GCA_903888055.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903871585.1/GCA_903871585.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824795.1/GCA_903824795.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907955.1/GCA_903907955.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903959825.1/GCA_903959825.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920605.1/GCA_903920605.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903928895.1/GCA_903928895.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945495.1/GCA_903945495.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869175.1/GCA_903869175.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902895.1/GCA_903902895.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903899515.1/GCA_903899515.1.map2genome
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903961455.1/GCA_903961455.1.map2genome
+ read genome_ID
+ cut -f2 test_extension/extension/extension//NEWDB/NEWDB.new.map2genome
+ sort -u
+ python mOTUs-extender/mOTUs-extender//map2genome.py test_extension/extension/extension//NEWDB/NEWDB.new.map2genome test_extension/extension/extension//NEWDB/NEWDB.new.map2genome_v5
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903843765.1/GCA_903843765.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903913555.1/GCA_903913555.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903842315.1/GCA_903842315.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903835685.1/GCA_903835685.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903841135.1/GCA_903841135.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903884585.1/GCA_903884585.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853455.1/GCA_903853455.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903937725.1/GCA_903937725.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903883715.1/GCA_903883715.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903914905.1/GCA_903914905.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903894295.1/GCA_903894295.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903936775.1/GCA_903936775.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854225.1/GCA_903854225.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903953565.1/GCA_903953565.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915855.1/GCA_903915855.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903857725.1/GCA_903857725.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824635.1/GCA_903824635.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903875175.1/GCA_903875175.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902845.1/GCA_903902845.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853495.1/GCA_903853495.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854255.1/GCA_903854255.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903934285.1/GCA_903934285.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903839445.1/GCA_903839445.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945225.1/GCA_903945225.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903918315.1/GCA_903918315.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903826635.1/GCA_903826635.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915725.1/GCA_903915725.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903924625.1/GCA_903924625.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903879295.1/GCA_903879295.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920995.1/GCA_903920995.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903939965.1/GCA_903939965.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824045.1/GCA_903824045.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903845665.1/GCA_903845665.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907865.1/GCA_903907865.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903900705.1/GCA_903900705.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869725.1/GCA_903869725.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903895255.1/GCA_903895255.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903935095.1/GCA_903935095.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903888055.1/GCA_903888055.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903871585.1/GCA_903871585.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824795.1/GCA_903824795.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907955.1/GCA_903907955.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903959825.1/GCA_903959825.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920605.1/GCA_903920605.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903928895.1/GCA_903928895.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945495.1/GCA_903945495.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869175.1/GCA_903869175.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902895.1/GCA_903902895.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903899515.1/GCA_903899515.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903961455.1/GCA_903961455.1.genes.len
+ read genome_ID
+ cat test_extension/extension/extension//NEWDB/NEWDB.new.len test_extension/extension/temp_db_folder//db_mOTU/db_mOTU_genes_length_NR
+ mkdir -p test_extension/extension/extension//NEWDB/sequences/
+ read line
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903843765.1/GCA_903843765.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903913555.1/GCA_903913555.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903842315.1/GCA_903842315.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903835685.1/GCA_903835685.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903841135.1/GCA_903841135.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903884585.1/GCA_903884585.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853455.1/GCA_903853455.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903937725.1/GCA_903937725.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903883715.1/GCA_903883715.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903914905.1/GCA_903914905.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903894295.1/GCA_903894295.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903936775.1/GCA_903936775.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854225.1/GCA_903854225.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903953565.1/GCA_903953565.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915855.1/GCA_903915855.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903857725.1/GCA_903857725.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824635.1/GCA_903824635.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903875175.1/GCA_903875175.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902845.1/GCA_903902845.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853495.1/GCA_903853495.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854255.1/GCA_903854255.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903934285.1/GCA_903934285.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903839445.1/GCA_903839445.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945225.1/GCA_903945225.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903918315.1/GCA_903918315.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903826635.1/GCA_903826635.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915725.1/GCA_903915725.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903924625.1/GCA_903924625.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903879295.1/GCA_903879295.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920995.1/GCA_903920995.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903939965.1/GCA_903939965.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824045.1/GCA_903824045.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903845665.1/GCA_903845665.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907865.1/GCA_903907865.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903900705.1/GCA_903900705.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869725.1/GCA_903869725.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903895255.1/GCA_903895255.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903935095.1/GCA_903935095.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903888055.1/GCA_903888055.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903871585.1/GCA_903871585.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824795.1/GCA_903824795.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907955.1/GCA_903907955.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903959825.1/GCA_903959825.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920605.1/GCA_903920605.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903928895.1/GCA_903928895.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945495.1/GCA_903945495.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869175.1/GCA_903869175.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902895.1/GCA_903902895.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903899515.1/GCA_903899515.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903961455.1/GCA_903961455.1_markerGenes2/COG0012.notab.fna
+ read genome_ID
+ read line
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903843765.1/GCA_903843765.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903913555.1/GCA_903913555.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903842315.1/GCA_903842315.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903835685.1/GCA_903835685.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903841135.1/GCA_903841135.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903884585.1/GCA_903884585.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853455.1/GCA_903853455.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903937725.1/GCA_903937725.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903883715.1/GCA_903883715.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903914905.1/GCA_903914905.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903894295.1/GCA_903894295.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903936775.1/GCA_903936775.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854225.1/GCA_903854225.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903953565.1/GCA_903953565.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915855.1/GCA_903915855.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903857725.1/GCA_903857725.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824635.1/GCA_903824635.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903875175.1/GCA_903875175.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902845.1/GCA_903902845.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853495.1/GCA_903853495.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854255.1/GCA_903854255.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903934285.1/GCA_903934285.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903839445.1/GCA_903839445.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945225.1/GCA_903945225.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903918315.1/GCA_903918315.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903826635.1/GCA_903826635.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915725.1/GCA_903915725.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903924625.1/GCA_903924625.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903879295.1/GCA_903879295.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920995.1/GCA_903920995.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903939965.1/GCA_903939965.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824045.1/GCA_903824045.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903845665.1/GCA_903845665.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907865.1/GCA_903907865.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903900705.1/GCA_903900705.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869725.1/GCA_903869725.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903895255.1/GCA_903895255.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903935095.1/GCA_903935095.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903888055.1/GCA_903888055.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903871585.1/GCA_903871585.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824795.1/GCA_903824795.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907955.1/GCA_903907955.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903959825.1/GCA_903959825.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920605.1/GCA_903920605.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903928895.1/GCA_903928895.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945495.1/GCA_903945495.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869175.1/GCA_903869175.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902895.1/GCA_903902895.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903899515.1/GCA_903899515.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903961455.1/GCA_903961455.1_markerGenes2/COG0016.notab.fna
+ read genome_ID
+ read line
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903843765.1/GCA_903843765.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903913555.1/GCA_903913555.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903842315.1/GCA_903842315.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903835685.1/GCA_903835685.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903841135.1/GCA_903841135.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903884585.1/GCA_903884585.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853455.1/GCA_903853455.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903937725.1/GCA_903937725.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903883715.1/GCA_903883715.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903914905.1/GCA_903914905.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903894295.1/GCA_903894295.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903936775.1/GCA_903936775.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854225.1/GCA_903854225.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903953565.1/GCA_903953565.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915855.1/GCA_903915855.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903857725.1/GCA_903857725.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824635.1/GCA_903824635.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903875175.1/GCA_903875175.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902845.1/GCA_903902845.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853495.1/GCA_903853495.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854255.1/GCA_903854255.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903934285.1/GCA_903934285.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903839445.1/GCA_903839445.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945225.1/GCA_903945225.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903918315.1/GCA_903918315.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903826635.1/GCA_903826635.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915725.1/GCA_903915725.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903924625.1/GCA_903924625.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903879295.1/GCA_903879295.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920995.1/GCA_903920995.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903939965.1/GCA_903939965.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824045.1/GCA_903824045.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903845665.1/GCA_903845665.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907865.1/GCA_903907865.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903900705.1/GCA_903900705.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869725.1/GCA_903869725.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903895255.1/GCA_903895255.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903935095.1/GCA_903935095.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903888055.1/GCA_903888055.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903871585.1/GCA_903871585.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824795.1/GCA_903824795.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907955.1/GCA_903907955.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903959825.1/GCA_903959825.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920605.1/GCA_903920605.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903928895.1/GCA_903928895.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945495.1/GCA_903945495.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869175.1/GCA_903869175.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902895.1/GCA_903902895.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903899515.1/GCA_903899515.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903961455.1/GCA_903961455.1_markerGenes2/COG0018.notab.fna
+ read genome_ID
+ read line
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903843765.1/GCA_903843765.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903913555.1/GCA_903913555.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903842315.1/GCA_903842315.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903835685.1/GCA_903835685.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903841135.1/GCA_903841135.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903884585.1/GCA_903884585.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853455.1/GCA_903853455.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903937725.1/GCA_903937725.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903883715.1/GCA_903883715.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903914905.1/GCA_903914905.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903894295.1/GCA_903894295.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903936775.1/GCA_903936775.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854225.1/GCA_903854225.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903953565.1/GCA_903953565.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915855.1/GCA_903915855.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903857725.1/GCA_903857725.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824635.1/GCA_903824635.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903875175.1/GCA_903875175.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902845.1/GCA_903902845.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853495.1/GCA_903853495.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854255.1/GCA_903854255.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903934285.1/GCA_903934285.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903839445.1/GCA_903839445.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945225.1/GCA_903945225.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903918315.1/GCA_903918315.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903826635.1/GCA_903826635.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915725.1/GCA_903915725.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903924625.1/GCA_903924625.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903879295.1/GCA_903879295.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920995.1/GCA_903920995.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903939965.1/GCA_903939965.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824045.1/GCA_903824045.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903845665.1/GCA_903845665.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907865.1/GCA_903907865.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903900705.1/GCA_903900705.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869725.1/GCA_903869725.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903895255.1/GCA_903895255.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903935095.1/GCA_903935095.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903888055.1/GCA_903888055.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903871585.1/GCA_903871585.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824795.1/GCA_903824795.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907955.1/GCA_903907955.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903959825.1/GCA_903959825.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920605.1/GCA_903920605.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903928895.1/GCA_903928895.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945495.1/GCA_903945495.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869175.1/GCA_903869175.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902895.1/GCA_903902895.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903899515.1/GCA_903899515.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903961455.1/GCA_903961455.1_markerGenes2/COG0172.notab.fna
+ read genome_ID
+ read line
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903843765.1/GCA_903843765.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903913555.1/GCA_903913555.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903842315.1/GCA_903842315.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903835685.1/GCA_903835685.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903841135.1/GCA_903841135.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903884585.1/GCA_903884585.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853455.1/GCA_903853455.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903937725.1/GCA_903937725.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903883715.1/GCA_903883715.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903914905.1/GCA_903914905.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903894295.1/GCA_903894295.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903936775.1/GCA_903936775.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854225.1/GCA_903854225.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903953565.1/GCA_903953565.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915855.1/GCA_903915855.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903857725.1/GCA_903857725.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824635.1/GCA_903824635.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903875175.1/GCA_903875175.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902845.1/GCA_903902845.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853495.1/GCA_903853495.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854255.1/GCA_903854255.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903934285.1/GCA_903934285.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903839445.1/GCA_903839445.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945225.1/GCA_903945225.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903918315.1/GCA_903918315.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903826635.1/GCA_903826635.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915725.1/GCA_903915725.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903924625.1/GCA_903924625.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903879295.1/GCA_903879295.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920995.1/GCA_903920995.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903939965.1/GCA_903939965.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824045.1/GCA_903824045.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903845665.1/GCA_903845665.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907865.1/GCA_903907865.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903900705.1/GCA_903900705.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869725.1/GCA_903869725.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903895255.1/GCA_903895255.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903935095.1/GCA_903935095.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903888055.1/GCA_903888055.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903871585.1/GCA_903871585.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824795.1/GCA_903824795.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907955.1/GCA_903907955.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903959825.1/GCA_903959825.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920605.1/GCA_903920605.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903928895.1/GCA_903928895.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945495.1/GCA_903945495.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869175.1/GCA_903869175.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902895.1/GCA_903902895.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903899515.1/GCA_903899515.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903961455.1/GCA_903961455.1_markerGenes2/COG0215.notab.fna
+ read genome_ID
+ read line
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903843765.1/GCA_903843765.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903913555.1/GCA_903913555.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903842315.1/GCA_903842315.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903835685.1/GCA_903835685.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903841135.1/GCA_903841135.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903884585.1/GCA_903884585.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853455.1/GCA_903853455.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903937725.1/GCA_903937725.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903883715.1/GCA_903883715.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903914905.1/GCA_903914905.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903894295.1/GCA_903894295.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903936775.1/GCA_903936775.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854225.1/GCA_903854225.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903953565.1/GCA_903953565.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915855.1/GCA_903915855.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903857725.1/GCA_903857725.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824635.1/GCA_903824635.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903875175.1/GCA_903875175.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902845.1/GCA_903902845.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853495.1/GCA_903853495.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854255.1/GCA_903854255.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903934285.1/GCA_903934285.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903839445.1/GCA_903839445.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945225.1/GCA_903945225.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903918315.1/GCA_903918315.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903826635.1/GCA_903826635.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915725.1/GCA_903915725.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903924625.1/GCA_903924625.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903879295.1/GCA_903879295.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920995.1/GCA_903920995.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903939965.1/GCA_903939965.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824045.1/GCA_903824045.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903845665.1/GCA_903845665.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907865.1/GCA_903907865.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903900705.1/GCA_903900705.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869725.1/GCA_903869725.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903895255.1/GCA_903895255.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903935095.1/GCA_903935095.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903888055.1/GCA_903888055.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903871585.1/GCA_903871585.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824795.1/GCA_903824795.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907955.1/GCA_903907955.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903959825.1/GCA_903959825.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920605.1/GCA_903920605.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903928895.1/GCA_903928895.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945495.1/GCA_903945495.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869175.1/GCA_903869175.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902895.1/GCA_903902895.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903899515.1/GCA_903899515.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903961455.1/GCA_903961455.1_markerGenes2/COG0495.notab.fna
+ read genome_ID
+ read line
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903843765.1/GCA_903843765.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903913555.1/GCA_903913555.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903842315.1/GCA_903842315.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903835685.1/GCA_903835685.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903841135.1/GCA_903841135.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903884585.1/GCA_903884585.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853455.1/GCA_903853455.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903937725.1/GCA_903937725.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903883715.1/GCA_903883715.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903914905.1/GCA_903914905.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903894295.1/GCA_903894295.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903936775.1/GCA_903936775.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854225.1/GCA_903854225.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903953565.1/GCA_903953565.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915855.1/GCA_903915855.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903857725.1/GCA_903857725.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824635.1/GCA_903824635.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903875175.1/GCA_903875175.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902845.1/GCA_903902845.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853495.1/GCA_903853495.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854255.1/GCA_903854255.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903934285.1/GCA_903934285.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903839445.1/GCA_903839445.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945225.1/GCA_903945225.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903918315.1/GCA_903918315.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903826635.1/GCA_903826635.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915725.1/GCA_903915725.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903924625.1/GCA_903924625.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903879295.1/GCA_903879295.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920995.1/GCA_903920995.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903939965.1/GCA_903939965.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824045.1/GCA_903824045.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903845665.1/GCA_903845665.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907865.1/GCA_903907865.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903900705.1/GCA_903900705.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869725.1/GCA_903869725.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903895255.1/GCA_903895255.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903935095.1/GCA_903935095.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903888055.1/GCA_903888055.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903871585.1/GCA_903871585.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824795.1/GCA_903824795.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907955.1/GCA_903907955.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903959825.1/GCA_903959825.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920605.1/GCA_903920605.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903928895.1/GCA_903928895.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945495.1/GCA_903945495.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869175.1/GCA_903869175.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902895.1/GCA_903902895.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903899515.1/GCA_903899515.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903961455.1/GCA_903961455.1_markerGenes2/COG0525.notab.fna
+ read genome_ID
+ read line
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903843765.1/GCA_903843765.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903913555.1/GCA_903913555.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903842315.1/GCA_903842315.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903835685.1/GCA_903835685.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903841135.1/GCA_903841135.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903884585.1/GCA_903884585.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853455.1/GCA_903853455.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903937725.1/GCA_903937725.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903883715.1/GCA_903883715.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903914905.1/GCA_903914905.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903894295.1/GCA_903894295.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903936775.1/GCA_903936775.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854225.1/GCA_903854225.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903953565.1/GCA_903953565.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915855.1/GCA_903915855.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903857725.1/GCA_903857725.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824635.1/GCA_903824635.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903875175.1/GCA_903875175.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902845.1/GCA_903902845.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853495.1/GCA_903853495.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854255.1/GCA_903854255.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903934285.1/GCA_903934285.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903839445.1/GCA_903839445.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945225.1/GCA_903945225.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903918315.1/GCA_903918315.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903826635.1/GCA_903826635.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915725.1/GCA_903915725.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903924625.1/GCA_903924625.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903879295.1/GCA_903879295.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920995.1/GCA_903920995.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903939965.1/GCA_903939965.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824045.1/GCA_903824045.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903845665.1/GCA_903845665.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907865.1/GCA_903907865.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903900705.1/GCA_903900705.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869725.1/GCA_903869725.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903895255.1/GCA_903895255.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903935095.1/GCA_903935095.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903888055.1/GCA_903888055.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903871585.1/GCA_903871585.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824795.1/GCA_903824795.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907955.1/GCA_903907955.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903959825.1/GCA_903959825.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920605.1/GCA_903920605.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903928895.1/GCA_903928895.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945495.1/GCA_903945495.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869175.1/GCA_903869175.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902895.1/GCA_903902895.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903899515.1/GCA_903899515.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903961455.1/GCA_903961455.1_markerGenes2/COG0533.notab.fna
+ read genome_ID
+ read line
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903843765.1/GCA_903843765.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903913555.1/GCA_903913555.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903842315.1/GCA_903842315.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903835685.1/GCA_903835685.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903841135.1/GCA_903841135.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903884585.1/GCA_903884585.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853455.1/GCA_903853455.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903937725.1/GCA_903937725.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903883715.1/GCA_903883715.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903914905.1/GCA_903914905.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903894295.1/GCA_903894295.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903936775.1/GCA_903936775.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854225.1/GCA_903854225.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903953565.1/GCA_903953565.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915855.1/GCA_903915855.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903857725.1/GCA_903857725.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824635.1/GCA_903824635.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903875175.1/GCA_903875175.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902845.1/GCA_903902845.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853495.1/GCA_903853495.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854255.1/GCA_903854255.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903934285.1/GCA_903934285.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903839445.1/GCA_903839445.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945225.1/GCA_903945225.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903918315.1/GCA_903918315.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903826635.1/GCA_903826635.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915725.1/GCA_903915725.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903924625.1/GCA_903924625.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903879295.1/GCA_903879295.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920995.1/GCA_903920995.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903939965.1/GCA_903939965.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824045.1/GCA_903824045.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903845665.1/GCA_903845665.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907865.1/GCA_903907865.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903900705.1/GCA_903900705.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869725.1/GCA_903869725.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903895255.1/GCA_903895255.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903935095.1/GCA_903935095.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903888055.1/GCA_903888055.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903871585.1/GCA_903871585.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824795.1/GCA_903824795.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907955.1/GCA_903907955.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903959825.1/GCA_903959825.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920605.1/GCA_903920605.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903928895.1/GCA_903928895.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945495.1/GCA_903945495.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869175.1/GCA_903869175.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902895.1/GCA_903902895.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903899515.1/GCA_903899515.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903961455.1/GCA_903961455.1_markerGenes2/COG0541.notab.fna
+ read genome_ID
+ read line
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903843765.1/GCA_903843765.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903913555.1/GCA_903913555.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903842315.1/GCA_903842315.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903835685.1/GCA_903835685.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903841135.1/GCA_903841135.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903884585.1/GCA_903884585.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853455.1/GCA_903853455.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903937725.1/GCA_903937725.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903883715.1/GCA_903883715.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903914905.1/GCA_903914905.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903894295.1/GCA_903894295.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903936775.1/GCA_903936775.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854225.1/GCA_903854225.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903953565.1/GCA_903953565.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915855.1/GCA_903915855.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903857725.1/GCA_903857725.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824635.1/GCA_903824635.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903875175.1/GCA_903875175.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902845.1/GCA_903902845.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853495.1/GCA_903853495.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854255.1/GCA_903854255.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903934285.1/GCA_903934285.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903839445.1/GCA_903839445.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945225.1/GCA_903945225.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903918315.1/GCA_903918315.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903826635.1/GCA_903826635.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915725.1/GCA_903915725.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903924625.1/GCA_903924625.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903879295.1/GCA_903879295.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920995.1/GCA_903920995.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903939965.1/GCA_903939965.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824045.1/GCA_903824045.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903845665.1/GCA_903845665.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907865.1/GCA_903907865.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903900705.1/GCA_903900705.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869725.1/GCA_903869725.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903895255.1/GCA_903895255.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903935095.1/GCA_903935095.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903888055.1/GCA_903888055.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903871585.1/GCA_903871585.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824795.1/GCA_903824795.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907955.1/GCA_903907955.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903959825.1/GCA_903959825.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920605.1/GCA_903920605.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903928895.1/GCA_903928895.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945495.1/GCA_903945495.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869175.1/GCA_903869175.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902895.1/GCA_903902895.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903899515.1/GCA_903899515.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903961455.1/GCA_903961455.1_markerGenes2/COG0552.notab.fna
+ read genome_ID
+ read line
+ mkdir -p test_extension/extension/extension//NEWDB/vsearch/
+ read line
+ python mOTUs-extender/mOTUs-extender//dereplicate_sequences.py test_extension/extension/extension//NEWDB/sequences/COG0012.fna test_extension/extension/extension//NEWDB/sequences/COG0012.derep100.fna test_extension/extension/extension//NEWDB/sequences/COG0012.derep100.clstr
+ vsearch --sortbylength test_extension/extension/extension//NEWDB/sequences/COG0012.derep100.fna --output test_extension/extension/extension//NEWDB/sequences/COG0012.derep100.sorted.fna
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/extension//NEWDB/sequences/COG0012.derep100.fna 100%
32601 nt in 30 seqs, min 792, max 1119, avg 1087
Getting lengths 100%
Sorting 100%
Median length: 1100
Writing output 100%
+ vsearch --threads 32 --allpairs_global test_extension/extension/extension//NEWDB/sequences/COG0012.derep100.sorted.fna --id 0.6 --mincols 20 --blast6out test_extension/extension/extension//NEWDB//vsearch//COG0012.derep100.m8
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/extension//NEWDB/sequences/COG0012.derep100.sorted.fna 100%
32601 nt in 30 seqs, min 792, max 1119, avg 1087
Masking 100%
Aligning 100%
Matching query sequences: 24 of 30 (80.00%)
+ python mOTUs-extender/mOTUs-extender//rereplicate_alignments.py test_extension/extension/extension//NEWDB/sequences/COG0012.derep100.clstr test_extension/extension/extension//NEWDB//vsearch//COG0012.derep100.m8 test_extension/extension/extension//NEWDB//vsearch//COG0012.m8
+ read line
+ python mOTUs-extender/mOTUs-extender//dereplicate_sequences.py test_extension/extension/extension//NEWDB/sequences/COG0016.fna test_extension/extension/extension//NEWDB/sequences/COG0016.derep100.fna test_extension/extension/extension//NEWDB/sequences/COG0016.derep100.clstr
+ vsearch --sortbylength test_extension/extension/extension//NEWDB/sequences/COG0016.derep100.fna --output test_extension/extension/extension//NEWDB/sequences/COG0016.derep100.sorted.fna
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/extension//NEWDB/sequences/COG0016.derep100.fna 100%
33423 nt in 32 seqs, min 900, max 1146, avg 1044
Getting lengths 100%
Sorting 100%
Median length: 1040
Writing output 100%
+ vsearch --threads 32 --allpairs_global test_extension/extension/extension//NEWDB/sequences/COG0016.derep100.sorted.fna --id 0.6 --mincols 20 --blast6out test_extension/extension/extension//NEWDB//vsearch//COG0016.derep100.m8
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/extension//NEWDB/sequences/COG0016.derep100.sorted.fna 100%
33423 nt in 32 seqs, min 900, max 1146, avg 1044
Masking 100%
Aligning 100%
Matching query sequences: 22 of 32 (68.75%)
+ python mOTUs-extender/mOTUs-extender//rereplicate_alignments.py test_extension/extension/extension//NEWDB/sequences/COG0016.derep100.clstr test_extension/extension/extension//NEWDB//vsearch//COG0016.derep100.m8 test_extension/extension/extension//NEWDB//vsearch//COG0016.m8
+ read line
+ python mOTUs-extender/mOTUs-extender//dereplicate_sequences.py test_extension/extension/extension//NEWDB/sequences/COG0018.fna test_extension/extension/extension//NEWDB/sequences/COG0018.derep100.fna test_extension/extension/extension//NEWDB/sequences/COG0018.derep100.clstr
+ vsearch --sortbylength test_extension/extension/extension//NEWDB/sequences/COG0018.derep100.fna --output test_extension/extension/extension//NEWDB/sequences/COG0018.derep100.sorted.fna
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/extension//NEWDB/sequences/COG0018.derep100.fna 100%
53550 nt in 30 seqs, min 1641, max 2013, avg 1785
Getting lengths 100%
Sorting 100%
Median length: 1767
Writing output 100%
+ vsearch --threads 32 --allpairs_global test_extension/extension/extension//NEWDB/sequences/COG0018.derep100.sorted.fna --id 0.6 --mincols 20 --blast6out test_extension/extension/extension//NEWDB//vsearch//COG0018.derep100.m8
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/extension//NEWDB/sequences/COG0018.derep100.sorted.fna 100%
53550 nt in 30 seqs, min 1641, max 2013, avg 1785
Masking 100%
Aligning 100%
Matching query sequences: 16 of 30 (53.33%)
+ python mOTUs-extender/mOTUs-extender//rereplicate_alignments.py test_extension/extension/extension//NEWDB/sequences/COG0018.derep100.clstr test_extension/extension/extension//NEWDB//vsearch//COG0018.derep100.m8 test_extension/extension/extension//NEWDB//vsearch//COG0018.m8
+ read line
+ python mOTUs-extender/mOTUs-extender//dereplicate_sequences.py test_extension/extension/extension//NEWDB/sequences/COG0172.fna test_extension/extension/extension//NEWDB/sequences/COG0172.derep100.fna test_extension/extension/extension//NEWDB/sequences/COG0172.derep100.clstr
+ vsearch --sortbylength test_extension/extension/extension//NEWDB/sequences/COG0172.derep100.fna --output test_extension/extension/extension//NEWDB/sequences/COG0172.derep100.sorted.fna
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/extension//NEWDB/sequences/COG0172.derep100.fna 100%
42252 nt in 33 seqs, min 1266, max 1332, avg 1280
Getting lengths 100%
Sorting 100%
Median length: 1275
Writing output 100%
+ vsearch --threads 32 --allpairs_global test_extension/extension/extension//NEWDB/sequences/COG0172.derep100.sorted.fna --id 0.6 --mincols 20 --blast6out test_extension/extension/extension//NEWDB//vsearch//COG0172.derep100.m8
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/extension//NEWDB/sequences/COG0172.derep100.sorted.fna 100%
42252 nt in 33 seqs, min 1266, max 1332, avg 1280
Masking 100%
Aligning 100%
Matching query sequences: 26 of 33 (78.79%)
+ python mOTUs-extender/mOTUs-extender//rereplicate_alignments.py test_extension/extension/extension//NEWDB/sequences/COG0172.derep100.clstr test_extension/extension/extension//NEWDB//vsearch//COG0172.derep100.m8 test_extension/extension/extension//NEWDB//vsearch//COG0172.m8
+ read line
+ python mOTUs-extender/mOTUs-extender//dereplicate_sequences.py test_extension/extension/extension//NEWDB/sequences/COG0215.fna test_extension/extension/extension//NEWDB/sequences/COG0215.derep100.fna test_extension/extension/extension//NEWDB/sequences/COG0215.derep100.clstr
+ vsearch --sortbylength test_extension/extension/extension//NEWDB/sequences/COG0215.derep100.fna --output test_extension/extension/extension//NEWDB/sequences/COG0215.derep100.sorted.fna
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/extension//NEWDB/sequences/COG0215.derep100.fna 100%
45273 nt in 32 seqs, min 1152, max 1674, avg 1415
Getting lengths 100%
Sorting 100%
Median length: 1404
Writing output 100%
+ vsearch --threads 32 --allpairs_global test_extension/extension/extension//NEWDB/sequences/COG0215.derep100.sorted.fna --id 0.6 --mincols 20 --blast6out test_extension/extension/extension//NEWDB//vsearch//COG0215.derep100.m8
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/extension//NEWDB/sequences/COG0215.derep100.sorted.fna 100%
45273 nt in 32 seqs, min 1152, max 1674, avg 1415
Masking 100%
Aligning 100%
Matching query sequences: 17 of 32 (53.12%)
+ python mOTUs-extender/mOTUs-extender//rereplicate_alignments.py test_extension/extension/extension//NEWDB/sequences/COG0215.derep100.clstr test_extension/extension/extension//NEWDB//vsearch//COG0215.derep100.m8 test_extension/extension/extension//NEWDB//vsearch//COG0215.m8
+ read line
+ python mOTUs-extender/mOTUs-extender//dereplicate_sequences.py test_extension/extension/extension//NEWDB/sequences/COG0495.fna test_extension/extension/extension//NEWDB/sequences/COG0495.derep100.fna test_extension/extension/extension//NEWDB/sequences/COG0495.derep100.clstr
+ vsearch --sortbylength test_extension/extension/extension//NEWDB/sequences/COG0495.derep100.fna --output test_extension/extension/extension//NEWDB/sequences/COG0495.derep100.sorted.fna
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/extension//NEWDB/sequences/COG0495.derep100.fna 100%
85038 nt in 32 seqs, min 2238, max 3405, avg 2657
Getting lengths 100%
Sorting 100%
Median length: 2601
Writing output 100%
+ vsearch --threads 32 --allpairs_global test_extension/extension/extension//NEWDB/sequences/COG0495.derep100.sorted.fna --id 0.6 --mincols 20 --blast6out test_extension/extension/extension//NEWDB//vsearch//COG0495.derep100.m8
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/extension//NEWDB/sequences/COG0495.derep100.sorted.fna 100%
85038 nt in 32 seqs, min 2238, max 3405, avg 2657
Masking 100%
Aligning 100%
Matching query sequences: 21 of 32 (65.62%)
+ python mOTUs-extender/mOTUs-extender//rereplicate_alignments.py test_extension/extension/extension//NEWDB/sequences/COG0495.derep100.clstr test_extension/extension/extension//NEWDB//vsearch//COG0495.derep100.m8 test_extension/extension/extension//NEWDB//vsearch//COG0495.m8
+ read line
+ python mOTUs-extender/mOTUs-extender//dereplicate_sequences.py test_extension/extension/extension//NEWDB/sequences/COG0525.fna test_extension/extension/extension//NEWDB/sequences/COG0525.derep100.fna test_extension/extension/extension//NEWDB/sequences/COG0525.derep100.clstr
+ vsearch --sortbylength test_extension/extension/extension//NEWDB/sequences/COG0525.derep100.fna --output test_extension/extension/extension//NEWDB/sequences/COG0525.derep100.sorted.fna
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/extension//NEWDB/sequences/COG0525.derep100.fna 100%
90765 nt in 34 seqs, min 2517, max 2907, avg 2670
Getting lengths 100%
Sorting 100%
Median length: 2664
Writing output 100%
+ vsearch --threads 32 --allpairs_global test_extension/extension/extension//NEWDB/sequences/COG0525.derep100.sorted.fna --id 0.6 --mincols 20 --blast6out test_extension/extension/extension//NEWDB//vsearch//COG0525.derep100.m8
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/extension//NEWDB/sequences/COG0525.derep100.sorted.fna 100%
90765 nt in 34 seqs, min 2517, max 2907, avg 2670
Masking 100%
Aligning 100%
Matching query sequences: 23 of 34 (67.65%)
+ python mOTUs-extender/mOTUs-extender//rereplicate_alignments.py test_extension/extension/extension//NEWDB/sequences/COG0525.derep100.clstr test_extension/extension/extension//NEWDB//vsearch//COG0525.derep100.m8 test_extension/extension/extension//NEWDB//vsearch//COG0525.m8
+ read line
+ python mOTUs-extender/mOTUs-extender//dereplicate_sequences.py test_extension/extension/extension//NEWDB/sequences/COG0533.fna test_extension/extension/extension//NEWDB/sequences/COG0533.derep100.fna test_extension/extension/extension//NEWDB/sequences/COG0533.derep100.clstr
+ vsearch --sortbylength test_extension/extension/extension//NEWDB/sequences/COG0533.derep100.fna --output test_extension/extension/extension//NEWDB/sequences/COG0533.derep100.sorted.fna
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/extension//NEWDB/sequences/COG0533.derep100.fna 100%
32484 nt in 31 seqs, min 1005, max 1152, avg 1048
Getting lengths 100%
Sorting 100%
Median length: 1044
Writing output 100%
+ vsearch --threads 32 --allpairs_global test_extension/extension/extension//NEWDB/sequences/COG0533.derep100.sorted.fna --id 0.6 --mincols 20 --blast6out test_extension/extension/extension//NEWDB//vsearch//COG0533.derep100.m8
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/extension//NEWDB/sequences/COG0533.derep100.sorted.fna 100%
32484 nt in 31 seqs, min 1005, max 1152, avg 1048
Masking 100%
Aligning 100%
Matching query sequences: 17 of 31 (54.84%)
+ python mOTUs-extender/mOTUs-extender//rereplicate_alignments.py test_extension/extension/extension//NEWDB/sequences/COG0533.derep100.clstr test_extension/extension/extension//NEWDB//vsearch//COG0533.derep100.m8 test_extension/extension/extension//NEWDB//vsearch//COG0533.m8
+ read line
+ python mOTUs-extender/mOTUs-extender//dereplicate_sequences.py test_extension/extension/extension//NEWDB/sequences/COG0541.fna test_extension/extension/extension//NEWDB/sequences/COG0541.derep100.fna test_extension/extension/extension//NEWDB/sequences/COG0541.derep100.clstr
+ vsearch --sortbylength test_extension/extension/extension//NEWDB/sequences/COG0541.derep100.fna --output test_extension/extension/extension//NEWDB/sequences/COG0541.derep100.sorted.fna
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/extension//NEWDB/sequences/COG0541.derep100.fna 100%
42996 nt in 31 seqs, min 1320, max 1515, avg 1387
Getting lengths 100%
Sorting 100%
Median length: 1374
Writing output 100%
+ vsearch --threads 32 --allpairs_global test_extension/extension/extension//NEWDB/sequences/COG0541.derep100.sorted.fna --id 0.6 --mincols 20 --blast6out test_extension/extension/extension//NEWDB//vsearch//COG0541.derep100.m8
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/extension//NEWDB/sequences/COG0541.derep100.sorted.fna 100%
42996 nt in 31 seqs, min 1320, max 1515, avg 1387
Masking 100%
Aligning 100%
Matching query sequences: 25 of 31 (80.65%)
+ python mOTUs-extender/mOTUs-extender//rereplicate_alignments.py test_extension/extension/extension//NEWDB/sequences/COG0541.derep100.clstr test_extension/extension/extension//NEWDB//vsearch//COG0541.derep100.m8 test_extension/extension/extension//NEWDB//vsearch//COG0541.m8
+ read line
+ python mOTUs-extender/mOTUs-extender//dereplicate_sequences.py test_extension/extension/extension//NEWDB/sequences/COG0552.fna test_extension/extension/extension//NEWDB/sequences/COG0552.derep100.fna test_extension/extension/extension//NEWDB/sequences/COG0552.derep100.clstr
+ vsearch --sortbylength test_extension/extension/extension//NEWDB/sequences/COG0552.derep100.fna --output test_extension/extension/extension//NEWDB/sequences/COG0552.derep100.sorted.fna
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/extension//NEWDB/sequences/COG0552.derep100.fna 100%
28803 nt in 30 seqs, min 849, max 1341, avg 960
Getting lengths 100%
Sorting 100%
Median length: 927
Writing output 100%
+ vsearch --threads 32 --allpairs_global test_extension/extension/extension//NEWDB/sequences/COG0552.derep100.sorted.fna --id 0.6 --mincols 20 --blast6out test_extension/extension/extension//NEWDB//vsearch//COG0552.derep100.m8
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/extension//NEWDB/sequences/COG0552.derep100.sorted.fna 100%
28803 nt in 30 seqs, min 849, max 1341, avg 960
Masking 100%
Aligning 100%
Matching query sequences: 21 of 30 (70.00%)
+ python mOTUs-extender/mOTUs-extender//rereplicate_alignments.py test_extension/extension/extension//NEWDB/sequences/COG0552.derep100.clstr test_extension/extension/extension//NEWDB//vsearch//COG0552.derep100.m8 test_extension/extension/extension//NEWDB//vsearch//COG0552.m8
+ read line
+ read line
+ IFS=,
+ set COG0012 94.7
+ python mOTUs-extender/mOTUs-extender//normalizePercentIdentity.py test_extension/extension/extension//NEWDB/vsearch/COG0012.m8 94.7
+ read line
+ IFS=,
+ set COG0016 96.2
+ python mOTUs-extender/mOTUs-extender//normalizePercentIdentity.py test_extension/extension/extension//NEWDB/vsearch/COG0016.m8 96.2
+ read line
+ IFS=,
+ set COG0018 95.5
+ python mOTUs-extender/mOTUs-extender//normalizePercentIdentity.py test_extension/extension/extension//NEWDB/vsearch/COG0018.m8 95.5
+ read line
+ IFS=,
+ set COG0172 95.2
+ python mOTUs-extender/mOTUs-extender//normalizePercentIdentity.py test_extension/extension/extension//NEWDB/vsearch/COG0172.m8 95.2
+ read line
+ IFS=,
+ set COG0215 95.7
+ python mOTUs-extender/mOTUs-extender//normalizePercentIdentity.py test_extension/extension/extension//NEWDB/vsearch/COG0215.m8 95.7
+ read line
+ IFS=,
+ set COG0495 96.0
+ python mOTUs-extender/mOTUs-extender//normalizePercentIdentity.py test_extension/extension/extension//NEWDB/vsearch/COG0495.m8 96.0
+ read line
+ IFS=,
+ set COG0525 95.0
+ python mOTUs-extender/mOTUs-extender//normalizePercentIdentity.py test_extension/extension/extension//NEWDB/vsearch/COG0525.m8 95.0
+ read line
+ IFS=,
+ set COG0533 94.9
+ python mOTUs-extender/mOTUs-extender//normalizePercentIdentity.py test_extension/extension/extension//NEWDB/vsearch/COG0533.m8 94.9
+ read line
+ IFS=,
+ set COG0541 95.9
+ python mOTUs-extender/mOTUs-extender//normalizePercentIdentity.py test_extension/extension/extension//NEWDB/vsearch/COG0541.m8 95.9
+ read line
+ IFS=,
+ set COG0552 95.3
+ python mOTUs-extender/mOTUs-extender//normalizePercentIdentity.py test_extension/extension/extension//NEWDB/vsearch/COG0552.m8 95.3
+ read line
+ read line
+ python mOTUs-extender/mOTUs-extender//filterBlastByAlignLength_Fraction.py test_extension/extension/extension//NEWDB/vsearch//COG0012.normalized.m8 test_extension/extension/extension//NEWDB/NEWDB.all.len -l 20 -f 0.50 -r test_extension/extension/extension//NEWDB/vsearch//COG0012.normalized.excludedPairs
+ read line
+ python mOTUs-extender/mOTUs-extender//filterBlastByAlignLength_Fraction.py test_extension/extension/extension//NEWDB/vsearch//COG0016.normalized.m8 test_extension/extension/extension//NEWDB/NEWDB.all.len -l 20 -f 0.50 -r test_extension/extension/extension//NEWDB/vsearch//COG0016.normalized.excludedPairs
+ read line
+ python mOTUs-extender/mOTUs-extender//filterBlastByAlignLength_Fraction.py test_extension/extension/extension//NEWDB/vsearch//COG0018.normalized.m8 test_extension/extension/extension//NEWDB/NEWDB.all.len -l 20 -f 0.50 -r test_extension/extension/extension//NEWDB/vsearch//COG0018.normalized.excludedPairs
+ read line
+ python mOTUs-extender/mOTUs-extender//filterBlastByAlignLength_Fraction.py test_extension/extension/extension//NEWDB/vsearch//COG0172.normalized.m8 test_extension/extension/extension//NEWDB/NEWDB.all.len -l 20 -f 0.50 -r test_extension/extension/extension//NEWDB/vsearch//COG0172.normalized.excludedPairs
+ read line
+ python mOTUs-extender/mOTUs-extender//filterBlastByAlignLength_Fraction.py test_extension/extension/extension//NEWDB/vsearch//COG0215.normalized.m8 test_extension/extension/extension//NEWDB/NEWDB.all.len -l 20 -f 0.50 -r test_extension/extension/extension//NEWDB/vsearch//COG0215.normalized.excludedPairs
+ read line
+ python mOTUs-extender/mOTUs-extender//filterBlastByAlignLength_Fraction.py test_extension/extension/extension//NEWDB/vsearch//COG0495.normalized.m8 test_extension/extension/extension//NEWDB/NEWDB.all.len -l 20 -f 0.50 -r test_extension/extension/extension//NEWDB/vsearch//COG0495.normalized.excludedPairs
+ read line
+ python mOTUs-extender/mOTUs-extender//filterBlastByAlignLength_Fraction.py test_extension/extension/extension//NEWDB/vsearch//COG0525.normalized.m8 test_extension/extension/extension//NEWDB/NEWDB.all.len -l 20 -f 0.50 -r test_extension/extension/extension//NEWDB/vsearch//COG0525.normalized.excludedPairs
+ read line
+ python mOTUs-extender/mOTUs-extender//filterBlastByAlignLength_Fraction.py test_extension/extension/extension//NEWDB/vsearch//COG0533.normalized.m8 test_extension/extension/extension//NEWDB/NEWDB.all.len -l 20 -f 0.50 -r test_extension/extension/extension//NEWDB/vsearch//COG0533.normalized.excludedPairs
+ read line
+ python mOTUs-extender/mOTUs-extender//filterBlastByAlignLength_Fraction.py test_extension/extension/extension//NEWDB/vsearch//COG0541.normalized.m8 test_extension/extension/extension//NEWDB/NEWDB.all.len -l 20 -f 0.50 -r test_extension/extension/extension//NEWDB/vsearch//COG0541.normalized.excludedPairs
+ read line
+ python mOTUs-extender/mOTUs-extender//filterBlastByAlignLength_Fraction.py test_extension/extension/extension//NEWDB/vsearch//COG0552.normalized.m8 test_extension/extension/extension//NEWDB/NEWDB.all.len -l 20 -f 0.50 -r test_extension/extension/extension//NEWDB/vsearch//COG0552.normalized.excludedPairs
+ read line
+ cat test_extension/extension/extension//NEWDB/vsearch/COG0012.normalized.excludedPairs test_extension/extension/extension//NEWDB/vsearch/COG0016.normalized.excludedPairs test_extension/extension/extension//NEWDB/vsearch/COG0018.normalized.excludedPairs test_extension/extension/extension//NEWDB/vsearch/COG0172.normalized.excludedPairs test_extension/extension/extension//NEWDB/vsearch/COG0215.normalized.excludedPairs test_extension/extension/extension//NEWDB/vsearch/COG0495.normalized.excludedPairs test_extension/extension/extension//NEWDB/vsearch/COG0525.normalized.excludedPairs test_extension/extension/extension//NEWDB/vsearch/COG0533.normalized.excludedPairs test_extension/extension/extension//NEWDB/vsearch/COG0541.normalized.excludedPairs test_extension/extension/extension//NEWDB/vsearch/COG0552.normalized.excludedPairs
+ ls test_extension/extension/extension//NEWDB/vsearch//COG0012.normalized.50p.20n.m8 test_extension/extension/extension//NEWDB/vsearch//COG0016.normalized.50p.20n.m8 test_extension/extension/extension//NEWDB/vsearch//COG0018.normalized.50p.20n.m8 test_extension/extension/extension//NEWDB/vsearch//COG0172.normalized.50p.20n.m8 test_extension/extension/extension//NEWDB/vsearch//COG0215.normalized.50p.20n.m8 test_extension/extension/extension//NEWDB/vsearch//COG0495.normalized.50p.20n.m8 test_extension/extension/extension//NEWDB/vsearch//COG0525.normalized.50p.20n.m8 test_extension/extension/extension//NEWDB/vsearch//COG0533.normalized.50p.20n.m8 test_extension/extension/extension//NEWDB/vsearch//COG0541.normalized.50p.20n.m8 test_extension/extension/extension//NEWDB/vsearch//COG0552.normalized.50p.20n.m8
+ cut -f1 test_extension/extension/extension//NEWDB/NEWDB.new.len
+ sort
+ uniq -d
+ '[' -s test_extension/extension/extension//NEWDB/test_uniq1.txt ']'
+ cat test_extension/extension/extension//NEWDB/NEWDB.new.genomeIDs
+ sort
+ uniq -d
+ '[' -s test_extension/extension/extension//NEWDB/test_uniq2.txt ']'
+ cut -f1 test_extension/extension/extension//NEWDB/NEWDB.new.map2genome
+ sort
+ uniq -d
+ '[' -s test_extension/extension/extension//NEWDB/test_uniq3.txt ']'
+ python mOTUs-extender/mOTUs-extender//combineDistances_5.py -d 55.0 -c 3 test_extension/extension/extension//NEWDB/vsearch/files_normalized.txt test_extension/extension/extension//NEWDB/NEWDB.new.len test_extension/extension/extension//NEWDB/NEWDB.new.genomeIDs test_extension/extension/extension//NEWDB/NEWDB.new.map2genome_v5 test_extension/extension/extension//NEWDB/vsearch/AllCOGs.normalized.excludedPairs test_extension/extension/extension//NEWDB/vsearch/combined.normalized.new.m8
Average OG lengths calculated
COG0525	2664.734693877551
COG0012	1086.84
COG0495	2635.0408163265306
COG0016	1045.0408163265306
COG0533	1056.0
COG0215	1413.3
COG0172	1279.56
COG0541	1386.4285714285713
COG0018	1793.76
COG0552	978.84
reading: test_extension/extension/extension//NEWDB/vsearch//COG0012.normalized.50p.20n.m8
reading: test_extension/extension/extension//NEWDB/vsearch//COG0016.normalized.50p.20n.m8
reading: test_extension/extension/extension//NEWDB/vsearch//COG0018.normalized.50p.20n.m8
reading: test_extension/extension/extension//NEWDB/vsearch//COG0172.normalized.50p.20n.m8
reading: test_extension/extension/extension//NEWDB/vsearch//COG0215.normalized.50p.20n.m8
reading: test_extension/extension/extension//NEWDB/vsearch//COG0495.normalized.50p.20n.m8
reading: test_extension/extension/extension//NEWDB/vsearch//COG0525.normalized.50p.20n.m8
reading: test_extension/extension/extension//NEWDB/vsearch//COG0533.normalized.50p.20n.m8
reading: test_extension/extension/extension//NEWDB/vsearch//COG0541.normalized.50p.20n.m8
reading: test_extension/extension/extension//NEWDB/vsearch//COG0552.normalized.50p.20n.m8
+ python mOTUs-extender/mOTUs-extender//clusterUsingDistsCutoff_4_conComp_3.py test_extension/extension/extension//NEWDB/vsearch/combined.normalized.new.m8 test_extension/extension/extension//NEWDB/NEWDB.new.genomeIDs 99.0 test_extension/extension/extension//NEWDB/vsearch/NEWDB.clustering -d ID -l average
0.01
50
50
Preclustering results: Connected Components: 29	Singletons: 0
Singletons: 16
+ sed -i 's/^/NEWDB_/' test_extension/extension/extension//NEWDB/vsearch/NEWDB.clustering
+ cut -f1 test_extension/extension/extension//NEWDB/vsearch/NEWDB.clustering
+ read line
+ paste test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp1 test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp1
+ sed 's/^/COG0012./'
+ read line
+ paste test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp1 test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp1
+ sed 's/^/COG0016./'
+ read line
+ paste test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp1 test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp1
+ sed 's/^/COG0018./'
+ read line
+ paste test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp1 test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp1
+ sed 's/^/COG0172./'
+ read line
+ paste test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp1 test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp1
+ sed 's/^/COG0215./'
+ read line
+ paste test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp1 test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp1
+ sed 's/^/COG0495./'
+ read line
+ paste test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp1 test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp1
+ sed 's/^/COG0525./'
+ read line
+ paste test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp1 test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp1
+ sed 's/^/COG0533./'
+ read line
+ paste test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp1 test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp1
+ sed 's/^/COG0541./'
+ read line
+ paste test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp1 test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp1
+ sed 's/^/COG0552./'
+ read line
+ awk -F '	' ' { t = $1; $1 = $2; $2 = t; print; } ' 'OFS=	' test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv
+ python mOTUs-extender/mOTUs-extender//reformatMapping2Clustering.py test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.line.tsv
+ rm test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv.temp
+ python mOTUs-extender/mOTUs-extender//reformatClustering2Mapping_3.py test_extension/extension/extension//NEWDB/vsearch/NEWDB.clustering test_extension/extension/extension//NEWDB/vsearch/NEWDB.clustering.map.temp
+ sed 's/\t/,/g' test_extension/extension/extension//NEWDB/vsearch/NEWDB.clustering.map.temp
+ read line
+ IFS=,
+ set NEWDB_1 GCA_903824045.1
+ sed 's/GCA_903824045.1$/NEWDB_1/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_1
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_1 GCA_903888055.1
+ sed 's/GCA_903888055.1$/NEWDB_1/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_1
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_1 GCA_903914905.1
+ sed 's/GCA_903914905.1$/NEWDB_1/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_1
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_2 GCA_903824635.1
+ sed 's/GCA_903824635.1$/NEWDB_2/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_2
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_2 GCA_903857725.1
+ sed 's/GCA_903857725.1$/NEWDB_2/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_2
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_2 GCA_903902845.1
+ sed 's/GCA_903902845.1$/NEWDB_2/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_2
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_2 GCA_903907955.1
+ sed 's/GCA_903907955.1$/NEWDB_2/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_2
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_2 GCA_903937725.1
+ sed 's/GCA_903937725.1$/NEWDB_2/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_2
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_3 GCA_903824795.1
+ sed 's/GCA_903824795.1$/NEWDB_3/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_3
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_3 GCA_903945495.1
+ sed 's/GCA_903945495.1$/NEWDB_3/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_3
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_4 GCA_903826635.1
+ sed 's/GCA_903826635.1$/NEWDB_4/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_4
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_4 GCA_903884585.1
+ sed 's/GCA_903884585.1$/NEWDB_4/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_4
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_5 GCA_903835685.1
+ sed 's/GCA_903835685.1$/NEWDB_5/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_5
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_5 GCA_903879295.1
+ sed 's/GCA_903879295.1$/NEWDB_5/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_5
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_5 GCA_903913555.1
+ sed 's/GCA_903913555.1$/NEWDB_5/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_5
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_6 GCA_903839445.1
+ sed 's/GCA_903839445.1$/NEWDB_6/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_6
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_6 GCA_903907865.1
+ sed 's/GCA_903907865.1$/NEWDB_6/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_6
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_6 GCA_903920995.1
+ sed 's/GCA_903920995.1$/NEWDB_6/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_6
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_7 GCA_903841135.1
+ sed 's/GCA_903841135.1$/NEWDB_7/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_7
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_7 GCA_903915725.1
+ sed 's/GCA_903915725.1$/NEWDB_7/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_7
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_8 GCA_903842315.1
+ sed 's/GCA_903842315.1$/NEWDB_8/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_8
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_8 GCA_903928895.1
+ sed 's/GCA_903928895.1$/NEWDB_8/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_8
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_9 GCA_903845665.1
+ sed 's/GCA_903845665.1$/NEWDB_9/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_9
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_9 GCA_903900705.1
+ sed 's/GCA_903900705.1$/NEWDB_9/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_9
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_10 GCA_903854225.1
+ sed 's/GCA_903854225.1$/NEWDB_10/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_10
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_10 GCA_903869725.1
+ sed 's/GCA_903869725.1$/NEWDB_10/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_10
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_11 GCA_903854255.1
+ sed 's/GCA_903854255.1$/NEWDB_11/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_11
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_11 GCA_903924625.1
+ sed 's/GCA_903924625.1$/NEWDB_11/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_11
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_12 GCA_903875175.1
+ sed 's/GCA_903875175.1$/NEWDB_12/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_12
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_12 GCA_903920605.1
+ sed 's/GCA_903920605.1$/NEWDB_12/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_12
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_13 GCA_903935095.1
+ sed 's/GCA_903935095.1$/NEWDB_13/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_13
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_13 GCA_903939965.1
+ sed 's/GCA_903939965.1$/NEWDB_13/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_13
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_13 GCA_903945225.1
+ sed 's/GCA_903945225.1$/NEWDB_13/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_13
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_13 GCA_903959825.1
+ sed 's/GCA_903959825.1$/NEWDB_13/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_13
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_14 GCA_903843765.1
+ sed 's/GCA_903843765.1$/NEWDB_14/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_14
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_15 GCA_903853455.1
+ sed 's/GCA_903853455.1$/NEWDB_15/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_15
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_16 GCA_903853495.1
+ sed 's/GCA_903853495.1$/NEWDB_16/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_16
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_17 GCA_903869175.1
+ sed 's/GCA_903869175.1$/NEWDB_17/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_17
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_18 GCA_903871585.1
+ sed 's/GCA_903871585.1$/NEWDB_18/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_18
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_19 GCA_903883715.1
+ sed 's/GCA_903883715.1$/NEWDB_19/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_19
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_20 GCA_903894295.1
+ sed 's/GCA_903894295.1$/NEWDB_20/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_20
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_21 GCA_903895255.1
+ sed 's/GCA_903895255.1$/NEWDB_21/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_21
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_22 GCA_903899515.1
+ sed 's/GCA_903899515.1$/NEWDB_22/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_22
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_23 GCA_903902895.1
+ sed 's/GCA_903902895.1$/NEWDB_23/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_23
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_24 GCA_903915855.1
+ sed 's/GCA_903915855.1$/NEWDB_24/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_24
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_25 GCA_903918315.1
+ sed 's/GCA_903918315.1$/NEWDB_25/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_25
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_26 GCA_903934285.1
+ sed 's/GCA_903934285.1$/NEWDB_26/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_26
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_27 GCA_903936775.1
+ sed 's/GCA_903936775.1$/NEWDB_27/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_27
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_28 GCA_903953565.1
+ sed 's/GCA_903953565.1$/NEWDB_28/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_28
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ IFS=,
+ set NEWDB_29 GCA_903961455.1
+ sed 's/GCA_903961455.1$/NEWDB_29/g' test_extension/extension/extension//NEWDB/NEWDB.new.map.temp
+ grep NEWDB_29
+ awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
+ read line
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903843765.1/GCA_903843765.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903913555.1/GCA_903913555.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903842315.1/GCA_903842315.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903835685.1/GCA_903835685.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903841135.1/GCA_903841135.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903884585.1/GCA_903884585.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853455.1/GCA_903853455.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903937725.1/GCA_903937725.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903883715.1/GCA_903883715.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903914905.1/GCA_903914905.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903894295.1/GCA_903894295.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903936775.1/GCA_903936775.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854225.1/GCA_903854225.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903953565.1/GCA_903953565.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915855.1/GCA_903915855.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903857725.1/GCA_903857725.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824635.1/GCA_903824635.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903875175.1/GCA_903875175.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902845.1/GCA_903902845.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853495.1/GCA_903853495.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854255.1/GCA_903854255.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903934285.1/GCA_903934285.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903839445.1/GCA_903839445.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945225.1/GCA_903945225.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903918315.1/GCA_903918315.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903826635.1/GCA_903826635.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915725.1/GCA_903915725.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903924625.1/GCA_903924625.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903879295.1/GCA_903879295.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920995.1/GCA_903920995.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903939965.1/GCA_903939965.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824045.1/GCA_903824045.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903845665.1/GCA_903845665.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907865.1/GCA_903907865.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903900705.1/GCA_903900705.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869725.1/GCA_903869725.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903895255.1/GCA_903895255.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903935095.1/GCA_903935095.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903888055.1/GCA_903888055.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903871585.1/GCA_903871585.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824795.1/GCA_903824795.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907955.1/GCA_903907955.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903959825.1/GCA_903959825.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920605.1/GCA_903920605.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903928895.1/GCA_903928895.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945495.1/GCA_903945495.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869175.1/GCA_903869175.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902895.1/GCA_903902895.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903899515.1/GCA_903899515.1.padded.fasta
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903961455.1/GCA_903961455.1.padded.fasta
+ read genome_ID
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903843765.1/GCA_903843765.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903913555.1/GCA_903913555.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903842315.1/GCA_903842315.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903835685.1/GCA_903835685.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903841135.1/GCA_903841135.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903884585.1/GCA_903884585.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853455.1/GCA_903853455.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903937725.1/GCA_903937725.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903883715.1/GCA_903883715.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903914905.1/GCA_903914905.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903894295.1/GCA_903894295.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903936775.1/GCA_903936775.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854225.1/GCA_903854225.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903953565.1/GCA_903953565.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915855.1/GCA_903915855.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903857725.1/GCA_903857725.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824635.1/GCA_903824635.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903875175.1/GCA_903875175.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902845.1/GCA_903902845.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903853495.1/GCA_903853495.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903854255.1/GCA_903854255.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903934285.1/GCA_903934285.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903839445.1/GCA_903839445.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945225.1/GCA_903945225.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903918315.1/GCA_903918315.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903826635.1/GCA_903826635.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903915725.1/GCA_903915725.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903924625.1/GCA_903924625.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903879295.1/GCA_903879295.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920995.1/GCA_903920995.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903939965.1/GCA_903939965.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824045.1/GCA_903824045.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903845665.1/GCA_903845665.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907865.1/GCA_903907865.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903900705.1/GCA_903900705.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869725.1/GCA_903869725.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903895255.1/GCA_903895255.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903935095.1/GCA_903935095.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903888055.1/GCA_903888055.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903871585.1/GCA_903871585.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903824795.1/GCA_903824795.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903907955.1/GCA_903907955.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903959825.1/GCA_903959825.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903920605.1/GCA_903920605.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903928895.1/GCA_903928895.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903945495.1/GCA_903945495.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903869175.1/GCA_903869175.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903902895.1/GCA_903902895.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903899515.1/GCA_903899515.1.padded.coords
+ read genome_ID
+ cat test_extension/extension/extension//dbs/GCA_903961455.1/GCA_903961455.1.padded.coords
+ read genome_ID
+ python mOTUs-extender/mOTUs-extender//postprocess_min1.py test_extension/extension/temp_db_folder// test_extension/extension/extension//NEWDB/ mOTUs-extender/mOTUs-extender//cutoffs_fscore_specIAsRef.csv 32
vsearch --threads 32 --usearch_global test_extension/extension/temp_db_folder///min1_db_mOTU/db_mOTU_DB_NR_unpadded.fasta --db test_extension/extension/extension//NEWDB/NEWDB.padded --id 0.8  --maxaccepts 100 --maxrejects 100 --mincols 20 --alnout test_extension/extension/extension//NEWDB/NEWDB.padded_min1.m8.aln --blast6out test_extension/extension/extension//NEWDB/NEWDB.padded_min1.m8
vsearch v2.14.1_linux_x86_64, 2011.3GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file test_extension/extension/extension//NEWDB/NEWDB.padded 100%
857877 nt in 496 seqs, min 909, max 3605, avg 1730
Masking 100%
Counting k-mers 100%
Creating k-mer index 100%
Searching 4%
Searching 14%
Searching 100%
Matching unique query sequences: 999 of 222782 (0.45%)
17 MGCs will be removed from the -1 mOTU
17 MGs will be removed from the -1 mOTU
Writing test_extension/extension/extension//NEWDB//updated_min1_db_mOTU/db_mOTU_MAP_genes_to_MGCs.tsv
Wrote 222765 MGs. Removed 17 MGs
################################################################
Writing test_extension/extension/extension//NEWDB//updated_min1_db_mOTU/db_mOTU_MAP_MGCs_to_mOTUs.tsv
Wrote 213635 MGCs
################################################################
Writing test_extension/extension/extension//NEWDB//updated_min1_db_mOTU/db_mOTU_MAP_MGCs_to_mOTUs_in-line.tsv
################################################################
Writing test_extension/extension/extension//NEWDB//updated_min1_db_mOTU/db_mOTU_padding_coordinates_NR.tsv
Written MGs 222765
Removed MGs 17
################################################################
Writing test_extension/extension/extension//NEWDB//updated_min1_db_mOTU/db_mOTU_genes_length_NR
Written MGs 222765
Removed MGs 17
################################################################
Writing test_extension/extension/extension//NEWDB//updated_min1_db_mOTU/db_mOTU_DB_NR.fasta
Written MGs 222765
Removed MGs 17
################################################################
Writing test_extension/extension/extension//NEWDB//updated_min1_db_mOTU//db_mOTU_padding_coordinates_CEN.tsv
Written MGs 213635
Removed MGs 17 --> Numbers are different to the previous ones as these are centroid sequences
################################################################
Writing test_extension/extension/extension//NEWDB//updated_min1_db_mOTU//db_mOTU_DB_CEN.fasta.annotations
Written MGs 213635
Removed MGs 17 --> Numbers are different to the previous ones as these are centroid sequences
################################################################
Writing test_extension/extension/temp_db_folder///min1_db_mOTU//db_mOTU_DB_CEN.fasta
Writing test_extension/extension/extension//NEWDB//updated_min1_db_mOTU//db_mOTU_DB_CEN.fasta
Written MGs 213635
Removed MGs 17
################################################################
Succesfully finished
+ python mOTUs-extender/mOTUs-extender//extend_mOTUs_getClusterTaxonomy.py test_extension/extension/extension//NEWDB/vsearch/NEWDB.clustering mOTUs-extender/mOTUs-extender/test/genomes.tax
+ cp test_extension/extension/extension//NEWDB/vsearch/NEWDB.clustering test_extension/extension/extension//NEWDB/
+ cp test_extension/extension/extension//NEWDB/vsearch/NEWDB.taxonomy test_extension/extension/extension//NEWDB/
+ cp test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.tsv test_extension/extension/extension//NEWDB/
+ mkdir -p test_extension/extension/extension//NEWDB/db_mOTU/
+ cp test_extension/extension/temp_db_folder//db_mOTU/db_mOTU_taxonomy_ref-mOTUs_short_names.tsv test_extension/extension/extension//NEWDB/db_mOTU/
+ cp test_extension/extension/temp_db_folder//db_mOTU/db_mOTU_taxonomy_ref-mOTUs.tsv test_extension/extension/extension//NEWDB/db_mOTU/
+ cp mOTUs-extender/mOTUs-extender//versions test_extension/extension/extension//NEWDB/db_mOTU/db_mOTU_versions
+ cat test_extension/extension/temp_db_folder//db_mOTU/db_mOTU_DB_CEN.fasta test_extension/extension/extension//NEWDB/updated_min1_db_mOTU/db_mOTU_DB_CEN.fasta
+ cat test_extension/extension/temp_db_folder//db_mOTU/db_mOTU_DB_CEN.fasta.annotations test_extension/extension/extension//NEWDB/updated_min1_db_mOTU/db_mOTU_DB_CEN.fasta.annotations
+ cat test_extension/extension/temp_db_folder//db_mOTU/db_mOTU_DB_NR.fasta test_extension/extension/extension//NEWDB/updated_min1_db_mOTU/db_mOTU_DB_NR.fasta test_extension/extension/extension//NEWDB/NEWDB.padded
+ cat test_extension/extension/temp_db_folder//db_mOTU/db_mOTU_genes_length_NR test_extension/extension/extension//NEWDB/updated_min1_db_mOTU/db_mOTU_genes_length_NR test_extension/extension/extension//NEWDB/NEWDB.new.len
+ cat test_extension/extension/temp_db_folder//db_mOTU/db_mOTU_MAP_genes_to_MGCs.tsv test_extension/extension/extension//NEWDB/updated_min1_db_mOTU/db_mOTU_MAP_genes_to_MGCs.tsv test_extension/extension/extension//NEWDB/NEWDB.map
+ cat test_extension/extension/temp_db_folder//db_mOTU/db_mOTU_MAP_MGCs_to_mOTUs_in-line.tsv test_extension/extension/extension//NEWDB/updated_min1_db_mOTU/db_mOTU_MAP_MGCs_to_mOTUs_in-line.tsv test_extension/extension/extension//NEWDB/vsearch/NEWDB.mOTU-LG.map.line.tsv
+ cat test_extension/extension/temp_db_folder//db_mOTU/db_mOTU_MAP_MGCs_to_mOTUs.tsv test_extension/extension/extension//NEWDB/updated_min1_db_mOTU/db_mOTU_MAP_MGCs_to_mOTUs.tsv test_extension/extension/extension//NEWDB//NEWDB.mOTU-LG.map.tsv
+ cat test_extension/extension/temp_db_folder//db_mOTU/db_mOTU_padding_coordinates_CEN.tsv test_extension/extension/extension//NEWDB/updated_min1_db_mOTU/db_mOTU_padding_coordinates_CEN.tsv
+ cat test_extension/extension/temp_db_folder//db_mOTU/db_mOTU_padding_coordinates_NR.tsv test_extension/extension/extension//NEWDB/updated_min1_db_mOTU/db_mOTU_padding_coordinates_NR.tsv test_extension/extension/extension//NEWDB/NEWDB.padded.coords
+ cat test_extension/extension/temp_db_folder//db_mOTU/db_mOTU_taxonomy_meta-mOTUs.tsv test_extension/extension/extension//NEWDB//NEWDB.taxonomy
+ bwa index -b 1000000000 test_extension/extension/extension//NEWDB/db_mOTU/db_mOTU_DB_CEN.fasta
[bwa_index] Pack FASTA... 6.64 sec
[bwa_index] Construct BWT for the packed sequence...
[BWTIncCreate] textLength=1610445818, availableWord=1313234988
[bwt_gen] Finished constructing BWT in 8 iterations.
[bwa_index] 743.16 seconds elapse.
[bwa_index] Update BWT... 4.73 sec
[bwa_index] Pack forward-only FASTA... 3.66 sec
[bwa_index] Construct SA from BWT and Occ... 222.22 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index -b 1000000000 test_extension/extension/extension//NEWDB/db_mOTU/db_mOTU_DB_CEN.fasta
[main] Real time: 1075.250 sec; CPU: 980.404 sec
+ bwa index -b 1000000000 test_extension/extension/extension//NEWDB/db_mOTU/db_mOTU_DB_NR.fasta
[bwa_index] Pack FASTA... 15.90 sec
[bwa_index] Construct BWT for the packed sequence...
[BWTIncCreate] textLength=3880994364, availableWord=1472883968
[BWTIncConstructFromPacked] 10 iterations done. 2429666412 characters processed.
[bwt_gen] Finished constructing BWT in 17 iterations.
[bwa_index] 2184.21 seconds elapse.
[bwa_index] Update BWT... 13.07 sec
[bwa_index] Pack forward-only FASTA... 8.41 sec
[bwa_index] Construct SA from BWT and Occ... 494.89 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index -b 1000000000 test_extension/extension/extension//NEWDB/db_mOTU/db_mOTU_DB_NR.fasta
[main] Real time: 2734.442 sec; CPU: 2716.497 sec
+ bwa mem test_extension/extension/extension//NEWDB/db_mOTU/db_mOTU_DB_NR.fasta test_extension/extension/temp_db_folder//db_mOTU/db_mOTU_test/test1_single.fastq
+ grep '^@SQ'
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 2495 sequences (411507 bp)...
[M::mem_process_seqs] Processed 2495 reads in 1.291 CPU sec, 1.322 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem test_extension/extension/extension//NEWDB/db_mOTU/db_mOTU_DB_NR.fasta test_extension/extension/temp_db_folder//db_mOTU/db_mOTU_test/test1_single.fastq
[main] Real time: 4.820 sec; CPU: 4.341 sec
+ bwa mem test_extension/extension/extension//NEWDB/db_mOTU/db_mOTU_DB_CEN.fasta test_extension/extension/temp_db_folder//db_mOTU/db_mOTU_test/test1_single.fastq
+ grep '^@SQ'
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 2495 sequences (411507 bp)...
[M::mem_process_seqs] Processed 2495 reads in 0.867 CPU sec, 0.870 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem test_extension/extension/extension//NEWDB/db_mOTU/db_mOTU_DB_CEN.fasta test_extension/extension/temp_db_folder//db_mOTU/db_mOTU_test/test1_single.fastq
[main] Real time: 2.400 sec; CPU: 2.151 sec
+ python mOTUs-extender/mOTUs-extender//create_cami_tax.py test_extension/extension/temp_db_folder//db_mOTU/db_mOTU_taxonomy_CAMI.tsv test_extension/extension/extension//NEWDB/db_mOTU/db_mOTU_taxonomy_meta-mOTUs.tsv test_extension/extension/extension//NEWDB/db_mOTU/db_mOTU_taxonomy_CAMI.tsv NEWDB
+ echo 'New database deposited in: test_extension/extension/extension//NEWDB/db_mOTU/'
New database deposited in: test_extension/extension/extension//NEWDB/db_mOTU/
2022-08-12 16:06:53,695 INFO: Executing: sed -i 's/manual_update/3.0.1_NEWDB/g' test_extension/extension/extension/NEWDB/db_mOTU/db_mOTU_versions
2022-08-12 16:06:53,736 INFO: mOTUs-extender finished successfully

```
</details>

The extended database (`test_extension/extension/extension/NEWDB/db_mOTU/`) is now ready to profile metagenomic sequencing data using the `motus profile` routine (see [here](https://github.com/motu-tool/mOTUs/wiki))

Genome to mOTU associations can be found here:

```
cat test_extension/extension/extension/NEWDB/NEWDB.clustering

NEWDB_1		GCA_903824045.1;GCA_903888055.1;GCA_903914905.1
NEWDB_2		GCA_903824635.1;GCA_903857725.1;GCA_903902845.1;GCA_903907955.1;GCA_903937725.1
NEWDB_3		GCA_903824795.1;GCA_903945495.1
NEWDB_4		GCA_903826635.1;GCA_903884585.1
NEWDB_5		GCA_903835685.1;GCA_903879295.1;GCA_903913555.1
NEWDB_6		GCA_903839445.1;GCA_903907865.1;GCA_903920995.1
NEWDB_7		GCA_903841135.1;GCA_903915725.1
NEWDB_8		GCA_903842315.1;GCA_903928895.1
NEWDB_9		GCA_903845665.1;GCA_903900705.1
NEWDB_10	GCA_903854225.1;GCA_903869725.1
NEWDB_11	GCA_903854255.1;GCA_903924625.1
NEWDB_12	GCA_903875175.1;GCA_903920605.1
NEWDB_13	GCA_903935095.1;GCA_903939965.1;GCA_903945225.1;GCA_903959825.1
NEWDB_14	GCA_903843765.1
NEWDB_15	GCA_903853455.1
NEWDB_16	GCA_903853495.1
NEWDB_17	GCA_903869175.1
NEWDB_18	GCA_903871585.1
NEWDB_19	GCA_903883715.1
NEWDB_20	GCA_903894295.1
NEWDB_21	GCA_903895255.1
NEWDB_22	GCA_903899515.1
NEWDB_23	GCA_903902895.1
NEWDB_24	GCA_903915855.1
NEWDB_25	GCA_903918315.1
NEWDB_26	GCA_903934285.1
NEWDB_27	GCA_903936775.1
NEWDB_28	GCA_903953565.1
NEWDB_29	GCA_903961455.1
```