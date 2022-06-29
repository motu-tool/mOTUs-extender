import glob
import pprint
import os

script_folder = os.path.normpath(config['scriptfolder']) + '/'
input_path = os.path.normpath(config['infolder']) + '/'
out_path = os.path.normpath(config['outfolder']) + '/'
temp_db_folder_path = os.path.normpath(config['temp_db_folder']) + '/'
dest_path = out_path + 'genes/'
dest_path_extension = out_path + 'extension/'



maxcount = 100000000000000



donefiles = set(glob.glob(dest_path + '*prodigal_fetchMGs.done'))
genomes = sorted(glob.glob(input_path + '*fa'))
genomes = genomes[:maxcount]

prodigal_fetchMGs_marker_files = []
for cnt, genome in enumerate(genomes):
    destpath = genome.replace(input_path, dest_path).rsplit('/', 1)[0] + '/' + genome.split('/')[-1][:-3]
    samplename = genome.split('/')[-1][:-3]
    marker = destpath + '/' + samplename + '.prodigal_fetchMGs.done'
    if marker not in donefiles:
        prodigal_fetchMGs_marker_files.append(marker)

rule all:
    input:
        prodigal_fetchMGs_marker_files

rule prodigalfetchMgs:
    input:
        fa = input_path + '{sample}.fa'
    output:
        marker = touch(dest_path + '{sample}/{sample}.prodigal_fetchMGs.done')
    params:
        outdir = dest_path + '{sample}/',
        faa = dest_path + '{sample}/prodigal/{sample}.faa',
        fna = dest_path + '{sample}/prodigal/{sample}.fna',
        gff = dest_path + '{sample}/prodigal/{sample}.gff',
        
    log:
        log = dest_path + '{sample}/{sample}.prodigal_fetchMGs.log',
        command = dest_path + '{sample}/{sample}.prodigal_fetchMGs.command'
    benchmark:
        dest_path + '{sample}/{sample}.prodigal_fetchMGs.benchmark'
    threads:
        1
    shell:
        '''
        command="
        mkdir -p {params.outdir}genome/
        mkdir -p {params.outdir}prodigal/
        mkdir -p {params.outdir}fetchMGs/
        rsync {input.fa} {params.outdir}genome/
       	prodigal -i {input.fa} -a {params.faa} -d {params.fna} -f gff -o {params.gff} -c -m -g 11 -p single
        rm -rf {params.outdir}fetchMGs/{wildcards.sample}-bestMGs/
        fetchMGs.pl -m extraction -v -i -d {params.fna} -o {params.outdir}fetchMGs/{wildcards.sample}-bestMGs/ {params.faa} -t {threads}
        mkdir -p {params.outdir}fetchMGs/{wildcards.sample}-bestMGs-renamed/
        python {script_folder}/parse_fetchMGs.py {params.outdir}fetchMGs/{wildcards.sample}-bestMGs/ {params.outdir}fetchMGs/{wildcards.sample}-bestMGs-renamed/ {wildcards.sample}
        rm -r {params.outdir}fetchMGs/{wildcards.sample}-bestMGs/hmmResults/
        rm -r {params.outdir}fetchMGs/{wildcards.sample}-bestMGs/temp/
        if [ -f {params.outdir}fetchMGs/{wildcards.sample}-bestMGs-renamed/motus.mgs.count.ok ]; then
            sh {script_folder}/extend_mOTUs_addMarkerGenes.sh {input.fa} {wildcards.sample} {dest_path_extension} {script_folder} {temp_db_folder_path} {params.outdir}fetchMGs/{wildcards.sample}-bestMGs-renamed/
        fi
        ";
        echo "$command" > {log.command};
        eval "$command" &> {log.log}
        '''
