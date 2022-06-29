set -euxo pipefail
#filter genomes from mOTU DB


fileListGenomes=$1
newDBName=$2
taxonomyFile=$3
new_database_folder=$4
scriptDir=$5
mOTU_folder=$6
threads=$7
cutoffsFile="$scriptDir/cutoffs_fscore_specIAsRef.csv"
mOTU_MG_file="$scriptDir/10RefMGs.IDs"
cutoff=99.0


fileListGenomes2="$new_database_folder/$newDBName.newGenomes.filteredMissingMGs.txt"
touch $fileListGenomes2 && rm $fileListGenomes2
touch $fileListGenomes2

while read genome_ID;
do
if [ -s "$new_database_folder/dbs/${genome_ID}/${genome_ID}.map" ]
then
echo $genome_ID >> $fileListGenomes2
else
echo "Excluding ${genome_ID} because it's lacking all 10 MGs used for mOTUs"
fi
done < $fileListGenomes

if [ -s "$fileListGenomes2" ]
then
echo "Adding genomes to database."
else
echo "No genomes to be added. Exiting."
exit 0
fi

mkdir -p $new_database_folder/$newDBName/vsearch

while read genome_ID;
do
cat $new_database_folder/dbs/${genome_ID}/vsearch/combined.normalized.distances_vs_db.m8
done < $fileListGenomes2 >> $new_database_folder/$newDBName/vsearch/combined.normalized.distances_vs_db.m8


#script to filter out genomes that are too close to existing mOTUs (just average best hit distance to mOTU)
python $scriptDir/extend_mOTUs_filterGenomesByDistance.py --cutoff $cutoff $new_database_folder/$newDBName/vsearch/combined.normalized.distances_vs_db.m8 $fileListGenomes2 > $new_database_folder/$newDBName/genomes.filtered.list

if [ -s $new_database_folder/$newDBName/genomes.filtered.list ]
then
echo "Adding genomes listed in $new_database_folder/$newDBName/genomes.filtered.list."
else
echo "All genomes are already represented in the mOTU database. Exiting."
exit 1
fi


while read genome_ID;
do
cat $new_database_folder/dbs/${genome_ID}/vsearch/combined.normalized.distances_vs_db.m8
done < $new_database_folder/$newDBName/genomes.filtered.list > $new_database_folder/$newDBName/vsearch/combined.normalized.distances_vs_db.m8

while read genome_ID;
do
cat $new_database_folder/dbs/${genome_ID}/${genome_ID}.map
done < $new_database_folder/$newDBName/genomes.filtered.list > $new_database_folder/$newDBName/${newDBName}.new.map.temp


while read genome_ID;
do
cat $new_database_folder/dbs/${genome_ID}/${genome_ID}.map2genome
done < $new_database_folder/$newDBName/genomes.filtered.list > $new_database_folder/$newDBName/${newDBName}.new.map2genome
cut -f2 $new_database_folder/$newDBName/${newDBName}.new.map2genome | sort -u > $new_database_folder/$newDBName/${newDBName}.new.genomeIDs

while read genome_ID;
do
cat $new_database_folder/dbs/${genome_ID}/${genome_ID}.genes.len
done < $new_database_folder/$newDBName/genomes.filtered.list > $new_database_folder/$newDBName/${newDBName}.new.len

cat $new_database_folder/$newDBName/${newDBName}.new.len $mOTU_folder/db_mOTU/db_mOTU_genes_length_NR > $new_database_folder/$newDBName/${newDBName}.all.len

mkdir -p $new_database_folder/$newDBName/sequences/
while read line;
do
while read genome_ID;
do
cat $new_database_folder/dbs/${genome_ID}/${genome_ID}_markerGenes2/$line.notab.fna
done < $new_database_folder/$newDBName/genomes.filtered.list > $new_database_folder/$newDBName/sequences/$line.fna
done < $mOTU_MG_file



mkdir -p $new_database_folder/$newDBName/vsearch/
while read line;
do
python $scriptDir/dereplicate_sequences.py $new_database_folder/$newDBName/sequences/$line.fna $new_database_folder/$newDBName/sequences/$line.derep100.fna $new_database_folder/$newDBName/sequences/$line.derep100.clstr
vsearch --sortbylength $new_database_folder/$newDBName/sequences/$line.derep100.fna --output $new_database_folder/$newDBName/sequences/$line.derep100.sorted.fna
vsearch --threads ${threads} --allpairs_global $new_database_folder/$newDBName/sequences/$line.derep100.sorted.fna --id 0.0 --mincols 20 --blast6out $new_database_folder/$newDBName//vsearch//$line.derep100.m8
python $scriptDir/rereplicate_alignments.py $new_database_folder/$newDBName/sequences/$line.derep100.clstr  $new_database_folder/$newDBName//vsearch//$line.derep100.m8 $new_database_folder/$newDBName//vsearch//$line.m8
done < $mOTU_MG_file


while read line;
do
IFS=","; set $line
python $scriptDir/normalizePercentIdentity.py $new_database_folder/$newDBName/vsearch/$1.m8 $2 > $new_database_folder/$newDBName/vsearch//$1.normalized.m8
done < $cutoffsFile

while read line;
do
#echo $line
python $scriptDir/filterBlastByAlignLength_Fraction.py $new_database_folder/$newDBName/vsearch//$line.normalized.m8 $new_database_folder/$newDBName/${newDBName}.all.len -l 20 -f 0.50 -r $new_database_folder/$newDBName/vsearch//$line.normalized.excludedPairs > $new_database_folder/$newDBName/vsearch//$line.normalized.50p.20n.m8
done < $mOTU_MG_file

cat $new_database_folder/$newDBName/vsearch/COG*.normalized.excludedPairs  > $new_database_folder/$newDBName/vsearch/AllCOGs.normalized.excludedPairs
ls $new_database_folder/$newDBName/vsearch//COG*.normalized.50p.20n.m8 > $new_database_folder/$newDBName/vsearch/files_normalized.txt

cut -f1 $new_database_folder/$newDBName/${newDBName}.new.len | sort | uniq -d > $new_database_folder/$newDBName/test_uniq1.txt
if [ -s $new_database_folder/$newDBName/test_uniq1.txt ]
then
echo "Duplicated entry in $new_database_folder/$newDBName/${newDBName}.new.len! Please check your input files (also see $new_database_folder/$newDBName/test_uniq1.txt)."
exit 1
fi

cat $new_database_folder/$newDBName/${newDBName}.new.genomeIDs | sort | uniq -d > $new_database_folder/$newDBName/test_uniq2.txt
if [ -s $new_database_folder/$newDBName/test_uniq2.txt ]
then
echo "Duplicated entry in $new_database_folder/$newDBName/${newDBName}.new.genomeIDs! Please check your input files (also see $new_database_folder/$newDBName/test_uniq2.txt)."
exit 1
fi

cut -f1 $new_database_folder/$newDBName/${newDBName}.new.map2genome | sort | uniq -d > $new_database_folder/$newDBName/test_uniq3.txt
if [ -s $new_database_folder/$newDBName/test_uniq3.txt ]
then
echo "Duplicated entry in new_database_folder/$newDBName/${newDBName}.new.map2genome! Please check your input files (also see $new_database_folder/$newDBName/test_uniq3.txt)."
exit 1
fi

python $scriptDir/combineDistances_4_useMap.py $new_database_folder/$newDBName/vsearch/files_normalized.txt $new_database_folder/$newDBName/${newDBName}.new.len $new_database_folder/$newDBName/${newDBName}.new.genomeIDs $new_database_folder/$newDBName/${newDBName}.new.map2genome $new_database_folder/$newDBName/vsearch/AllCOGs.normalized.excludedPairs $new_database_folder/$newDBName/vsearch/combined.normalized.new.m8


python $scriptDir/clusterUsingDistsCutoff_4_conComp_3.py $new_database_folder/$newDBName/vsearch/combined.normalized.new.m8 $new_database_folder/$newDBName/${newDBName}.new.genomeIDs $cutoff $new_database_folder/$newDBName/vsearch/$newDBName.clustering -d ID -l average

sed -i "s/^/${newDBName}_/" $new_database_folder/$newDBName/vsearch/$newDBName.clustering

cut -f1 $new_database_folder/$newDBName/vsearch/$newDBName.clustering > $new_database_folder/$newDBName/vsearch/$newDBName.mOTU-LG.map.tsv.temp1

while read line;
do
paste $new_database_folder/$newDBName/vsearch/$newDBName.mOTU-LG.map.tsv.temp1 $new_database_folder/$newDBName/vsearch/$newDBName.mOTU-LG.map.tsv.temp1 | sed "s/^/${line}./"
done < $mOTU_MG_file > $new_database_folder/$newDBName/vsearch/$newDBName.mOTU-LG.map.tsv

awk -F $'\t' ' { t = $1; $1 = $2; $2 = t; print; } ' OFS=$'\t' $new_database_folder/$newDBName/vsearch/$newDBName.mOTU-LG.map.tsv > $new_database_folder/$newDBName/vsearch/$newDBName.mOTU-LG.map.tsv.temp

python $scriptDir/reformatMapping2Clustering.py $new_database_folder/$newDBName/vsearch/$newDBName.mOTU-LG.map.tsv.temp $new_database_folder/$newDBName/vsearch/$newDBName.mOTU-LG.map.line.tsv

rm $new_database_folder/$newDBName/vsearch/$newDBName.mOTU-LG.map.tsv.temp

python $scriptDir/reformatClustering2Mapping_3.py $new_database_folder/$newDBName/vsearch/$newDBName.clustering $new_database_folder/$newDBName/vsearch/$newDBName.clustering.map.temp
sed "s/\t/,/g" $new_database_folder/$newDBName/vsearch/$newDBName.clustering.map.temp > $new_database_folder/$newDBName/vsearch/$newDBName.clustering.map.temp2


#### make taxonomy file (need script, choices: consensus, majority, all)

#this should be recoded
while read line;
do
IFS=","; set $line
sed "s/$2$/$1/g" $new_database_folder/$newDBName/${newDBName}.new.map.temp | grep $1 | awk '{print $1"\t"$2"\t"$3"\t"$3"."$4}'
done < $new_database_folder/$newDBName/vsearch/$newDBName.clustering.map.temp2 > $new_database_folder/$newDBName/${newDBName}.map


####copy all parts together into a new DB folder

while read genome_ID;
do
cat $new_database_folder/dbs/$genome_ID/$genome_ID.padded.fasta
done < $new_database_folder/$newDBName/genomes.filtered.list > $new_database_folder/$newDBName/$newDBName.padded


while read genome_ID;
do
cat $new_database_folder/dbs/$genome_ID/$genome_ID.padded.coords
done < $new_database_folder/$newDBName/genomes.filtered.list  > $new_database_folder/$newDBName/$newDBName.padded.coords



python $scriptDir/postprocess_min1.py $mOTU_folder/ $new_database_folder/$newDBName/ $scriptDir/cutoffs_fscore_specIAsRef.csv $threads

python $scriptDir/extend_mOTUs_getClusterTaxonomy.py $new_database_folder/$newDBName/vsearch/$newDBName.clustering $taxonomyFile > $new_database_folder/$newDBName/vsearch/$newDBName.taxonomy

cp $new_database_folder/$newDBName/vsearch/$newDBName.clustering $new_database_folder/$newDBName/
cp $new_database_folder/$newDBName/vsearch/$newDBName.taxonomy $new_database_folder/$newDBName/
cp $new_database_folder/$newDBName/vsearch/$newDBName.mOTU-LG.map.tsv $new_database_folder/$newDBName/

mkdir -p $new_database_folder/$newDBName/db_mOTU/


cp $mOTU_folder/db_mOTU/db_mOTU_taxonomy_ref-mOTUs_short_names.tsv $new_database_folder/$newDBName/db_mOTU/
cp $mOTU_folder/db_mOTU/db_mOTU_taxonomy_ref-mOTUs.tsv $new_database_folder/$newDBName/db_mOTU/
cp $scriptDir/versions $new_database_folder/$newDBName/db_mOTU/db_mOTU_versions








#TODO: centroids of extended db

cat $mOTU_folder/db_mOTU/db_mOTU_DB_CEN.fasta $new_database_folder/$newDBName/updated_min1_db_mOTU/db_mOTU_DB_CEN.fasta > $new_database_folder/$newDBName/db_mOTU/db_mOTU_DB_CEN.fasta
cat $mOTU_folder/db_mOTU/db_mOTU_DB_CEN.fasta.annotations $new_database_folder/$newDBName/updated_min1_db_mOTU/db_mOTU_DB_CEN.fasta.annotations > $new_database_folder/$newDBName/db_mOTU/db_mOTU_DB_CEN.fasta.annotations
cat $mOTU_folder/db_mOTU/db_mOTU_DB_NR.fasta $new_database_folder/$newDBName/updated_min1_db_mOTU/db_mOTU_DB_NR.fasta $new_database_folder/$newDBName/$newDBName.padded > $new_database_folder/$newDBName/db_mOTU/db_mOTU_DB_NR.fasta
cat $mOTU_folder/db_mOTU/db_mOTU_genes_length_NR $new_database_folder/$newDBName/updated_min1_db_mOTU/db_mOTU_genes_length_NR $new_database_folder/$newDBName/${newDBName}.new.len > $new_database_folder/$newDBName/db_mOTU/db_mOTU_genes_length_NR
cat $mOTU_folder/db_mOTU/db_mOTU_MAP_genes_to_MGCs.tsv $new_database_folder/$newDBName/updated_min1_db_mOTU/db_mOTU_MAP_genes_to_MGCs.tsv $new_database_folder/$newDBName/${newDBName}.map > $new_database_folder/$newDBName/db_mOTU/db_mOTU_MAP_genes_to_MGCs.tsv
cat $mOTU_folder/db_mOTU/db_mOTU_MAP_MGCs_to_mOTUs_in-line.tsv $new_database_folder/$newDBName/updated_min1_db_mOTU/db_mOTU_MAP_MGCs_to_mOTUs_in-line.tsv $new_database_folder/$newDBName/vsearch/$newDBName.mOTU-LG.map.line.tsv > $new_database_folder/$newDBName/db_mOTU/db_mOTU_MAP_MGCs_to_mOTUs_in-line.tsv
cat $mOTU_folder/db_mOTU/db_mOTU_MAP_MGCs_to_mOTUs.tsv $new_database_folder/$newDBName/updated_min1_db_mOTU/db_mOTU_MAP_MGCs_to_mOTUs.tsv $new_database_folder/$newDBName//$newDBName.mOTU-LG.map.tsv > $new_database_folder/$newDBName/db_mOTU/db_mOTU_MAP_MGCs_to_mOTUs.tsv
cat $mOTU_folder/db_mOTU/db_mOTU_padding_coordinates_CEN.tsv $new_database_folder/$newDBName/updated_min1_db_mOTU/db_mOTU_padding_coordinates_CEN.tsv > $new_database_folder/$newDBName/db_mOTU/db_mOTU_padding_coordinates_CEN.tsv
cat $mOTU_folder/db_mOTU/db_mOTU_padding_coordinates_NR.tsv $new_database_folder/$newDBName/updated_min1_db_mOTU/db_mOTU_padding_coordinates_NR.tsv $new_database_folder/$newDBName/$newDBName.padded.coords > $new_database_folder/$newDBName/db_mOTU/db_mOTU_padding_coordinates_NR.tsv
cat $mOTU_folder/db_mOTU/db_mOTU_taxonomy_meta-mOTUs.tsv $new_database_folder/$newDBName//$newDBName.taxonomy > $new_database_folder/$newDBName/db_mOTU/db_mOTU_taxonomy_meta-mOTUs.tsv




bwa index $new_database_folder/$newDBName/db_mOTU/db_mOTU_DB_CEN.fasta
bwa index $new_database_folder/$newDBName/db_mOTU/db_mOTU_DB_NR.fasta

bwa mem $new_database_folder/$newDBName/db_mOTU/db_mOTU_DB_NR.fasta $mOTU_folder/db_mOTU/db_mOTU_test/test1_single.fastq | grep "^@SQ"  > $new_database_folder/$newDBName/db_mOTU/db_mOTU_bam_header_NR
bwa mem $new_database_folder/$newDBName/db_mOTU/db_mOTU_DB_CEN.fasta $mOTU_folder/db_mOTU/db_mOTU_test/test1_single.fastq | grep "^@SQ"  > $new_database_folder/$newDBName/db_mOTU/db_mOTU_bam_header_CEN
python $scriptDir/create_cami_tax.py $mOTU_folder/db_mOTU/db_mOTU_taxonomy_CAMI.tsv $new_database_folder/$newDBName/db_mOTU/db_mOTU_taxonomy_meta-mOTUs.tsv $new_database_folder/$newDBName/db_mOTU/db_mOTU_taxonomy_CAMI.tsv $newDBName
echo "New database deposited in: $new_database_folder/$newDBName/db_mOTU/"
