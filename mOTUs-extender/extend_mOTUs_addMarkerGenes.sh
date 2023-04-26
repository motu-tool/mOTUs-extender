#set -euxo pipefail
#add genomes to mOTU DB

sequence_file=$1
genome_ID=$2
new_database_folder=$3
scriptDir=$4
mOTU_folder=$5
fetchMGs_folder=$6
mapfile=$7
cutoffsFile="$scriptDir/cutoffs_fscore_specIAsRef.csv"
mOTU_MG_file="$scriptDir/10RefMGs.IDs"


db_mOTU_genes_length_NR=$mOTU_folder/db_mOTU/db_mOTU_genes_length_NR
unpadded_mOTUseqsFile="$mOTU_folder/db_mOTU/db_mOTU_DB_NR.fasta.udb"




########
# run program

# TODO allow to add also MGs and Genome instead of only Genome. That way we can avoid the 2 fetchmgs + 1 prodigal call
mkdir -p $new_database_folder/dbs/$genome_ID
#python extend_mOTUs_DB/SCRIPTS/load_and_pad_10MGs.py  temp_db_folder/db_mOTU/db_mOTU_DB_NR.fasta.udb

python $scriptDir/load_and_pad_10MGs.py $sequence_file $new_database_folder/dbs/$genome_ID $fetchMGs_folder $genome_ID $unpadded_mOTUseqsFile


#------------------------
#make global len file
cat $new_database_folder/dbs/${genome_ID}/${genome_ID}.genes.len $db_mOTU_genes_length_NR > $new_database_folder/dbs/${genome_ID}/${genome_ID}.genes.all.len

#make global map file

cat $new_database_folder/dbs/${genome_ID}/${genome_ID}.map2genome $mapfile > $new_database_folder/dbs/${genome_ID}/${genome_ID}.all.map2genome
python $scriptDir/map2genome.py $new_database_folder/dbs/${genome_ID}/${genome_ID}.all.map2genome $new_database_folder/dbs/${genome_ID}/${genome_ID}.all.map2genome_v5

#make global genome IDs file
cut -f2 $new_database_folder/dbs/${genome_ID}/${genome_ID}.map2genome > $new_database_folder/dbs/$genome_ID/${genome_ID}.genomeIDs
cut -f2 $new_database_folder/dbs/${genome_ID}/${genome_ID}.all.map2genome > $new_database_folder/dbs/$genome_ID/${genome_ID}.all.genomeIDs



mkdir -p $new_database_folder/dbs/$genome_ID/vsearch/normalized/
while read line;
do
IFS=","; set $line
python $scriptDir/normalizePercentIdentity.py $new_database_folder/dbs/$genome_ID/vsearch/$1.distances_vs_db.m8 $2 > $new_database_folder/dbs/$genome_ID/vsearch/normalized/$1.distances_vs_db.m8
done < $cutoffsFile


while read line;
do
python $scriptDir/filterBlastByAlignLength_Fraction.py $new_database_folder/dbs/$genome_ID/vsearch/normalized/$line.distances_vs_db.m8 $new_database_folder/dbs/${genome_ID}/${genome_ID}.genes.all.len -l 20 -f 0.50 -r $new_database_folder/dbs/$genome_ID/vsearch/normalized/$line.distances_vs_db.excludedPairs > $new_database_folder/dbs/$genome_ID/vsearch/normalized/$line.distances_vs_db.50p.20n.m8 
done < $mOTU_MG_file

cat $new_database_folder/dbs/$genome_ID/vsearch/normalized//COG*.distances_vs_db.excludedPairs > $new_database_folder/dbs/${genome_ID}/vsearch/normalized/AllCOGs.distances_vs_db.excludedPairs

ls $new_database_folder/dbs/${genome_ID}/vsearch/normalized/COG*.distances_vs_db.50p.20n.m8  > $new_database_folder/dbs/${genome_ID}/vsearch/normalized//files_normalized.txt

python $scriptDir/combineDistances_5.py -d 55.0 -c 3 $new_database_folder/dbs/${genome_ID}/vsearch/normalized//files_normalized.txt $new_database_folder/dbs/${genome_ID}/${genome_ID}.genes.all.len $new_database_folder/dbs/${genome_ID}/${genome_ID}.all.genomeIDs $new_database_folder/dbs/${genome_ID}/${genome_ID}.all.map2genome_v5 $new_database_folder/dbs/${genome_ID}/vsearch/normalized/AllCOGs.distances_vs_db.excludedPairs $new_database_folder/dbs/${genome_ID}/vsearch/combined.normalized.distances_vs_db.m8

#python $scriptDir/combineDistances_4_useMap.py $new_database_folder/dbs/${genome_ID}/vsearch/normalized//files_normalized.txt $new_database_folder/dbs/${genome_ID}/${genome_ID}.genes.all.len $new_database_folder/dbs/${genome_ID}/${genome_ID}.all.genomeIDs $new_database_folder/dbs/${genome_ID}/${genome_ID}.all.map2genome $new_database_folder/dbs/${genome_ID}/vsearch/normalized/AllCOGs.distances_vs_db.excludedPairs $new_database_folder/dbs/${genome_ID}/vsearch/combined.normalized.distances_vs_db.m8



