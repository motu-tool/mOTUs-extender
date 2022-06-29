set -euxo pipefail
#add genomes to mOTU DB

sequence_file=$1
genome_ID=$2
new_database_folder=$3
scriptDir=$4
mOTU_folder=$5
cutoffsFile="$scriptDir/cutoffs_fscore_specIAsRef.csv"
mOTU_MG_file="$scriptDir/10RefMGs.IDs"


db_mOTU_genes_length_NR=$mOTU_folder/db_mOTU/db_mOTU_genes_length_NR
unpadded_mOTUseqsFile="$mOTU_folder/db_mOTU/db_mOTU_DB_NR.fasta"



#if [ ! -f "${unpadded_mOTUseqsFile}.udb" ]; then
#    echo "${unpadded_mOTUseqsFile}.udb not found. Indexing database"
#	vsearch --makeudb_usearch ${unpadded_mOTUseqsFile} --output ${unpadded_mOTUseqsFile}.udb
#fi

########
# run program

mkdir -p $new_database_folder/dbs/
python2 $scriptDir/makePaddedMGdb.py $sequence_file $new_database_folder/dbs/$genome_ID $genome_ID -f $scriptDir/fetchMGs/ -m single --OGType mOTU

#fetchMGs.pl -m extraction -t 1 -o $new_database_folder/dbs/${genome_ID}/${genome_ID}_markerGenes2 -d $new_database_folder/dbs/${genome_ID}/${genome_ID}.all.genes.fasta $new_database_folder/dbs/${genome_ID}/${genome_ID}.all.proteins.fasta
$scriptDir/fetchMGs/fetchMGs.pl -i -v -m extraction -t 1 -o $new_database_folder/dbs/${genome_ID}/${genome_ID}_markerGenes2 -d $new_database_folder/dbs/${genome_ID}/${genome_ID}.all.genes.fasta $new_database_folder/dbs/${genome_ID}/${genome_ID}.all.proteins.fasta



while read line;
do
sed -i "s/^>.*$/>${genome_ID}.$line/g" $new_database_folder/dbs/${genome_ID}/${genome_ID}_markerGenes2/$line.fna 
sed -i "s/^>.*$/>${genome_ID}.$line/g" $new_database_folder/dbs/${genome_ID}/${genome_ID}_markerGenes2/$line.faa
done < $mOTU_MG_file


while read line;
do
sed 's/\t/ /g' $new_database_folder/dbs/${genome_ID}/${genome_ID}_markerGenes2/$line.fna > $new_database_folder/dbs/${genome_ID}/${genome_ID}_markerGenes2/$line.notab.fna
done < $mOTU_MG_file

while read line;
do
seqtk comp $new_database_folder/dbs/${genome_ID}/${genome_ID}_markerGenes2/$line.notab.fna | cut -f1
done < $mOTU_MG_file | sed "s/$/\t${genome_ID}/" > $new_database_folder/dbs/${genome_ID}/${genome_ID}.map2genome

while read line;
do
seqtk comp $new_database_folder/dbs/${genome_ID}/${genome_ID}_markerGenes2/$line.notab.fna | cut -f1,2 | sed "s/$/\t${line}\t${genome_ID}/"
done < $mOTU_MG_file > $new_database_folder/dbs/${genome_ID}/${genome_ID}.map

#make global len file
cat $new_database_folder/dbs/${genome_ID}/${genome_ID}.genes.len $db_mOTU_genes_length_NR > $new_database_folder/dbs/${genome_ID}/${genome_ID}.genes.all.len

#make global map file
cat $new_database_folder/dbs/${genome_ID}/${genome_ID}.map2genome $scriptDir/mOTU.v2.5b.mOTU-LG.map > $new_database_folder/dbs/${genome_ID}/${genome_ID}.all.map2genome

#make global genome IDs file
cut -f2 $new_database_folder/dbs/${genome_ID}/${genome_ID}.map2genome > $new_database_folder/dbs/$genome_ID/${genome_ID}.genomeIDs
cut -f2 $new_database_folder/dbs/${genome_ID}/${genome_ID}.all.map2genome > $new_database_folder/dbs/$genome_ID/${genome_ID}.all.genomeIDs

mkdir -p $new_database_folder/dbs/$genome_ID/vsearch
while read line;
do
if [ -s "$new_database_folder/dbs/${genome_ID}/${genome_ID}_markerGenes2/$line.notab.fna" ]
then
vsearch --threads 1 --usearch_global $new_database_folder/dbs/${genome_ID}/${genome_ID}_markerGenes2/$line.notab.fna --db ${unpadded_mOTUseqsFile}.udb --id 0.8  --maxaccepts 100000 --mincols 20 --blast6out $new_database_folder/dbs/$genome_ID/vsearch/$line.distances_vs_db.m8
else
touch $new_database_folder/dbs/$genome_ID/vsearch/$line.distances_vs_db.m8
fi
done < $mOTU_MG_file

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

python $scriptDir/combineDistances_4_useMap.py $new_database_folder/dbs/${genome_ID}/vsearch/normalized//files_normalized.txt $new_database_folder/dbs/${genome_ID}/${genome_ID}.genes.all.len $new_database_folder/dbs/${genome_ID}/${genome_ID}.all.genomeIDs $new_database_folder/dbs/${genome_ID}/${genome_ID}.all.map2genome $new_database_folder/dbs/${genome_ID}/vsearch/normalized/AllCOGs.distances_vs_db.excludedPairs $new_database_folder/dbs/${genome_ID}/vsearch/combined.normalized.distances_vs_db.m8



