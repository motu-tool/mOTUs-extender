import os
import sys
import argparse
from collections import defaultdict
from collections import Counter

def parseClustering(clusteringFileName, prefix = ""):   
    dictCluster2Members = defaultdict(list) 
    dictID2Cluster = {}
    listClusterNames = []
    
    clusteringFile = open(clusteringFileName, "r")
    
    for strInLine in clusteringFile:
        strInLine = strInLine.rstrip('\n')
        arrInLine = strInLine.split('\t')   
        strClusterName = arrInLine[0]
        
        #cluster name can be modified using a prefix,
        strClusterName = prefix + strClusterName
        
        listClusterNames.append(strClusterName)
        
        strClusterMembers = arrInLine[1]
        arrClusterMembers = strClusterMembers.split(';')
        
        for strClusterMember in arrClusterMembers:
            if (not strClusterMember in dictID2Cluster):
                dictID2Cluster[strClusterMember] = strClusterName
                dictCluster2Members[strClusterName].append(strClusterMember)    
            else:
                print(strClusterMember + " " + strClusterName + " assigned to multiple taxonomy ranks in " + clusteringFileName + ".")          
        
    return(dictCluster2Members, dictID2Cluster, listClusterNames)

def getConsensus(listInput, mode):

    listInput2 = [x for x in listInput if x != "NA NA"]
    counts = Counter(listInput2)
    newTax = "NA NA"
    if(mode == "consensus"):
        if (len(counts) == 1):
            newTax = listInput2[0]
            
    elif(mode == "majority"):
        max = 0
        maxTax = "NA NA"
        for ID in listInput2:
            if(counts[ID] > max):
                maxTax = ID
                max = counts[ID]
        newTax = maxTax
        
    elif(mode == "all"):
        if (len(listInput2)>0):
            newTax = ";".join(listInput2)
        else:
            newTax = "NA NA"
        
    return(newTax)
    
#### make taxonomy file (need script, choices: consensus, majority, all)
#meta-mOTU_v2_ID kingdom phylum  class   order   family  genus   mOTU

def getNotMatchingGenomes(taxonomyFileName, dictCluster2Members):
    
    dictGenome2Taxonomy = defaultdict(list) 
    
    taxonomyFile = open(taxonomyFileName, "r")
    for strLine in taxonomyFile: 
        strLine = strLine.strip()
        arrLine = strLine.split('\t')       
        strGenomeID = arrLine[0]
        listTaxonomy = arrLine[1:]
        
        dictGenome2Taxonomy[strGenomeID] = listTaxonomy
    
    for cluster in dictCluster2Members:
        listClusterMembers = dictCluster2Members[cluster]
        listTaxKingdom = []
        listTaxPhylum = []
        listTaxClass = []
        listTaxOrder = []
        listTaxFamily = []
        listTaxGenus = []
        listTaxSpecies = []
        
        for genome in listClusterMembers:
            #print(genome)
            if (genome in dictGenome2Taxonomy):
                listGenomeTax = dictGenome2Taxonomy[genome]
                listTaxKingdom.append(listGenomeTax[0])
                listTaxPhylum.append(listGenomeTax[1])
                listTaxClass.append(listGenomeTax[2])
                listTaxOrder.append(listGenomeTax[3])
                listTaxFamily.append(listGenomeTax[4])
                listTaxGenus.append(listGenomeTax[5])
                listTaxSpecies.append(listGenomeTax[6])
        taxKingdom = getConsensus(listTaxKingdom, "majority")
        taxPhylum = getConsensus(listTaxPhylum, "majority")
        taxClass = getConsensus(listTaxClass, "majority")
        taxOrder = getConsensus(listTaxOrder, "majority")
        taxFamily = getConsensus(listTaxFamily, "majority")
        taxGenus = getConsensus(listTaxGenus, "majority")
        taxSpecies = getConsensus(listTaxSpecies, "majority")
    
        writeLine = "\t".join([cluster, taxKingdom, taxPhylum, taxClass, taxOrder, taxFamily, taxGenus, taxSpecies])
        print(writeLine)

        
def main(argv=None):
    
    parser = argparse.ArgumentParser(description='This program generates a taxonomy file for a genome clustering from a per genome taxonomy file with columns "genomeID kingdom phylum  class   order   family  genus   species"', add_help = True)
    parser.add_argument('clusteringFileName', action="store", help='file in m8 blast format. Also ok if this is a three column file.')
    parser.add_argument('taxonomyFile', action="store", help='file with lengths of sequences used in blast.')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')      
    args = parser.parse_args()
    
    (dictCluster2Members, dictID2Cluster, listClusterNames) = parseClustering(args.clusteringFileName, prefix = "")
    getNotMatchingGenomes(args.taxonomyFile, dictCluster2Members)
    
    return 0        # success
    
if __name__ == '__main__':
    status = main()
    sys.exit(status)    
                
