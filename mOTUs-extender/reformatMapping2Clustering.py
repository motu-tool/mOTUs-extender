import os
import sys

from collections import defaultdict


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


def reformatClustering2Mapping(dictCluster2Members, outfileName):

	outfile = open(outfileName, "w")
	
	for cluster in dictCluster2Members:
		listMembers = dictCluster2Members[cluster]
		strMembers = ';'.join(listMembers)
		
		strWriteLine = cluster + "\t" + strMembers
		outfile.write(strWriteLine + "\n")

########################	
def printUsage():
	print("Usage: reformatClustering2Mapping_2.py infile outfile")
	print("infile: clusteringfile")
	print("outfile: name of the outfile")


def main(argv=None):
	if(not argv):
		argv = sys.argv[1:]

	if (len(argv) == 2):
		infileName = argv[0]
		outfileName = argv[1]
	else:
		printUsage()
		sys.exit()
	
	(dictCluster2Members, dictID2Cluster, listClusterNames) = parseClustering(infileName)
	reformatClustering2Mapping(dictCluster2Members, outfileName)
	
	return 0        # success
	
if __name__ == '__main__':
    status = main()
    sys.exit(status)		
