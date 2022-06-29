import os
import sys

def reformatClustering2Mapping(infileName, outfilename):

	infile = open(infileName, "r")
	outfile = open(outfilename, "w")
	
	for strInLine in infile:
		strInLine = strInLine.rstrip('\n')
		arrInLine = strInLine.split('\t')
		
		strClusterNumber = arrInLine[0]
		strClusters = arrInLine[1]
		arrClusters = strClusters.split(';')	
		
		for strGeneID in arrClusters:
			strWriteLine = strClusterNumber + "\t" + strGeneID
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
	
	reformatClustering2Mapping(infileName, outfileName)
	
	return 0        # success
	
if __name__ == '__main__':
    status = main()
    sys.exit(status)		
