import os
import sys
import argparse
from collections import defaultdict

def parseOneColumnFileToList(filename):
	listFileContents = []
	
	currFile = open(filename, "r")

	for strfileLine in currFile:
		strfileLine = strfileLine.rstrip('\n')
		listFileContents.append(strfileLine)
		
	return(listFileContents)

def getNotMatchingGenomes(blastFileName, IDsList, cutoff):
	
	IDsListSet = set(IDsList)

	blastFile = open(blastFileName, 'r')
	for strLine in blastFile:
		strLine = strLine.strip()
		strLine2 = str(strLine)
		arrLine = strLine.split('\t')
		
		strQuery = arrLine[0]
		strSubject = arrLine[1]
		#print(arrLine)
		floatPercentID = float(arrLine[2])
		
		if (floatPercentID >= cutoff):
			if (strQuery in IDsListSet):
				IDsListSet.remove(strQuery)
			if (strSubject in IDsListSet):
				IDsListSet.remove(strSubject)
				
	IDsList = list(IDsListSet)
	for ID in IDsList:
		print(ID)
		

def main(argv=None):
	
	parser = argparse.ArgumentParser(description='This program calculates mOTU abundances for one sample', add_help = True)
	parser.add_argument('blastfile', action="store", help='file in m8 blast format. Also ok if this is a three column file.')
	parser.add_argument('genomeIDsFileName', action="store", help='file with lengths of sequences used in blast.')
	parser.add_argument('--cutoff', '-c', action="store", dest='cutoff', type=float, default=99.0, help='maximal percent ID for genomes to be kept, 99% used for normalized percent IDs.')
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')		
	args = parser.parse_args()
	
	genomeIDs = parseOneColumnFileToList(args.genomeIDsFileName)
	getNotMatchingGenomes(args.blastfile, genomeIDs, args.cutoff)

	return 0        # success
	
if __name__ == '__main__':
	status = main()
	sys.exit(status)	
				