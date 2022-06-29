import os
import sys
import argparse
from collections import defaultdict

def parse2columnFile_int(infileName):
	dictIn = defaultdict(int)
	
	infile = open(infileName, "r")

	for strLine in infile:
	
		strLine = strLine.rstrip('\n')
		arrLine = strLine.split('\t')
		
		strColumn_1_entry = arrLine[0]
		strColumn_2_entry = int(arrLine[1])

		dictIn[strColumn_1_entry] = strColumn_2_entry
		
	return(dictIn)
	

def filterAlignLength_Fraction(blastFileName, intMinAlignLength, dictGene2Length, minAlignFraction, removedDistFileName):

	blastFile = open(blastFileName, 'r')
	writeOutRemoved = False
	if (removedDistFileName != ""):
		removedDistFile = open(removedDistFileName, 'w')
		writeOutRemoved = True
		
	for strLine in blastFile:
		strLine = strLine.strip()
		strLine2 = str(strLine)
		arrLine = strLine.split('\t')
		
		strQuery = arrLine[0]
		strSubject = arrLine[1]
		#print(arrLine)
		intAlignLength = int(arrLine[3])
		
		queryLen = dictGene2Length[strQuery]
		subjectLen = dictGene2Length[strSubject]
		
		minSeqLength = min(queryLen, subjectLen)
		alignFraction = float(intAlignLength)/float(minSeqLength)
		
		if (intAlignLength > intMinAlignLength and alignFraction > minAlignFraction):
			print(strLine)
		elif(writeOutRemoved):
			removedLine = strQuery + "\t" + strSubject
			removedDistFile.write(removedLine + "\n")
		
		
		

def main(argv=None):
	
	parser = argparse.ArgumentParser(description='This program calculates mOTU abundances for one sample', add_help = True)
	parser.add_argument('blastfile', action="store", help='file in m8 blast format.')
	parser.add_argument('lengthfile', action="store", help='file with lengths of sequences used in blast.')
	parser.add_argument('--minAlignLength', '-l', action="store", dest='intMinAlignLength', type=int, default=100, help='minimal length of the alignment as found in column 4 in the m8 format.')
	parser.add_argument('--minAlignFraction', '-f', action="store", dest='minAlignFraction', type=float, default=0.5, help='minimal fraction overlapping between the 2 sequences (as fraction of the shorter sequence).')
	parser.add_argument('--removedDistFile', '-r', action="store", dest='removedDistFile', type=str, default="", help='File to store removed distances.')
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')		
	args = parser.parse_args()
	
	dictGene2Length = parse2columnFile_int(args.lengthfile)
	
	filterAlignLength_Fraction(args.blastfile, args.intMinAlignLength, dictGene2Length, args.minAlignFraction, args.removedDistFile)

	return 0        # success
	
if __name__ == '__main__':
	status = main()
	sys.exit(status)	
				