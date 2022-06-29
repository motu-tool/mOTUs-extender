import os
import sys
import argparse
from collections import defaultdict


def normalizePercentIdentity(blastFileName, percentID_NormFactor):

	distanceNormFactor = 100.0 - float(percentID_NormFactor)
	blastFile = open(blastFileName, 'r')
		
	for strLine in blastFile:
		strLine = strLine.strip()
		arrLine = strLine.split('\t')
		
		strQuery = arrLine[0]
		strSubject = arrLine[1]
		#print(arrLine)
		
		floatPercID = float(arrLine[2])
		restOfLine = ""
		if(len(arrLine) > 3):
			restOfLine = "\t" + "\t".join(arrLine[3:])
		
		floatDistance = 100.0 - floatPercID
		floatDistance_norm = floatDistance/distanceNormFactor
		#print(floatDistance)
		#print(distanceNormFactor)
		#print(floatDistance_norm)
		floatPercID_norm = 100.0 - floatDistance_norm
		if (floatPercID_norm < 0):
			floatPercID_norm = 0
		
		strWriteLine = strQuery + "\t" + strSubject + "\t" + str(floatPercID_norm) + restOfLine
		print(strWriteLine)
		
#100 = 100
#cutoff = 0.95		
		

def main(argv=None):
	
	parser = argparse.ArgumentParser(description='This program calculates mOTU abundances for one sample', add_help = True)
	parser.add_argument('blastfile', action="store", help='file in m8 blast format.')
	parser.add_argument('percentID_NormFactor', action="store", type=float, default=95.0, help='normaliation factor that is applied to the distances/percentIDs')
	#parser.add_argument('--minAlignLength', '-l', action="store", dest='intMinAlignLength', type=int, default=100, help='minimal length of the alignment as found in column 4 in the m8 format.')

	parser.add_argument('--version', action='version', version='%(prog)s 1.0')		
	args = parser.parse_args()
	
	normalizePercentIdentity(args.blastfile, args.percentID_NormFactor)

	return 0        # success
	
if __name__ == '__main__':
	status = main()
	sys.exit(status)	
				