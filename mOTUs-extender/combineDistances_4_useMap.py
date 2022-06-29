#!/usr/bin/env python

import os
import sys
import re
import os.path
from collections import defaultdict
import argparse

import numpy

#import scipy.sparse
#import scipy

#this script will only use weightedmean, but one can use mean if a length file with only 1s is provided

#helper function:
def parseOneColumnFileToList(filename):
	listFileContents = []
	
	currFile = open(filename, "r")

	for strfileLine in currFile:
		strfileLine = strfileLine.rstrip('\n')
		listFileContents.append(strfileLine)
		
	return(listFileContents)

def parse2columnFileConnections(infileName):
	from collections import defaultdict
	dictConnections = defaultdict(lambda: defaultdict(int))
	
	infile = open(infileName, "r")

	for strInLine in infile:
		strInLine = strInLine.rstrip('\n')
		arrInLine = strInLine.split('\t')
		
		#print(arrInLine[0])
		#print(arrInLine[1])
		#print("####")
		strColumn_1_entry = arrInLine[0]
		strColumn_2_entry = arrInLine[1]

		dictConnections[strColumn_1_entry][strColumn_2_entry] = 1
		dictConnections[strColumn_2_entry][strColumn_1_entry] = 1 

	return(dictConnections)		
	
	
#	
def parse2columnFileFloat(infileName):
	from collections import defaultdict
	dictIn ={}
	
	infile = open(infileName, "r")

	for strInLine in infile:
		strInLine = strInLine.rstrip('\n')
		arrInLine = strInLine.split('\t')
		
		strColumn_1_entry = arrInLine[0]
		strColumn_2_entry = float(arrInLine[1])

		if (not strColumn_1_entry in dictIn):
			dictIn[strColumn_1_entry] = strColumn_2_entry
		else:
			print(strColumn_1_entry + " assigned multiple times in " + infileName + ".")
		
	return(dictIn)	

def parse2columnFileStr(infileName):
	from collections import defaultdict
	dictIn ={}
	
	infile = open(infileName, "r")

	for strInLine in infile:
		strInLine = strInLine.rstrip('\n')
		arrInLine = strInLine.split('\t')
		
		strColumn_1_entry = arrInLine[0]
		strColumn_2_entry = arrInLine[1]

		if (not strColumn_1_entry in dictIn):
			dictIn[strColumn_1_entry] = strColumn_2_entry
		else:
			print(strColumn_1_entry + " " + strColumn_2_entry + " assigned multiple times in " + infileName + ".")
		
	return(dictIn)		
	

#(distanceMatrixResults, distanceMatrixWeigthSum, dictKey2Position) = initMatrices(listKeys)
def initMatrices(listKeys):
	numKeys = len(listKeys)
	#distanceMatrixResults = scipy.sparse.dok_matrix((numKeys,numKeys), dtype=float32)
	#distanceMatrixWeigthSum = scipy.sparse.dok_matrix((numKeys,numKeys), dtype=float32)
	
	#distanceMatrixResults = numpy.zeros((numKeys,numKeys), dtype=numpy.float16)
	#distanceMatrixWeigthSum = numpy.zeros((numKeys,numKeys), dtype=numpy.float16)

	distanceMatrixResults = defaultdict(float)
	distanceMatrixWeigthSum = defaultdict(float)
	distanceMatrixCount = defaultdict(int)
	
	dictKey2Position = defaultdict(int)
	#keyCount = 0
	for i, key in enumerate(listKeys):
		dictKey2Position[key] = i
		#keyCount += 1

#	keyCount = 0
#	for i in listKeys:
#		dictKey2Position[i] = keyCount
#		keyCount += 1
			
	return(distanceMatrixResults, distanceMatrixWeigthSum, distanceMatrixCount, dictKey2Position)
	
	
#####
#####
#####
#real functions	
#####
#####

#fillDistMatrix(distanceFile, distGeneID2group, dictGeneID2Length, distanceMatrixResults, distanceMatrixWeigthSum, distanceMatrixCount, dictKey2Position, setMatrixEntries, dictRemovedConnections)
def fillDistMatrix(distFileName, distGeneID2group, dictGeneLengths, distanceMatrixResults, distanceMatrixWeigthSum, distanceMatrixCount, dictKey2Position, setMatrixEntries, dictRemovedConnections):
	setMatrixEntriesFile = set()

	distFile = open(distFileName, "r")
	
#this can be way faster as I read...needs to read 40*15 gb so it should be fast.
	for strLine in distFile:
		strLine = strLine.strip()
		arrLine = strLine.split('\t')
		
		#print arrLine
		groupFound = True
		foundPartner = True
		position_1 = -1
		position_2 = -1
		
		strGeneID_1 = arrLine[0]
		strGeneID_2 = arrLine[1]	

		strDist = arrLine[2]
		floatDist = float(strDist)

		if (" " in strGeneID_1):
			arrGeneID_1 = strGeneID_1.split(' ')
			strGeneID_1 = arrGeneID_1[0]
		if ("\t" in strGeneID_1):
			arrGeneID_1 = strGeneID_1.split('\t')
			strGeneID_1 = arrGeneID_1[0]
			
		if (" " in strGeneID_2):
			arrGeneID_2 = strGeneID_2.split(' ')
			strGeneID_2 = arrGeneID_2[0]
		if ("\t" in strGeneID_2):
			arrGeneID_2 = strGeneID_2.split('\t')
			strGeneID_2 = arrGeneID_2[0]

		if (strGeneID_1 in distGeneID2group):
			strTaxProjID_1 = distGeneID2group[strGeneID_1]
		else:
			groupFound = False
			strTaxProjID_1 = "None"
			
		if (strGeneID_2 in distGeneID2group):
			strTaxProjID_2 = distGeneID2group[strGeneID_2]	
		else:
			groupFound = False
			strTaxProjID_2 = "None"
			

		floatWeight_1 = 0.0
		if (strGeneID_1 in dictGeneLengths):
			floatWeight_1 = dictGeneLengths[strGeneID_1]
		
		floatWeight_2 = 0.0
		if (strGeneID_2 in dictGeneLengths):
			floatWeight_2 = dictGeneLengths[strGeneID_2]

		floatWeight = min(floatWeight_1, floatWeight_2)
		floatWeightedDist = floatDist*floatWeight

		if (strTaxProjID_1 in dictKey2Position):
			position_1 = dictKey2Position[strTaxProjID_1]
		else:
			foundPartner = False

		if (strTaxProjID_2 in dictKey2Position):
			position_2 = dictKey2Position[strTaxProjID_2]
		else:
			foundPartner = False
		if (strGeneID_1 in dictRemovedConnections and strGeneID_2 in dictRemovedConnections[strGeneID_1] and (dictRemovedConnections[strGeneID_1][strGeneID_2] == 1)):
			foundPartner = False
		
		if (position_2 <= position_1):
			position_temp = position_1
			position_1 = position_2
			position_2 = position_temp
		
		tuplePositions = (position_1, position_2)
		if not (tuplePositions in setMatrixEntriesFile):
			if (foundPartner):
				#distanceMatrixResults[position_1, position_2] += floatWeightedDist
				#distanceMatrixWeigthSum[position_1,position_2] += floatWeight
				distanceMatrixResults[tuplePositions] += floatWeightedDist
				distanceMatrixWeigthSum[tuplePositions] += floatWeight
				distanceMatrixCount[tuplePositions] += 1
				setMatrixEntriesFile.add(tuplePositions)
				setMatrixEntries.add(tuplePositions)
		#else:
			#print(strGeneID_1 + "\t" + strGeneID_2 + "\tfound more than once")
			#print(tuplePositions)
			
	return(distanceMatrixResults, distanceMatrixWeigthSum, distanceMatrixCount, setMatrixEntries)
	
	
def calcResultsMatrixAndPrint(distanceMatrixResults, distanceMatrixWeigthSum, distanceMatrixCount, setMatrixEntries, listKeys, countCutoff, outfileName):

	outfile = open(outfileName, "w")
	
	setMatrixEntries = sorted(list(setMatrixEntries))

	for tuplePositions in setMatrixEntries:
		(position_1, position_2) = tuplePositions
		ID_1 = listKeys[position_1]
		ID_2 = listKeys[position_2]

		#weightedDistance = distanceMatrixResults[position_1, position_2]/distanceMatrixWeigthSum[position_1, position_2]
		weightedDistance = distanceMatrixResults[tuplePositions]/distanceMatrixWeigthSum[tuplePositions]
		if (distanceMatrixCount[tuplePositions] >= countCutoff):
			strOutLine = str(ID_1) + "\t" + str(ID_2) + "\t" + str(weightedDistance) + "\n"
			outfile.write(strOutLine)


def runFunctions(inFileNameList, fileIDsList, lengthFileName, mappingFileName, removedFileName, countCutoff, outfileName):

	listGroupIDs = parseOneColumnFileToList(fileIDsList)
	distGeneID2group = parse2columnFileStr(mappingFileName)
	dictGeneID2Length = parse2columnFileFloat(lengthFileName)
	(distanceMatrixResults, distanceMatrixWeigthSum, distanceMatrixCount, dictKey2Position) = initMatrices(listGroupIDs)
	setMatrixEntries = set()
	
	dictRemovedConnections = parse2columnFileConnections(removedFileName)
	
	listFileNames = parseOneColumnFileToList(inFileNameList)
	#countCutoff = int(len(listFileNames)/2)
	#countCutoff = int(len(listFileNames)/4)*2
	
	for distanceFile in listFileNames:
		print("reading: " + distanceFile)
		(distanceMatrixResults, distanceMatrixWeigthSum, distanceMatrixCount, setMatrixEntries) = fillDistMatrix(distanceFile, distGeneID2group, dictGeneID2Length, distanceMatrixResults, distanceMatrixWeigthSum, distanceMatrixCount, dictKey2Position, setMatrixEntries, dictRemovedConnections)
		#print(len(setMatrixEntries))

	calcResultsMatrixAndPrint(distanceMatrixResults, distanceMatrixWeigthSum, distanceMatrixCount, setMatrixEntries, listGroupIDs, countCutoff, outfileName)
		

###change this from here ... just 
def printUsage():
	print("Usage: combineDistances_3.py inFileNameList lengthFileName fileIDsList outfile")
	print("inFileNameList: file with a list of files with distances or similarities in format ID1\tID2\tdist/sim, the files should be combineable i.e. IDs should be in format XXX.YYY.identifier and XXX and YYY will be used to combine them. the whole XXX.YYY.identifier string should be unique so the weights can be assigned")
	print("lengthFileName: file in format ID\tlength; the ID should be XXX.YYY.identifier (see above)")
	print("removedFileName: file in format ID\tID; the ID should be XXX.YYY.identifier (see above)")
	print("fileIDsList: file with a list of ID that should be in the output file. Each ID should be the XXX.YYY in the XXX.YYY.identifier IDs in the single files")
	print("outfile")
		

def main(argv=None):
	argv = sys.argv[1:]
	
	if (len(argv) == 5):
		inFileNameList = argv[0]
		lengthFileName = argv[1]
		fileIDsList = argv[2]
		removedFileName = argv[3]
		mappingFileName = argv[4]
		outfileName = argv[5]
	else:
		printUsage()
		printUsage()
		sys.exit()
	
	countCutoff = 0
	#runFunctions(inFileNameList, fileIDsList, lengthFileName, removedFileName, outfileName)
	runFunctions(inFileNameList, fileIDsList, lengthFileName, mappingFileName, removedFileName, countCutoff, outfileName)
	return 0        # success

	

def main(argv=None):
	
	parser = argparse.ArgumentParser(description='This program calculates mOTU abundances for one sample', add_help = True)
	parser.add_argument('blastfile', action="store", help='file in m8 blast format.')
	parser.add_argument('lengthFile', action="store", help='file with lengths of sequences used in blast.')
	parser.add_argument('fileIDsList', action="store", help='file providing a list of genomes/IDs that should be kept')
	parser.add_argument('mappingFile', action="store", help='File mapping genes to genome ID of origin')
	parser.add_argument('removedFile', action="store", help='File that provides previously removed distances.')
	parser.add_argument('outfileName', action="store", help='Output file.')
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')		
	args = parser.parse_args()
	
	#dictGene2Length = parse2columnFile_int(args.lengthfile)
	countCutoff = 0.0
	runFunctions(args.blastfile, args.fileIDsList, args.lengthFile, args.mappingFile, args.removedFile, countCutoff, args.outfileName)
	return 0        # success
	
if __name__ == '__main__':
	status = main()
	sys.exit(status)	
				