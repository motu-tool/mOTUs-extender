#!/usr/bin/env python

import sys
import collections
import argparse
import gc
#import objgraph
#from memory_profiler import profile
# this script will only use weightedmean, but one can use mean if a length file with only 1s is provided

#from guppy import hpy
#h = hpy()

# helper function:
def parseOneColumnFileToList(filename):
    listFileContents = []

    with open(filename, "r") as currFile:
        for strfileLine in currFile:
            strfileLine = strfileLine.strip()
            listFileContents.append(strfileLine)

    return (listFileContents)


def parse2columnFileConnections(infileName):

    dictConnections = collections.defaultdict(lambda: collections.defaultdict(int))

    with open(infileName, "r") as infile:
        for strInLine in infile:
            strInLine = strInLine.strip()
            arrInLine = strInLine.split('\t')
            strColumn_1_entry = arrInLine[0]
            strColumn_2_entry = arrInLine[1]
            dictConnections[strColumn_1_entry][strColumn_2_entry] = 1
            dictConnections[strColumn_2_entry][strColumn_1_entry] = 1

    return (dictConnections)


#column 1:gene name
#column 2:genome
#column 3:OG
#(dictMapping, dictBackMapping_1, dictBackMapping_2) = parse3columnFileMap(infileName)
def parse3columnFileMap(infileName):

    dictMapping = collections.defaultdict(lambda: collections.defaultdict(str))
    dictBackMapping_1 = collections.defaultdict(str)
    dictBackMapping_2 = collections.defaultdict(str)
    
    with open(infileName, "r") as infile:
        for strInLine in infile:
            strInLine = strInLine.strip()
            arrInLine = strInLine.split('\t')
            strColumn_1_entry = arrInLine[0]
            strColumn_2_entry = arrInLine[1]
            strColumn_3_entry = arrInLine[2]
            dictMapping[strColumn_2_entry][strColumn_3_entry] = strColumn_1_entry
            dictBackMapping_1[strColumn_1_entry] = strColumn_2_entry
            dictBackMapping_2[strColumn_1_entry] = strColumn_3_entry
            
    return (dictMapping, dictBackMapping_1, dictBackMapping_2)


def parse3columnFileStr(infileName):
	from collections import defaultdict
	dictIn ={}
	
	infile = open(infileName, "r")

	for strInLine in infile:
		strInLine = strInLine.rstrip('\n')
		arrInLine = strInLine.split('\t')
		
		strColumn_1_entry = arrInLine[0]
		strColumn_2_entry = str(arrInLine[1])
		strColumn_3_entry = str(arrInLine[2])

		if (not strColumn_1_entry in dictIn):
			dictIn[strColumn_1_entry] = (strColumn_2_entry, strColumn_3_entry)
		else:
			sys.stderr.write(strColumn_1_entry + " " + strColumn_2_entry + " assigned multiple times in " + infileName + ".")
		
	return(dictIn)	


#
def parse2columnFileFloat(infileName):
    dictIn = {}

    with open(infileName, "r") as infile:
        for strInLine in infile:
            strInLine = strInLine.strip()
            arrInLine = strInLine.split('\t')

            strColumn_1_entry = arrInLine[0]
            strColumn_2_entry = float(arrInLine[1])

            if (strColumn_1_entry not in dictIn):
                dictIn[strColumn_1_entry] = strColumn_2_entry
            else:
                sys.stderr.write('{} {} assigned multiple times in {}.'.format(strColumn_1_entry, strColumn_2_entry, infileName))

    return (dictIn)


# (distanceMatrixResults, distanceMatrixWeigthSum, dictKey2Position) = initMatrices(listKeys)
def initMatrices(listKeys):
    # numKeys = len(listKeys)
    # distanceMatrixResults = scipy.sparse.dok_matrix((numKeys,numKeys), dtype=numpy.float32)
    # distanceMatrixWeigthSum = scipy.sparse.dok_matrix((numKeys,numKeys), dtype=numpy.float32)

    # distanceMatrixResults = numpy.zeros((numKeys,numKeys), dtype=numpy.float16)
    # distanceMatrixWeigthSum = numpy.zeros((numKeys,numKeys), dtype=numpy.float16)

    distanceMatrixResults = collections.defaultdict(float)
    distanceMatrixWeigthSum = collections.defaultdict(float)
    distanceMatrixOGs = collections.defaultdict(list)
    
    dictKey2Position = collections.defaultdict(int)
    # keyCount = 0
    for i, key in enumerate(listKeys):
        dictKey2Position[key] = i
    # keyCount += 1

    #    keyCount = 0
    #    for i in listKeys:
    #        dictKey2Position[i] = keyCount
    #        keyCount += 1

    return (distanceMatrixResults, distanceMatrixWeigthSum, distanceMatrixOGs, dictKey2Position)

def calculateOGaverageLen(dictGeneID2Length, geneID2OG):

    dictOG_totalLen_helper = collections.defaultdict(int)
    dictOG_totalCount_helper = collections.defaultdict(float)
    
    for geneID in dictGeneID2Length:
        
        gene_length = dictGeneID2Length[geneID]
        OG = geneID2OG[geneID]
        dictOG_totalLen_helper[OG] += int(gene_length)
        dictOG_totalCount_helper[OG] +=  1.0
    
    dictOG_averageLen = collections.defaultdict(float)
    for OG in dictOG_totalLen_helper:
        dictOG_averageLen[OG] = dictOG_totalLen_helper[OG]/dictOG_totalCount_helper[OG]
    
    return(dictOG_averageLen)


def getGenomes2OGList(dictMapGenome2OG2gene):
    
    dictGenomes2OGList = collections.defaultdict(list)
    
    for genome_ID in dictMapGenome2OG2gene:
        for OG in dictMapGenome2OG2gene[genome_ID]:
            dictGenomes2OGList[genome_ID].append(OG)

    return(dictGenomes2OGList)
#####
#####
#####
# real functions
#####
#####

#@profile


def fillDistMatrix(distFileName, dictGene2Genome, dictGene2OG, dictGeneLengths, distanceMatrixResults, distanceMatrixWeigthSum, distanceMatrixOGs, 
                   dictKey2Position, setMatrixEntries, dictRemovedConnections):
                   
    # distFile = open(distFileName, "r")
    # this can be way faster as I read...needs to read 40*15 gb so it should be fast.
    with open(distFileName, "r") as distFile:
        setMatrixEntriesFile = set()
        
        for strLine in distFile:
            # for strLine in distFile:
            strLine = strLine.strip()
            arrLine = strLine.split('\t')

            # print arrLine
            foundPartner = True
            position_1 = -1
            position_2 = -1

            strGeneID_1 = arrLine[0]
            strGeneID_2 = arrLine[1]

            strDist = arrLine[2]
            floatDist = float(strDist)
               
            #XXX|YYY|identifier|taxid
            if (" " in strGeneID_1):
                arrGeneID_1 = strGeneID_1.split(' ')
                strGeneID_1 = arrGeneID_1[0]

            if (" " in strGeneID_2):
                arrGeneID_2 = strGeneID_2.split(' ')
                strGeneID_2 = arrGeneID_2[0]

            genomeID_1 = dictGene2Genome[strGeneID_1]
            genomeID_2 = dictGene2Genome[strGeneID_2]
               
            OG_1 = dictGene2OG[strGeneID_1]
            OG_2 = dictGene2OG[strGeneID_2]
            
            if (OG_1 != OG_2):
                sys.stderr.write(str(strGeneID_1) + " " + str(OG_1) + " " + str(strGeneID_2) + " " + str(OG_2) + " have a distance but are annotated to different OGs.\n")
            
#########
#            if (strGeneID_1.count("|") > 1):
#                arrTaxProjID_1 = strGeneID_1.split('|')
#                strTaxProjID_1 = "|".join(arrTaxProjID_1[0:2] + [arrTaxProjID_1[3]])
#            else:
#                strTaxProjID_1 = strGeneID_1
#
#            if (strGeneID_2.count("|") > 1):
#                arrTaxProjID_2 = strGeneID_2.split('|')
#                strTaxProjID_2 = "|".join(arrTaxProjID_2[0:2] + [arrTaxProjID_2[3]])
#            else:
#                strTaxProjID_2 = strGeneID_2
###########
            floatWeight_1 = 0.0
            if (strGeneID_1 in dictGeneLengths):
                floatWeight_1 = dictGeneLengths[strGeneID_1]

            floatWeight_2 = 0.0
            if (strGeneID_2 in dictGeneLengths):
                floatWeight_2 = dictGeneLengths[strGeneID_2]

            floatWeight = min(floatWeight_1, floatWeight_2)
            floatWeightedDist = floatDist * floatWeight

            if (genomeID_1 in dictKey2Position):
                position_1 = dictKey2Position[genomeID_1]
            else:
                foundPartner = False

            if (genomeID_2 in dictKey2Position):
                position_2 = dictKey2Position[genomeID_2]
            else:
                foundPartner = False
				
            if (strGeneID_1 in dictRemovedConnections and strGeneID_2 in dictRemovedConnections[strGeneID_1] and (dictRemovedConnections[strGeneID_1][strGeneID_2] == 1)):
                foundPartner = False
                sys.stderr.write("Removed distance found in file\n")

            if (position_2 <= position_1):
                position_temp = position_1
                position_1 = position_2
                position_2 = position_temp

            tuplePositions = (position_1, position_2)
            if not (tuplePositions in setMatrixEntriesFile):
                #print (tuplePositions)
                a=1
            #if(True):
                if (foundPartner):
                    # distanceMatrixResults[position_1, position_2] += floatWeightedDist
                    # distanceMatrixWeigthSum[position_1,position_2] += floatWeight
                    distanceMatrixResults[tuplePositions] += floatWeightedDist
                    distanceMatrixWeigthSum[tuplePositions] += floatWeight

                    distanceMatrixOGs[tuplePositions].append(OG_1)
                    setMatrixEntriesFile.add(tuplePositions)
                    setMatrixEntries.add(tuplePositions)

        # else:
        # print(strGeneID_1 + "\t" + strGeneID_2 + "\tfound more than once")
        # print(tuplePositions)
        #print(len(setMatrixEntriesFile))
        del(setMatrixEntriesFile)
        #print(len(dictRemovedConnections))
    #gc.collect()
    #print(len(setMatrixEntries))
    #print(len(distanceMatrixResults))
    #print(len(distanceMatrixWeigthSum))
    #print(len(distanceMatrixCount))
	
    #print(objgraph.show_most_common_types(limit=20))
	
    #print(gc.garbage)
    #return (distanceMatrixResults, distanceMatrixWeigthSum, distanceMatrixCount, distanceMatrixOGs, setMatrixEntries)
    return (distanceMatrixResults, distanceMatrixWeigthSum, distanceMatrixOGs, setMatrixEntries)

def calcResultsMatrixAndPrint(distanceMatrixResults, distanceMatrixWeigthSum, setMatrixEntries, 
                              distanceMatrixOGs, dictGenomes2OGList, dictOG_averageLen, listKeys, countCutoff, defaultDist, outfileName):
#def calcResultsMatrixAndPrint(distanceMatrixResults, distanceMatrixWeigthSum, distanceMatrixCount, setMatrixEntries, 
#                              distanceMatrixOGs, dictGenomes2OGList, dictOG_averageLen, listKeys, countCutoff, defaultDist, outfileName):                              
    with open(outfileName, "w") as outfile:
        setMatrixEntries = sorted(list(setMatrixEntries))
        for tuplePositions in setMatrixEntries:
            #print(tuplePositions)
            (position_1, position_2) = tuplePositions
            id_1 = listKeys[position_1]
            id_2 = listKeys[position_2]
            # weightedDistance = distanceMatrixResults[position_1, position_2]/distanceMatrixWeigthSum[position_1, position_2]
            #weightedDistance = distanceMatrixResults[tuplePositions] / distanceMatrixWeigthSum[tuplePositions]
            
            totalWDist = distanceMatrixResults[tuplePositions]
            totalSum = distanceMatrixWeigthSum[tuplePositions]
            
            listOGswithDist = set(distanceMatrixOGs[tuplePositions])
            numOGswithDist = len(listOGswithDist)
            
            listOGs_1 = dictGenomes2OGList[id_1]
            listOGs_2 = dictGenomes2OGList[id_2]
            setSharedOGs = set(listOGs_1).intersection(set(listOGs_2))
            numSharedOGs = len(setSharedOGs)
            #setSharedOGs - set(listOGswithDist)
            listMissingOGs = list(setSharedOGs.difference(listOGswithDist))
            lenMissingOGs = len(setSharedOGs)
            for OG in listMissingOGs:
                OGlength = dictOG_averageLen[OG]
                OGweight = OGlength * float(defaultDist)
                
                totalWDist += OGweight
                totalSum += OGlength
                
            weightedDistance = totalWDist / totalSum
            
            if (numOGswithDist >= countCutoff):
                strOutLine = '{}\t{}\t{}\t{}\t{}\n'.format(id_1, id_2, weightedDistance, numOGswithDist, numSharedOGs)
                outfile.write(strOutLine)
                #print(strOutLine)
            else:
                strOutLine = '{}\t{}\t{}\t{}\t{}\n'.format(id_1, id_2, defaultDist, numOGswithDist, numSharedOGs)
                outfile.write(strOutLine)
                #print(strOutLine)

def printDicttoErrNice(dictIn):
    
    for key in dictIn:
        value = dictIn[key]
        strOutLine = '{}\t{}\n'.format(key, value)
        sys.stderr.write(strOutLine)


def runFunctions(inFileNameList, geneMapFile, fileIDsList, lengthFileName, removedFileName, countCutoff, defaultDist, outfileName):
    
    listIDs = parseOneColumnFileToList(fileIDsList)

    #column 1:gene name
    #column 2:genome
    #column 3:OG
    (dictGenome2OG2gene, dictGene2Genome, dictGene2OG) = parse3columnFileMap(geneMapFile)
    
    dictGenomes2OGList = getGenomes2OGList(dictGenome2OG2gene)
    
    dictGeneID2Length = parse2columnFileFloat(lengthFileName)
    (distanceMatrixResults, distanceMatrixWeigthSum, distanceMatrixOGs, dictKey2Position) = initMatrices(listIDs)
    setMatrixEntries = set()
    
    dictOG_averageLen = calculateOGaverageLen(dictGeneID2Length, dictGene2OG)
    sys.stderr.write("Average OG lengths calculated\n")
    printDicttoErrNice(dictOG_averageLen)
    
    dictRemovedConnections = parse2columnFileConnections(removedFileName)

    listFileNames = parseOneColumnFileToList(inFileNameList)
    # countCutoff = int(len(listFileNames)/2)
    # countCutoff = int(len(listFileNames)/4)*2


    for distanceFile in listFileNames:
        sys.stderr.write("reading: {} \n".format(distanceFile))
        (distanceMatrixResults, distanceMatrixWeigthSum, distanceMatrixOGs, setMatrixEntries) = fillDistMatrix(
            distanceFile, dictGene2Genome, dictGene2OG, dictGeneID2Length, distanceMatrixResults, distanceMatrixWeigthSum, distanceMatrixOGs, 
            dictKey2Position, setMatrixEntries, dictRemovedConnections)
            
    #calcResultsMatrixAndPrint(distanceMatrixResults, distanceMatrixWeigthSum, distanceMatrixCount, setMatrixEntries, listIDs, countCutoff, outfileName)
    #still needed
    #distanceMatrixOGs, dictGenomes2OGList, defaultDist
    calcResultsMatrixAndPrint(distanceMatrixResults, distanceMatrixWeigthSum, setMatrixEntries, 
                              distanceMatrixOGs, dictGenomes2OGList, dictOG_averageLen, listIDs, countCutoff, defaultDist, outfileName)
    #print h.heap()
###change this from here ... just
def printUsage():
    print("Usage: combineDistances_3.py inFileNameList lengthFileName fileIDsList outfile")
    print(
        "inFileNameList: file with a list of files with distances or similarities in format ID1\tID2\tdist/sim, the files should be combineable i.e. IDs should be in format XXX.YYY.identifier and XXX and YYY will be used to combine them. the whole XXX.YYY.identifier string should be unique so the weights can be assigned")
    print("lengthFileName: file in format ID\tlength; the ID should be XXX.YYY.identifier (see above)")
    print("removedFileName: file in format ID\tID; the ID should be XXX.YYY.identifier (see above)")
    print(
        "fileIDsList: file with a list of ID that should be in the output file. Each ID should be the XXX.YYY in the XXX.YYY.identifier IDs in the single files")
    print("outfile")


def main_old(argv=None):
    argv = sys.argv[1:]

    if (len(argv) == 5):
        inFileNameList = argv[0]
        lengthFileName = argv[1]
        fileIDsList = argv[2]
        removedFileName = argv[3]
        outfileName = argv[4]
    else:
        printUsage()
        sys.exit()

    runFunctions(inFileNameList, fileIDsList, lengthFileName, removedFileName, outfileName)

    return 0  # success


def main(argv=None):
    parser = argparse.ArgumentParser(description='This program calculates mOTU abundances for one sample',
                                     add_help=True)
    parser.add_argument('inFileNameList', action="store",
                        help='File with a list of files with distances or similarities in format ID1\tID2\tdist/sim, which IDs belong together is specified in geneMapFile')
    parser.add_argument('lengthFile', action="store", help='file with lengths of sequences used in blast.')
    parser.add_argument('fileIDsList', action="store", help='file providing a list of genomes/IDs that should be kept')
    parser.add_argument('geneMapFile', action="store", help='file mapping gene name to genome and OG, tab delimited')
    # parser.add_argument('mappingFile', action="store", help='File mapping genes to genome ID of origin')
    parser.add_argument('removedFile', action="store", help='File that provides previously removed distances.')
    parser.add_argument('--countCutoff', '-c', action="store", dest='countCutoff', type=float, default=0.0,
                        help='Minimal number of Marker gene distance pairs needed to compute combined distance.')
    parser.add_argument('--defaultDist', '-d', action="store", dest='defaultDist', type=float, default=0.0,
                        help='Default distance for missing distances, should be set to or just below the cutoff for dist calculations.')
    parser.add_argument('outfileName', action="store", help='Output file.')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()

    #add option default distance (to be set around the dist calculation cutoff, maybe 5% below)
    #add option for gene 2 genome 2 OG file

    # dictGene2Length = parse2columnFile_int(args.lengthfile)
    # countCutoff = 0.0
    # runFunctions(args.inFileNameList, args.fileIDsList, args.lengthFile, args.mappingFile, args.removedFile, args.countCutoff, args.outfileName)
    runFunctions(args.inFileNameList, args.geneMapFile, args.fileIDsList, args.lengthFile, args.removedFile, args.countCutoff, args.defaultDist,
                 args.outfileName)  
    return 0  # success


if __name__ == '__main__':
    status = main()
    sys.exit(status)
