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

    dictConnections = set()

    with open(infileName, "r") as infile:
        for strInLine in infile:
            arrInLine = strInLine.strip().split('\t')
            strColumn_1_entry = arrInLine[0]
            strColumn_2_entry = arrInLine[1]
            dictConnections.add((strColumn_2_entry, strColumn_1_entry))
            dictConnections.add((strColumn_1_entry, strColumn_2_entry))

    return dictConnections


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
            
    return (dictMapping, dictBackMapping_1, dictBackMapping_2) # genome_2_geneIds, geneId_2_genome, geneId_2_cog


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
def parse2columnFileInt(infileName):
    dictIn = {}

    with open(infileName, "r") as infile:
        for strInLine in infile:
            strInLine = strInLine.strip()
            arrInLine = strInLine.split('\t')

            strColumn_1_entry = arrInLine[0]
            strColumn_2_entry = int(arrInLine[1])

            if (strColumn_1_entry not in dictIn):
                dictIn[strColumn_1_entry] = strColumn_2_entry
            else:
                sys.stderr.write('{} {} assigned multiple times in {}.'.format(strColumn_1_entry, strColumn_2_entry, infileName))

    return (dictIn)

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


# # (distanceMatrixResults, distanceMatrixWeigthSum, dictKey2Position) = initMatrices(listKeys)
# def initMatrices(listKeys):
#     # numKeys = len(listKeys)
#     # distanceMatrixResults = scipy.sparse.dok_matrix((numKeys,numKeys), dtype=numpy.float32)
#     # distanceMatrixWeigthSum = scipy.sparse.dok_matrix((numKeys,numKeys), dtype=numpy.float32)
#
#     # distanceMatrixResults = numpy.zeros((numKeys,numKeys), dtype=numpy.float16)
#     # distanceMatrixWeigthSum = numpy.zeros((numKeys,numKeys), dtype=numpy.float16)
#
#     distanceMatrixResults = collections.defaultdict(float)
#     distanceMatrixWeigthSum = collections.defaultdict(float)
#     distanceMatrixOGs = collections.defaultdict(list)
#
#     dictKey2Position = collections.defaultdict(int)
#     # keyCount = 0
#     for i, key in enumerate(listKeys):
#         dictKey2Position[key] = i
#     # keyCount += 1
#
#     #    keyCount = 0
#     #    for i in listKeys:
#     #        dictKey2Position[i] = keyCount
#     #        keyCount += 1
#
#     return (distanceMatrixResults, distanceMatrixWeigthSum, distanceMatrixOGs, dictKey2Position)

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
        dictOG_averageLen[OG] = int(dictOG_totalLen_helper[OG]/dictOG_totalCount_helper[OG])
    
    return(dictOG_averageLen)


def getGenomes2OGSet(dictMapGenome2OG2gene):
    
    dictGenomes2OGSet = collections.defaultdict(set)
    
    for genome_ID in dictMapGenome2OG2gene:
        for OG in dictMapGenome2OG2gene[genome_ID]:
            dictGenomes2OGSet[genome_ID].add(OG)

    return(dictGenomes2OGSet)
#####
#####
#####
# real functions
#####
#####

#@profile


def fillDistMatrix(distFileName, gene_2_genome, gene_2_cog, gene_2_length, dictRemovedConnections, alignments):
                   

    with open(distFileName, "r") as distFile:
        for line in distFile:

            splits = line.strip().split('\t', 3)


            geneId_1 = splits[0]
            geneId_2 = splits[1]
            intDist = int(float(splits[2]) * 100)

            genomeID_1 = gene_2_genome[geneId_1]
            genomeID_2 = gene_2_genome[geneId_2]
               
            cog_1 = gene_2_cog[geneId_1]
            cog_2 = gene_2_cog[geneId_2]
            
            if (cog_1 != cog_2):
                sys.stderr.write(str(geneId_1) + " " + str(genomeID_1) + " " + str(geneId_2) + " " + str(genomeID_2) + " have a distance but are annotated to different OGs.\n")
            

            gene_length_1 = gene_2_length[geneId_1] # The weight is basically the length
            gene_length_2 = gene_2_length[geneId_2]


            shorter_gene_length = min(gene_length_1, gene_length_2) # pick the shorter gene
            #weighted_distance = intDist * shorter_gene_length # length time percid

            # replace with a tuple/set search
            if ((geneId_1, geneId_2) in dictRemovedConnections):
                sys.stderr.write("Removed distance found in file\n")
            else:
                if genomeID_1 < genomeID_2:
                    alignments[(genomeID_1, genomeID_2)][cog_1] = (intDist, shorter_gene_length)
                else:
                    alignments[(genomeID_2, genomeID_1)][cog_1]  = (intDist, shorter_gene_length)
    return alignments


def get_weighted_sum_and_total_length(distances):
    total_length = 0
    total_weighted_sum = 0
    for _, (distance, length) in distances.items():
        total_length += length
        total_weighted_sum += length * distance
    return total_length, total_weighted_sum
def calcResultsMatrixAndPrint(alignments, outfileName, countCutoff, defaultDist, genome_2_set_of_cogs, cog_2_avglength):
    with open(outfileName, "w") as outfile:

        for (genomeId_1, genomeId_2), distances in alignments.items():
            if len(distances) < countCutoff: # not enough aligned
                continue
            genomeId_1_cogs = genome_2_set_of_cogs[genomeId_1]
            genomeId_2_cogs = genome_2_set_of_cogs[genomeId_2]
            total_cog_overlap = genomeId_1_cogs.intersection(genomeId_2_cogs)
            if len(total_cog_overlap) < countCutoff:
                continue
            totalSum, totalWDist = get_weighted_sum_and_total_length(distances)
            if len(distances) != len(total_cog_overlap):
                for cog in filter(lambda x: x not in distances, total_cog_overlap):
                    OGlength = cog_2_avglength[cog]
                    OGweight = defaultDist * OGlength
                    totalWDist += OGweight
                    totalSum += OGlength
            weightedDistance = float(totalWDist / totalSum)/ 100.0
            if weightedDistance < 60.0:
                continue
            weightedDistance = "{:.2f}".format(weightedDistance)
            strOutLine = f'{genomeId_1}\t{genomeId_2}\t{weightedDistance}\n'
            outfile.write(strOutLine)


def printDicttoErrNice(dictIn):
    
    for key in dictIn:
        value = dictIn[key]
        strOutLine = '{}\t{}\n'.format(key, int(value))
        sys.stderr.write(strOutLine)


def runFunctions(inFileNameList, geneMapFile, fileIDsList, lengthFileName, removedFileName, countCutoff, defaultDist, outfileName):
    
    genomeIds = parseOneColumnFileToList(fileIDsList)

    (genome_2_geneIds, geneId_2_genome, geneId_2_cog) = parse3columnFileMap(geneMapFile)
    
    genome_2_set_of_cogs = getGenomes2OGSet(genome_2_geneIds)

    geneId_2_length = parse2columnFileInt(lengthFileName)

    
    dictOG_averageLen = calculateOGaverageLen(geneId_2_length, geneId_2_cog)
    sys.stderr.write("Average OG lengths calculated\n")
    printDicttoErrNice(dictOG_averageLen)
    
    dictRemovedConnections = parse2columnFileConnections(removedFileName)

    listFileNames = parseOneColumnFileToList(inFileNameList)
    alignments = collections.defaultdict(lambda: collections.defaultdict(int))
    for distanceFile in listFileNames:

        sys.stderr.write("reading: {} \n".format(distanceFile))
        fillDistMatrix(distanceFile, geneId_2_genome, geneId_2_cog, geneId_2_length, dictRemovedConnections, alignments)

    defaultDist = int(defaultDist * 100.0)
    calcResultsMatrixAndPrint(alignments, outfileName, countCutoff, defaultDist, genome_2_set_of_cogs, dictOG_averageLen)

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



def main(argv=None):
    parser = argparse.ArgumentParser(description='This program calculates mOTU abundances for one sample',
                                     add_help=True)
    parser.add_argument('inFileNameList', action="store",
                        help='File with a list of files with distances or similarities in format ID1\tID2\tdist/sim, which IDs belong together is specified in geneMapFile')
    parser.add_argument('lengthFile', action="store", help='file with lengths of sequences used in blast.')
    parser.add_argument('fileIDsList', action="store", help='file providing a list of genomes/IDs that should be kept')
    parser.add_argument('geneMapFile', action="store", help='file mapping gene name to genome and OG, tab delimited')
    parser.add_argument('removedFile', action="store", help='File that provides previously removed distances.')
    parser.add_argument('--countCutoff', '-c', action="store", dest='countCutoff', type=int, default=0,
                        help='Minimal number of Marker gene distance pairs needed to compute combined distance.')
    parser.add_argument('--defaultDist', '-d', action="store", dest='defaultDist', type=float, default=0.0,
                        help='Default distance for missing distances, should be set to or just below the cutoff for dist calculations.')
    parser.add_argument('outfileName', action="store", help='Output file.')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()


    runFunctions(args.inFileNameList, args.geneMapFile, args.fileIDsList, args.lengthFile, args.removedFile, args.countCutoff, args.defaultDist,
                 args.outfileName)  
    return 0  # success


if __name__ == '__main__':
    status = main()
    sys.exit(status)
