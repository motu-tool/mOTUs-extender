#!/usr/bin/env python

# source /g/bork/mende/bin/virtualenv/py2.7_bio1/bin/activate

import os
import sys
import argparse
# sys.path.append("/g/bork/mende/lib/python/lib/python2.7/site-packages")
#import hcluster

import numpy
import scipy.cluster.hierarchy
import scipy.sparse.csgraph
import scipy.sparse
from collections import defaultdict

import scipy.cluster.hierarchy as module_cluster


#import fastcluster
#scipy.sparse.csgraph.connected_components(csgraph, directed=True, connection='weak', return_labels=True)

#import fastcluster

# read file with all ids
#listInfileData = readSingleColumnFileToList(infileName)
def readSingleColumnFileToList(infileName):

    infile = open(infileName, "r")
    listInfileData = []
    for strLine in infile:
        strLine = strLine.rstrip('\n')
        if (" " in strLine):
            arrLine = strLine.split(' ')
            strLine = arrLine[0]

        listInfileData.append(strLine)

    return(listInfileData)


def parse2columnFileDictLists(infileName):
    dictIn = defaultdict(list)
    dictInvOut = defaultdict(list)

    infile = open(infileName, "r")

    for strInLine in infile:

        strInLine = strInLine.rstrip('\n')
        arrInLine = strInLine.split('\t')

        strColumn_1_entry = arrInLine[0]
        strColumn_2_entry = arrInLine[1]

        dictIn[strColumn_1_entry].append(strColumn_2_entry)
        dictInvOut[strColumn_2_entry].append(strColumn_1_entry)

    return(dictIn, dictInvOut)


def reformatDict(dictID2content, dictOldIDs2NewIDsList, listIDs):
    dictOut = {}

    for ID in dictID2content:
        listNewID = dictOldIDs2NewIDsList[ID]
        for newID in listNewID:
            if (newID in listIDs):
                dictOut[newID] = dictID2content[ID]
                # print str(newID) + "\t" + str(dictID2content[ID])
    return(dictOut)


def reformatDict2(dictID2list, dictOldIDs2NewIDsList, listIDs):
    dictOut = defaultdict(list)

    for clusterID in dictID2list:
        listClusterIDs = dictID2list[clusterID]
        for ID in listClusterIDs:
            listNewID = dictOldIDs2NewIDsList[ID]
            for newID in listNewID:
                if (newID in listIDs):
                    dictOut[clusterID].append(newID)
                    # print str(clusterID) + "\t" + str(newID)

    return(dictOut)


def parseClustering(clusteringFileName, prefix=""):
    dictCluster2Members = defaultdict(list)
    dictID2Cluster = {}
    listClusterNames = []

    clusteringFile = open(clusteringFileName, "r")

    for strInLine in clusteringFile:
        strInLine = strInLine.rstrip('\n')
        arrInLine = strInLine.split('\t')
        strClusterName = arrInLine[0]

        # cluster name can be modified using a prefix,
        strClusterName = prefix + strClusterName

        listClusterNames.append(strClusterName)

        strClusterMembers = arrInLine[1]
        arrClusterMembers = strClusterMembers.split(';')

        for strClusterMember in arrClusterMembers:
            if (not strClusterMember in dictID2Cluster):
                dictID2Cluster[strClusterMember] = strClusterName
                dictCluster2Members[strClusterName].append(strClusterMember)
            else:
                print(strClusterMember + " " + strClusterName +
                      " assigned to multiple taxonomy ranks in " + clusteringFileName + ".")

    return(dictCluster2Members, dictID2Cluster, listClusterNames)


#(distanceMatrix, distanceMatrixEdges, dictKey2Position) = initMatrices(listKeys)
def initMatrices(listKeys):
    numKeys = len(listKeys)
    print(numKeys)
    distanceMatrix = scipy.sparse.dok_matrix(
        (numKeys, numKeys), dtype=numpy.float32)
    distanceMatrixEdges = scipy.sparse.dok_matrix(
        (numKeys, numKeys), dtype=int)

#	dictKey2Position = defaultdict(int)
#	keyCount = 0
#	for i in listKeys:
#		dictKey2Position[i] = keyCount
#		keyCount += 1

    dictKey2Position = {}
    for index, item in enumerate(listKeys):
        dictKey2Position[item] = index

    return(distanceMatrix, distanceMatrixEdges, dictKey2Position)

#(distanceMatrixEdges, distanceMatrix) = fillDistMatrix(distanceMatrixEdges, distanceMatrix, dictKey2Position)


def fillDistMatrix(distFileName, distanceMatrixEdges, distanceMatrix, dictKey2Position, dictID2Cluster, dictCluster2Members, type, listKeys, floatCutoff, forceClusterSeparation, hasPreClustering):

    setMatrixEntries = set()
    distFile = open(distFileName, "r")

    for strLine in distFile:
        strLine = strLine.strip()
        arrLine = strLine.split('\t')

        strID_1 = arrLine[0]
        strID_2 = arrLine[1]

        if (" " in strID_1):
            arrID_1 = strID_1.split(' ')
            strID_1 = arrID_1[0]

        if (" " in strID_2):
            arrID_2 = strID_2.split(' ')
            strID_2 = arrID_2[0]

        strDist = arrLine[2]
        floatDist = float(strDist)

        if (strID_1 in dictKey2Position and strID_2 in dictKey2Position):
            position_1 = dictKey2Position[strID_1]
            position_2 = dictKey2Position[strID_2]

            if (position_2 <= position_1):
                position_temp = position_1
                position_1 = position_2
                position_2 = position_temp

            tuplePositions = (position_1, position_2)

            if not (tuplePositions in setMatrixEntries):

                # reformat ID/dist to dist
                if (type == "ID"):
                    total = 100.0
                    newtotal = 1.0
                    floatDist = (float(total) - floatDist) * \
                        (float(newtotal)/float(total))

                if (type == "ID_1"):
                    total = 1.0
                    newtotal = 1.0
                    floatDist = (float(total) - floatDist) * \
                        (float(newtotal)/float(total))

                if(strID_1 in dictID2Cluster and strID_2 in dictID2Cluster and dictID2Cluster[strID_1] == dictID2Cluster[strID_2]):
                    floatDist = 0.0
                elif(floatDist <= 0.0001 and hasPreClustering):
                    floatDist = 0.0001

                # original distance is used for connected components computation
                if (floatDist <= (floatCutoff+0.001)):
                    distanceMatrixEdges[position_1, position_2] = 1
                    distanceMatrixEdges[position_2, position_1] = 1
                if (forceClusterSeparation):
                    if(strID_1 in dictID2Cluster and strID_2 in dictID2Cluster and dictID2Cluster[strID_1] != dictID2Cluster[strID_2]):
                        floatDist = 10000.0

                # distances for clustering are preprocessed to keep specI clusters as in the original clustering
                distanceMatrix[position_1, position_2] = 2.0 + floatDist
                distanceMatrix[position_2, position_1] = 2.0 + floatDist

                setMatrixEntries.add(tuplePositions)

    numKeys = len(listKeys)
    for i in range(numKeys):
        distanceMatrix[i, i] = 2.0
        distanceMatrixEdges[i, i] = 1

    for strCluster in dictCluster2Members:

        listIDs = dictCluster2Members[strCluster]
#		listPositions = [ dictKey2Position[x] for x in listIDs]
        # print(len(listPositions))
        for strID_1 in listIDs:
            for strID_2 in listIDs:
                if (strID_1 in dictKey2Position and strID_2 in dictKey2Position):
                    position_1 = dictKey2Position[strID_1]
                    position_2 = dictKey2Position[strID_2]
                    distanceMatrix[position_1, position_2] = 2.0
                    distanceMatrixEdges[position_1, position_2] = 1

#		for i in listPositions:
#			for j in listPositions:
#				distanceMatrix[i, j] = 1.0
#				distanceMatrixEdges[i, j] = 1

    return(distanceMatrixEdges, distanceMatrix)


#(dictComponent2IDsNumber, listSingletons) = getConnectedComponents(distanceMatrixEdges, dictKey2Position)
def getConnectedComponents(distanceMatrixEdges, listKeys):

    listSingletons = []
    dictComponent2IDsNumber = defaultdict(list)
    distanceMatrixEdges = scipy.sparse.csr_matrix(distanceMatrixEdges)

    #N_components, component_list = scipy.sparse.csgraph.connected_components(distanceMatrixEdges, directed=False)
    N_components, component_list = scipy.sparse.csgraph.connected_components(
        distanceMatrixEdges)

    # reformat component_list to access list of clusters
# >>> component_list
#array([-2,  0,  0,  0,  1,  1,  1, -2,  2,  2], dtype=int32)
    if (len(listKeys) != len(component_list)):
        print("Weird 1")

    print(len(component_list))
    # change to enumerate
    for index, item in enumerate(component_list):

        # for i in range(len(component_list)):
        #clusterNumber = component_list[index]
        clusterNumber = item
        # this could be different for scipy.sparse.csgraph.connected_components
        if (clusterNumber != -2):
            dictComponent2IDsNumber[clusterNumber].append(index)
        else:
            listSingletons.append(index)

    print("Preclustering results: Connected Components: " +
          str(N_components) + "\tSingletons: " + str(len(listSingletons)))

    return(dictComponent2IDsNumber, listSingletons)

# takes list of nodeNumbers from a distance matrix
#dictCluster2IDs = defaultdict(list)


def getSubClustering(listIDs_subcluster, setClusters, distanceMatrix, strLinkMethod, strCutoffMethod, floatCutoff, dictCluster2IDs, dictID2Cluster, listKeys, forceClusterSeparation):

    listIDs_subcluster = sorted(listIDs_subcluster)

    numberOfIDs = len(listIDs_subcluster)
    matrixInternalDist = numpy.ones(
        (numberOfIDs, numberOfIDs), dtype=numpy.float32)
    # print(numberOfIDs)
    # print(listIDs_subcluster)
    internalIndex1 = 0
    for index1 in listIDs_subcluster:
        key1 = listKeys[index1]
        internalIndex2 = 0
        for index2 in listIDs_subcluster:
            key2 = listKeys[index2]
            dist = distanceMatrix[index1, index2]
            if(dist >= 2.0):
                matrixInternalDist[internalIndex1][internalIndex2] = dist - 2.0

            if(key1 in dictID2Cluster and key2 in dictID2Cluster and dictID2Cluster[key1] == dictID2Cluster[key2]):
                matrixInternalDist[internalIndex1][internalIndex2] = 0.0
#			elif(dist <= 0.0001):
#				matrixInternalDist[internalIndex1][internalIndex2] = 0.0001

            if (forceClusterSeparation):
                if(key1 in dictID2Cluster and key2 in dictID2Cluster and dictID2Cluster[key1] != dictID2Cluster[key2]):
                    matrixInternalDist[internalIndex1][internalIndex2] = 10000.0

            internalIndex2 += 1
        internalIndex1 += 1


#	dictInternalKey2Position = defaultdict(int)
#	internalKeyCount = 0
#	for i in listIDs_subcluster:
#		dictInternalKey2Position[i] = internalKeyCount
#		internalKeyCount += 1

    dictInternalKey2Position = {}
    for index, item in enumerate(listIDs_subcluster):
        dictInternalKey2Position[item] = index

    vectInternalDistMatrix = scipy.spatial.distance.squareform(
        matrixInternalDist)
    linkDistMatrix = module_cluster.linkage(
        vectInternalDistMatrix, method=strLinkMethod)
    del vectInternalDistMatrix
    vectorClusters = scipy.cluster.hierarchy.fcluster(
        linkDistMatrix, floatCutoff, strCutoffMethod)

    maxClusterNumber = 0
    # this is a bit weird, but this ensures that clusters are continous
    if (len(setClusters) > 0):
        maxClusterNumber = max(setClusters)

    vectorClustersOffset = [x + maxClusterNumber for x in vectorClusters]

    for cluster in vectorClustersOffset:
        if (cluster in setClusters):
            print(cluster)
            print(maxClusterNumber)
            print(listIDs_subcluster)
            print(vectorClusters)

            vectorClustersOffset = [
                x + maxClusterNumber + 1 for x in vectorClusters]

    newClusterNeeded = True

    for strCurrentKeyName, intCluster in zip(listIDs_subcluster, vectorClustersOffset):
        #strCurrentKeyName = index2entry[i]
        if (newClusterNeeded and intCluster in dictCluster2IDs):
            print("Problem detected 1")
            print(strCurrentKeyName)
            print(intCluster)

        if (intCluster <= maxClusterNumber):
            print("Problem detected 2")
            print(strCurrentKeyName)
            print(intCluster)
            print(maxClusterNumber)

        dictCluster2IDs[intCluster].append(strCurrentKeyName)
        newClusterNeeded = False
        setClusters.add(intCluster)

    return(dictCluster2IDs, setClusters)


def runClustering(listIDsFileName, distFileName, floatCutoff, type, strLinkMethod, strCutoffMethod, noConnectedComponents, forceClusterSeparation, preClusteringFileName, preClusterMapFile, outfileName):

    outfile = open(outfileName, "w")
    setClusters = set()

    hasPreClustering = True

    if (type == "ID"):
        total = 100.0
        newtotal = 1.0
        floatCutoff = (float(total) - floatCutoff) * \
            (float(newtotal)/float(total))

    if (type == "ID_1"):
        total = 1.0
        newtotal = 1.0
        floatCutoff = (float(total) - floatCutoff) * \
            (float(newtotal)/float(total))

    print(str(floatCutoff))

    listIDs = readSingleColumnFileToList(listIDsFileName)

    if (preClusteringFileName != ""):
        (dictCluster2Members, dictID2Cluster, listClusterNames) = parseClustering(
            preClusteringFileName, prefix="")
        hasPreClustering = True
        if (preClusterMapFile != ""):
            (dictOldIDs2NewIDs, dictNewIDs2OldIDs) = parse2columnFileDictLists(
                preClusterMapFile)
            dictID2Cluster = reformatDict(
                dictID2Cluster, dictOldIDs2NewIDs, listIDs)
            dictCluster2Members = reformatDict2(
                dictCluster2Members, dictOldIDs2NewIDs, listIDs)
    else:
        dictID2Cluster = {}
        dictCluster2Members = {}
        hasPreClustering = False

    listIDs = readSingleColumnFileToList(listIDsFileName)
    (distanceMatrix, distanceMatrixEdges, dictKey2Position) = initMatrices(listIDs)

    (distanceMatrixEdges, distanceMatrix) = fillDistMatrix(distFileName, distanceMatrixEdges, distanceMatrix, dictKey2Position,
                                                           dictID2Cluster, dictCluster2Members, type, listIDs, floatCutoff, forceClusterSeparation, hasPreClustering)
    if (not noConnectedComponents):
        (dictComponent2IDsNumber, listSingletons) = getConnectedComponents(
            distanceMatrixEdges, dictKey2Position)
    else:
        listSingletons = []
        dictComponent2IDsNumber = {}
        dictComponent2IDsNumber["All"] = [
            index for index, item in enumerate(listIDs)]

    dictCluster2IDs = defaultdict(list)
    for component in dictComponent2IDsNumber:
        listIDsInComponent = dictComponent2IDsNumber[component]

        if (len(listIDsInComponent) > 1):
            #dictCluster2IDs, setClusters = getSubClustering(listIDsInComponent, setClusters, distanceMatrix, strLinkMethod, strCutoffMethod, floatCutoff, dictCluster2IDs)

            dictCluster2IDs, setClusters = getSubClustering(listIDsInComponent, setClusters, distanceMatrix, strLinkMethod,
                                                            strCutoffMethod, floatCutoff, dictCluster2IDs, dictID2Cluster, listIDs, forceClusterSeparation)

        elif(len(listIDsInComponent) == 1):
            listSingletons.append(listIDsInComponent[0])

    print("Singletons: " + str(len(listSingletons)))

    for cluster in dictCluster2IDs:

        listClusterIDs = dictCluster2IDs[cluster]
        listClusterIDsOriginalIDs = [str(listIDs[x]) for x in listClusterIDs]

        strWriteLine = str(cluster) + "\t" + \
            ";".join(listClusterIDsOriginalIDs) + "\n"
        outfile.write(strWriteLine)
    if (len(dictCluster2IDs) > 0):
        maxCluster = max(dictCluster2IDs.keys())
    else:
        maxCluster = 0

    for x in listSingletons:
        maxCluster += 1
        strID = str(listIDs[x])
        strWriteLine = str(maxCluster) + "\t" + str(strID) + "\n"
        outfile.write(strWriteLine)


########################
def main(argv=None):
    if(not argv):
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(
        description='This program is used to cluster sequences into clusters based on sequence identity values', add_help=True)
    parser.add_argument('distancefile', action="store",
                        help='file with distances or similarities in format ID1\tID2\tdist/sim')
    parser.add_argument('listIDsFileName', action="store",
                        help='file listing all IDs to be used for the clustering')
    parser.add_argument('cutoff', action="store", type=float,
                        help='cutoff for the clustering')
    parser.add_argument('outfile', action="store",
                        help='prefix for outfile(s)')
    parser.add_argument('--distOrID', '-d', action='store', dest='type', default="ID",
                        help="'dist' or 'ID'; if 'ID' then its transformed into a distance inside, defaults to 'ID' [ID_1]")
    parser.add_argument('--linkMethod', '-l', action='store', dest='linkMethod', default='average',
                        help="linkMethod: 'single','complete','average','weighted', defaults to average")
    parser.add_argument('--noConnectedComponents', '-n', action='store_true', default=False,
                        dest='noConnectedComponents', help='Set to not use the connected component preclustering')
    parser.add_argument('--forceClusterSeparation', '-s', action='store_true', default=False,
                        dest='forceClusterSeparation', help='Force preclusters to stay separate from each other')
    parser.add_argument('--preClusteringFile', '-p', action='store', dest='preClusteringFile',
                        default="", help='File with a preclustering e.g. specI clusters. The preclusters will be kept.')
    parser.add_argument('--preClusterMap', '-m', action='store', dest='preClusterMapFile', default="",
                        help='Tab delimited 2 column file mapping IDs in the preclustering to IDs in the distance file. E.g. specI clusters to geneID in the currently clustered OG')
    parser.add_argument('--fastcluster', '-f', action='store_true', dest='fastcluster',
                        default=False, help='Set to use fastcluster module for clustering')
    parser.add_argument('--version', '-v', action='version',
                        version='%(prog)s 1.0')
    args = parser.parse_args()

    if (args.fastcluster):
        try:
            import fastcluster as module_cluster
        except ImportError:
            import scipy.cluster.hierarchy as module_cluster

    strCutoffMethod = "distance"
    runClustering(args.listIDsFileName, args.distancefile, args.cutoff, args.type, args.linkMethod, strCutoffMethod,
                  args.noConnectedComponents, args.forceClusterSeparation, args.preClusteringFile, args.preClusterMapFile, args.outfile)
    return 0		# success


if __name__ == '__main__':
    status = main()
    sys.exit(status)


# fcluster: use 2 modes:
#	inconsistency
#	distance


# linkage use 4 modes:
# method='single'
# method='complete'
# method='average'
# method='weighted'
