# Author: Katherina Cortes
# Date: September 29, 2021
# Purpose: Create subnetworks and full network from FA loci genes and STRING database

import random, math

# @param numNetworks number of subnetworks to make from the locilists
# @param fullNetwork gene:gene interactions dictionary
# @param lociLists list of lists where each sublist is the genes at one loci
# @returns subnetworks list of the loci subnetwork dictionaries
def makeLociSubnetworks(numNetworks, fullNetwork, lociLists):
    subNetworks = []
    for i in range(numNetworks):
        subNetworks.append(createSubnetwork(fullNetwork, lociLists))
    return subNetworks


# @param fullNetwork: gene:gene intearctions dictionary
# @param lociLists: list of lists where each sublist is the genes at one loci
# @returns subNetwork: dictionary subnetwork of loci network with one random gene from each loci
def createSubnetwork(fullNetwork, lociLists):
    subNetwork = {}
    # for loci in network:
    # pick gene randomly
    for loci in lociLists:
        randGene = random.choice(list(loci))
        subNetwork[randGene] = {}

    # get edges
    for gene1 in subNetwork:
        # connect nodes with edges
        for gene2 in fullNetwork[gene1]:
            if gene2 in subNetwork:
                weight = fullNetwork[gene1][gene2]
                subNetwork[gene1][gene2] = weight

    return subNetwork



# limitations: once you get higher the bins get more empty
# @param fullNetwork: gene:gene intearctions dictionary
# @param numBins: number of bins to separate the bins into by edge density
# @return numBins: nodes from fullNetwork organized into bins
def makeFixedBins(fullNetwork, numBins):
    fullNetworkBins = []
    # dict should be node: num Edges
    # sort dict by values
    # 128 bins like in paper -> equally spaced

    # sort nodes by number of edges
    networkSorted = sorted(fullNetwork, key=lambda k: len(fullNetwork[k]), reverse=False)
    maxEdges = len(fullNetwork[networkSorted[-1]])
    binSize = round(maxEdges/numBins)

    for i in range(numBins):
        fullNetworkBins.append([])

    # put nodes in bins according to edge number
    # bins can be variable sizes
    for node in networkSorted:
        nodeDensity = len(fullNetwork[node])
        nodeBin = math.floor(nodeDensity/binSize)
        fullNetworkBins[nodeBin].append(node)

    return fullNetworkBins


# @param fullNetwork: gene:gene intearctions dictionary
# @param numBins: number of bins to separate the bins into by edge density
# @return numBins: nodes from fullNetwork organized into bins
def makeQuantileBins(fullNetwork, numBins):
    fullNetworkBins = []
    # dict should be node: num Edges
    # sort dict by values
    # 128 bins like in paper -> equally spaced

    # sort nodes by number of edges
    networkSorted = sorted(fullNetwork, key=lambda k: len(fullNetwork[k]), reverse=False)
    numNodes = len(fullNetwork)
    binSize = round(numNodes/numBins)

    # make each bin the same size
    # nodes in bins according to numbder of edges
    for i in range(numBins):
        binStart = i*binSize
        binEnd = (i+1)*binSize
        fullNetworkBins.append(networkSorted[binStart:binEnd])
    return fullNetworkBins


# assumptions: no node in the lociSubN will have over 10k edges
# create a subnetwork of the cofunctional network
# where nodes have similar density to nodes in the loci sub network
# @param fullNetwork: gene:gene interactions dictionary
# @param fullNetworkBins: nodes from fullNetwork organized into bins by edge density
# @param lociSubN: subnetworks made from the loci genes
# @return coFSubnetworks: subnetworks of fullNetwork where nodes random selected to be equivalent to nodes in
#   loci subnetworks
def makeCoFSubnetworks(fullNetwork, fullNetworkBins, lociSubN,):
    # for each node in lociSubN
    #   pick node from networkBins where numEdges is similar for both
    # make edges between picked nodes
    # if bin is empty go to the closest... look above and below

    coFSubnetworks = []
    lociD = []
    coFD = []

    for subNetwork in lociSubN:
        tempSubnetwork = {}
        for node in subNetwork:
            # check node density from the fullNetwork so that its equivalent
            if node in fullNetwork:
                nodeDensity = len(fullNetwork[node])
            else:
                nodeDensity = 0
            lociD.append(nodeDensity)

            # loop through bin
            # first bin that smaller than last element pick from that bin
            for bin in fullNetworkBins:
                if len(bin) > 0:
                    binMax = len(fullNetwork[bin[-1]])
                else:
                    binMax = -1
                if binMax >= nodeDensity:
                    # pick random from bin
                    # make sure random not already in subnetwork
                    node = random.choice(bin)
                    if node in tempSubnetwork:
                        node = random.choice(bin)
                    tempSubnetwork[node] = {}
                    coFD.append(len(fullNetwork[node]))
                    break

            # pick node with similar density
            # density is from
            # make sure node isn't already in network
            # if bin is empty go to nearest -> look above and below
        # make edges between picked nodes
        for tempNode in tempSubnetwork:
            for edge in fullNetwork[tempNode]:
                if edge in tempSubnetwork:
                    tempSubnetwork[tempNode][edge] = fullNetwork[tempNode][edge]

        coFSubnetworks.append(tempSubnetwork)
    return coFSubnetworks
