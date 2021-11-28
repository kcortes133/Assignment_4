# Author: Katherina Cortes
# Date: November 16, 2021
# Purpose: Take a tab delimited .gmz file of gene sets, a given tab-delimited STRING file
#   create a sub network of the gene interactions from the input file using the STRING file
#   get statistical significance

# //TODO
# START
# Generate the initial population
# Compute fitness
# REPEAT
#    Selection
#    Crossover
#    Mutation
#    Compute fitness
# UNTIL population has converged
# STOP
import random, statistics


def makeEdges(network, connections):
    b = len(network)
    for gene in network:
        for gene2 in connections[gene]:
            if gene2 in network:
                weight = connections[gene][gene2]
                network[gene] = {gene2:weight}
    e = len(network)
    if b != e:
        print('i hate me')
    return network

def mutation(subnetworks, lociLists, connections):
    newPop = []
    [newPop.append(dict.fromkeys(subnetworks[x], {})) for x in range(len(subnetworks))]

    for network in range(len(subnetworks)):
        b = len(subnetworks[network])
        temp = subnetworks[network]
        for gene in subnetworks[network]:
            # mutation at 5% chance
            subnetworks[network][gene] = {}
            mChance = random.randint(0, 100)
            if mChance < 5:
                for loci in lociLists:
                    if gene in loci:
                        mutatedG = random.choice(loci)
                        while mutatedG == gene:
                            mutatedG = random.choice(loci)

                        # replace gene with mutatedG
                        newPop[network][mutatedG] = {}
                        newPop[network].pop(gene)
                        e = len(newPop[network])
                        if b != e:
                            print(gene)
                            print(mutatedG)

        newPop[network] = makeEdges(newPop[network], connections)
        e = len(newPop[network])
        if len(newPop[network]) < 12:
            print(newPop[network])
        if b != e:
            print('fuck this')
            print(temp)
            print(newPop[network])

    return newPop


def calculateSelectionScores(subnetworks):
    scores = []
    for network in subnetworks:
        edgeDensity = statistics.calcEdgeDensity(network)
        scores.append(edgeDensity**3)
    return scores


def mating(subnetworks, lociLists, connections):
    origSNetworks = subnetworks
    selScores = calculateSelectionScores(subnetworks)
    selList = []
    networkIndex = 0
    for score in selScores:
        for s in range(int(score+1)):
            selList.append(networkIndex)
        networkIndex +=1

    for network in subnetworks:
        mateIndex = random.choice(selList)
        # get corresponding network to mate with network
        mateN = list(origSNetworks[mateIndex].keys())
        # cant change dict keys in a loop
        nodes = list(network.keys())
        # mate
        # genes must come from same loci
        for gene1 in range(len(nodes)):
            for gene2 in mateN:
                for loci in lociLists:
                    if gene1 in loci and gene2 in loci:
                        choice = random.randint(0,1)
                        if choice:
                            network.pop(gene1)
                            network[gene2] = {}
        makeEdges(network, connections)
    return subnetworks


def calcPopEdgeDensity(population):
    eDensity = 0
    for network in population:
        eDensity += statistics.calcEdgeDensity(network)
    return eDensity


def geneticAlg(subnetworks, lociLists, connections):
    change = 100
    generations = 0
    while change > 0.005:
        startingEDensity = calcPopEdgeDensity(subnetworks)
        newPop = mutation(subnetworks, lociLists, connections)
        newPop = mating(newPop, lociLists, connections)
        #caculate change
        endEDensity = calcPopEdgeDensity(newPop)

        change = abs((endEDensity - startingEDensity)/startingEDensity)
        subnetworks = newPop
        generations+=1
        print(generations)