# Author: Katherina Cortes
# Date: November 16, 2021
# Purpose: Take a tab delimited .gmz file of gene sets, a given tab-delimited STRING file
#   create a sub network of the gene interactions from the input file using the STRING file
#   get statistical significance

# START
# Generate the initial population
# Calculate fitness
# WHILE
#    Mutate
#    Selection
#    Mate
#    Compute fitness
# UNTIL population does not have significant change in density
# STOP
import copy
import random, statistics
import matplotlib.pyplot as plt

# @param network: dictionary of subnetwork of genes
# @param connections: dictionary of all genes and their interactions
# @returns network: dictionary of gene subnetwork with edges
def makeEdges(network, connections):
    for gene in network:
        for gene2 in connections[gene]:
            if gene2 in network:
                weight = connections[gene][gene2]
                network[gene] = {gene2:weight}
    return network


# Mutate at 5% chance each gene in all subnetworks in a list
# mutated genes are chosen from same loci as old gene
# @param subnetworks: list of subnetworks to mutate
# @param lociLists: list of list of loci and their genes
# @param connections: all genes and connections to other genes
# @returns newPop: list of mustated subnetworks
def mutation(subnetworks, lociLists, connections):
    newPop = copy.deepcopy(subnetworks)

    for network in range(len(subnetworks)):
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

        newPop[network] = makeEdges(newPop[network], connections)

    return newPop


# calculate selection scores for networks to use in mating
# @param subnetworks: list of subnetworks
# @returns scores: list of scores for corresponding subnetworks
def calculateSelectionScores(subnetworks):
    scores = []
    for network in subnetworks:
        edgeDensity = statistics.calcEdgeDensityW(network)
        scores.append(edgeDensity)
    return scores


# Mate each network with a randomly chosen one. Networks with higher selection scores
# are more likely to be randomly chosen for mating.
# @param subnetworks: list of subnetworks to mate
# @param lociLists: list of lists of loci and their genes
# @param connections: dictionary of genes and interactions
# @returns subnetworks: list of mated subnetworks
def mating(subnetworks, lociLists, connections):
    origSNetworks = copy.deepcopy(subnetworks)
    selScores = calculateSelectionScores(subnetworks)
    selList = []
    networkIndex = 0
    for score in selScores:
        for s in range(int((score*10)+1)):
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
            if gene1 in network:
                network[gene1] = {}


        makeEdges(network, connections)
    return subnetworks


# calculate edge density of networks using edge weights
# use to determine when to stop the genetic algorithm
# @param population: networks to calculate density of
# @returns eDensity: total edge density for all networks
def calcPopEdgeDensity(population):
    eDensity = 0
    for network in population:
        #eDensity += statistics.calcEdgeDensity(network)
        eDensity += statistics.calcEdgeDensityW(network)
    return eDensity


# Mutate
# Mate
# keep going until overall edge density has not improved by more than 0.5%
# @param subnetworks: list of subnetworks to mate
# @param lociLists: list of lists of loci and their genes
# @param connections: dictionary of genes and interactions
# @returns newPop: list of new subnetworks generated from the genetic algorithm
def geneticAlg(subnetworks, lociLists, connections):
    change = 100
    generation = 0
    changes = []
    generations = []
    while change > 0.005:
        # genetic algorithm
        startingEDensity = calcPopEdgeDensity(subnetworks)
        newPop = mutation(subnetworks, lociLists, connections)
        newPop = mating(newPop, lociLists, connections)

        #caculate change
        endEDensity = calcPopEdgeDensity(newPop)
        change = abs((endEDensity - startingEDensity)/startingEDensity)
        subnetworks = newPop
        generation += 1

        print(generation)
        print(change)
        changes.append(change)
        generations.append(generation)


    with open('GA_Generations_Stats', 'w') as f:
        for c in changes:
            f.write(str(c) + '\n')

    plt.figure(figsize=(9, 3))

    plt.plot(generations, changes)
    plt.show()


    return newPop