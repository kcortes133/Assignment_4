# Author: Katherina Cortes
# Date: September 29, 2021
# Purpose: Different possible statistical tests for networks and subnetworks

import matplotlib.pyplot as plt


# edge density defined in the paper as edge count
# bool avgDensity if want avg num of edges/node -> good if networks have diff num of nodes
#   default set to false
# @param network: network dictionary with gene nodes and interactions edges
# @returns density: number of edges in network
def calcEdgeDensity(network):
    density = 0

    # each edge is represented twice so
    for node in network:
        density += len(network[node])
    density = density/2

    return density


# edge density using edge weights instead of number of edges
# bool avgDensity if want avg num of edges/node -> good if networks have diff num of nodes
#   default set to false
# @param network: network dictionary with gene nodes and interactions edges
# @returns density: total edge density in network
def calcEdgeDensityW(network):
    density =0

    for node in network:
        for edge in network[node]:
            density += float(network[node][edge])
    density = density/2
    return density


# @param: densities
# @displays: histogram of densities
def histogram(densities):
    # Plot Histogram on x
    plt.hist(densities)
    plt.gca().set(title='Density Distribution Histogram', ylabel='Density');
    plt.show()
    return


# @param dens1: list of density distribution
# @param dens2: list of different density distribution
# @displays: histogram of the two densities
def overlappingHistogram(dens1, dens2):

    # Plot Histogram on x
    plt.hist([dens1, dens2])
    plt.gca().set(title='Edge Density Histogram', ylabel='Number of Nodes');
    plt.legend(loc='upper right')
    plt.show()


# @param lociSubN: list of loci subnetworks dictionaries
# @param coFSubN: list of cofunctional subnetworks dictionaries
# @returns pval: int pval for loci subnetworks with cofunctional subnetworks as null hypothesis
def empiricalPVal(lociSubN, coFPopDensities):
    # use avg density of final population in random trial
    # P value representing fraction of random trials producing final pop of subnetworks
    # with higher avg density than the avg density seen with true loci inputs

    # get average network density for all loci networks
    lociDensity = 0
    for subNet in lociSubN:
        # calculate the edge density
        tempLD = calcEdgeDensityW(subNet)
        lociDensity += tempLD
    lociDensity = lociDensity/len(lociSubN)

    print(lociDensity)

    # calculate p-val by using cof (null analysis) density distribution
    # and avg loci density
    coFDensities = sorted(coFPopDensities)
    pos = 0
    p=0
    densPos = 0
    for c in coFDensities:
        if c <=lociDensity:
            p = pos
            densPos = c
        pos +=1

    pval = (len(coFDensities) - p)/len(coFDensities)

    # plot a histogram of cof density distribution and avd loci density
    # put pval next to avg loci dashed line
    plt.hist(coFDensities)
    plt.axvline(densPos, color='k', linestyle='dashed', linewidth=1)
    plt.title('Empirical P-Value')
    min_ylim, max_ylim = plt.ylim()
    plt.text(lociDensity, max_ylim*0.9, 'Pval = '+ str(pval))
    plt.show()
    return pval
