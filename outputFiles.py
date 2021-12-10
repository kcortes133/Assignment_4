# Author: Katherina Cortes
# Date: December 4, 2021
# Purpose: output files describing gene rank, subnetwork score, p-value and statistics
#   about the genetic algorithm generations


import geneticAlgorithm


# assumptions: assumes top 10 networks dont have the same score
# @param pval: pval for genetic algorithm output
# @param networks:
# @param topNets:
# @outputs: file with
def outputNetworks(pval, networks, topNets):
    networkScores = geneticAlgorithm.calculateSelectionScores(networks)
    scoredNetworks = dict(zip(networkScores, networks))
    sortedScores = sorted(networkScores, reverse=True)

    fileN = 'Day3_Output_Network'

    for i in range(topNets):
        network = scoredNetworks[sortedScores[i]]
        fName = fileN + str(i) + '_pval'+ str(pval) + '.txt'
        with open(fName, 'w') as f:
            for gene1 in network:
                for gene2 in network:
                    if gene1 != gene2:
                        if gene2 in network[gene1]:
                            weight = network[gene1][gene2]
                            f.write(gene1 + '\t' + gene2 + '\t' + str(weight) +'\n')
                        else:
                            f.write(gene1 + '\t' + gene2 + '\t' + 'NA' +'\n')

    return


# @param lociLists:
# @param geneScores:
def outputGeneScoresinLoci(geneScores, lociLists):
    fileN = 'Day3_Output.gmt'
    with open(fileN, 'w') as f:
        for loci in lociLists:
            for gene in loci:
                f.write(gene + ':' + str(geneScores[gene]) + '\t')
            f.write('\n')

    return