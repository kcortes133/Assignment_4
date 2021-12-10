# Author: Katherina Cortes
# Date: October 28,2021
# Purpose: Calculate gene score of each gene in input file. Gene score is defined as the
#   amount of edges it contributes to a subnetwork that consists of one gene chosen from each loci
#   gene scores are then averaged to get one gene score. This score indicates relative connectivity to
#   genes in other loci


# Calculate gene scores for each gene for each subnetwork
# Assumption: each gene only appears once in one locus
# genes cant be connected to themselves
# @param lociSubNs: List of dictionaries of subnetworks
# @param lociLists: list of list of genes separated by loci
# @param interactions: dictionary of genes and their edges
# @returns geneScores: dictionary of genes and the list of gene scores calculated from each subnetwork
def getGeneScores(lociSubNs, lociL, interactions):
    geneScores = {}

    # get subNetwork
    for subN in lociSubNs:
        # pick a gene from a loci
        for lociGene in subN:
            # figure out which loci it comes from
            for l in lociL:
                if lociGene in l:
                    loci = l
                    break

            # get all the genes in the specified loci
            for gene in loci:
                # get number of connections in the subNetwork with chosen gene
                # make sure not looking for a connection with loci gene
                geneS = 0
                # check if there are edges with other genes in subnetwork
                for diffGene in subN:
                    # make sure not checking for edge with the gene from same loci
                    # gene score based on edge weights instead of number
                    if diffGene != lociGene and diffGene in interactions[gene]:
                        #geneS += 1
                        geneS = float(interactions[gene][diffGene])

                if gene in geneScores:
                    geneScores[gene].append(geneS)
                else:
                    geneScores[gene] = [geneS]

    return geneScores


# Average the gene Score lists to get one score for each gene
#
# @param geneScores: a dictionary of lists of gene scores for each gene
# @returns geneSAvg: dictionary of average gene score for each gene
def getGeneScoreAvg(geneScores):
    geneSAvg = {}
    for gene in geneScores:
        geneSAvg[gene] = sum(geneScores[gene])/len(geneScores[gene])
    return geneSAvg


# Get the genes with the highest gene scores from each loci. Number of genes defined
# by numGenes
#
# @param geneAvg: dictionary of genes and gene scores
# @param lociLists: list of list of genes separated by loci
# @param numGenes: number of genes to get from each loci
# @returns genes: list of numGenes from each loci
def getTopLociGenes(geneAvg, lociLists, numGenes):
    genes = []
    # account for numGenes being more than length of loci list
    for loci in lociLists:
        lociScores = {}
        for gene in loci:
            lociScores[gene] = geneAvg[gene]
        lociSorted = sorted(lociScores, key=lambda k: lociScores[k], reverse=True)

        n = numGenes
        if numGenes > len(lociSorted):
            n = len(lociSorted)
        genes.extend(lociSorted[:n])

    return genes