# Author: Katherina Cortes
# Date: September 29, 2021
# Purpose: get data from input and STRING files

import operator
from functools import reduce


# assumptions:
#   input is tab delimited
#   first two indices in line should be thrown away
#   gene in 'locus for gene' will be repeated in the list
#   locus number is unimportant for network
#   gene interactions should be examined across and within all gene sets
#
# limitations:
#   gene interactions only shown if in STRING network
#   only known gene interactions and given genes
#
#   Returns a list of genes.
#   Expects a .txt tab-delimited file as input.
#   First two columns in file rows are discarded
#
#   @param fileIn file with genes to be read in
#   @param asLoci bool of whether to have list of loci lists or just one list
#   @return genes list of genes
def readInput(fileIn, asLoci=True):
    # read in GMT file
    # separate file into list by tabs
    fI = open(fileIn, 'r')
    genes = fI.read().strip().split('\n')
    genes = [g.split('\t') for g in genes]

    # Delete first two columns
    for g in genes:
        if len(g) > 1:
            del g[0]
            del g[0]

    # condense list of lists into one-dimensional list
    if not asLoci:
        genes = reduce(operator.add, genes)

    return genes


# assumptions:
#   input is organized as protein'\t'protein'\t'weight'\n'
#   protein1:protein2:weight is also represented in the list as protein2:protein1:weight
#
#   limitations:
#   dictionary would not work if the protein:protein interaction was only represented once
#
#   Returns dictionary of protein interactions ordered protein:protein:weight
#   Expects tab-delimited STRING file as input
#
#   @param   stringFile STRING file of protein-protein interactions
#            organized as protein'\t'protein'\t'weight'\n'
#   @returns interactions dict of protein interactions
def makeInteractionNetwork(stringFile):
    # dictionary of dictionary 1st protein : 2nd protein: weight
    interactions = {}
    with open(stringFile, 'r') as f:
        for w, i in enumerate(f):
            gene1, gene2, weight = i.strip('\n').split('\t')
            if gene1 in interactions:
                interactions[gene1][gene2] = weight
            else:
                interactions[gene1] = {gene2: weight}
    return interactions


#   @param   genes list of genes to check interactions for (gotten from readInput(file))
#   @param   interactionsNetwork dict of protein-protein interactions from makeInteractionNetwork(stringFile)
#   @returns geneInteractions dict of gene interactions from input GMT file
def makeNetwork(genes, interactionsNetwork):
    geneInteractions = {}
    genes = reduce(operator.add, genes)
    for gene1 in genes:
        # input gene1 check if connected to any other gene
        geneInteractions[gene1] = {}
        if gene1 in interactionsNetwork:
            for gene2 in genes:
                if gene2 in interactionsNetwork[gene1]:
                    if gene2 in genes:
                        # make sure not to duplicate edges
                        #if gene1 not in geneInteractions[gene2]:
                        geneInteractions[gene1][gene2] = interactionsNetwork[gene1][gene2]
    return geneInteractions