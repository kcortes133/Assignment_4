# Author: Katherina Cortes
# Date: October 30, 2021
# Purpose: Take a tab delimited .gmz file of gene sets, a given tab-delimited STRING file
#   create a sub network of the gene interactions from the input file using the STRING file
#   get statistical significance

import argparse, logging, random, time
import networkCreation, fileParsing, geneScoring, geneticAlgorithm, statistics

# arguments:
#   - input file
#   - STRING file
#   - number of bins for cof network default=128
#   - number of networks to run
#   - fixed or quantile


parser = argparse.ArgumentParser(description='Get the statistical significance of the connection of genes vs random '
                                             'genes.')
parser.add_argument('genesFile', metavar='genes', type=str, default='Input.gmt.txt',
                    help='the input file of genes')
parser.add_argument('--interactionsFile', metavar='gene_interactions', type=str, default='STRING.txt',
                    help='the input file of gene interactions')
parser.add_argument('--geneOutFile', type=str, default='geneScores.txt', help='outfile that contains top genes in displayed '
                                                                    'network with calculated gene scores')
parser.add_argument('--networkOutFile', type=str, default='geneNetwork.txt', help='outfile that contains top genes in displayed '
                                                                              'network with calculated gene scores')
parser.add_argument('--calcPVal', type=bool, default=False, help='the number of bins to separate edge densities into')
parser.add_argument('--numBins', type=int, default=128, help='the number of bins to separate edge densities into')
parser.add_argument('--numSubnetworks', type=int, default=5000, help='the number of subnetworks to make')
parser.add_argument('--topGenes', type=bool, default=False, help='graph only top n genes regarless of loci')
parser.add_argument('--numGenes', type=int, default=3, help='number of genes from each loci or total genes'
                                                                 ' if topGenes is true')

args = parser.parse_args()

def main():
    start = time.time()
    random.seed(5)
    visualize = False

    # read in networks
    lociLists = fileParsing.readInput(args.genesFile)
    interactions = fileParsing.makeInteractionNetwork(args.interactionsFile)
    network = fileParsing.makeNetwork(lociLists, interactions)

    # make loci subnetworks
    lociSubN = networkCreation.makeLociSubnetworks(args.numSubnetworks, network, lociLists)

    # calculate gene scores and sort genes by score
    geneScores = geneScoring.getGeneScores(lociSubN, lociLists, network)
    geneAvg = geneScoring.getGeneScoreAvg(geneScores)
    networkSorted = sorted(geneAvg, key=lambda k: geneAvg[k], reverse=True)

    newPop = geneticAlgorithm.geneticAlg(lociSubN, lociLists, network)

    if args.calcPVal:
        numBins = args.numBins
        # make bins for coFunctional subnetwork creation
        qNetworkBins = networkCreation.makeQuantileBins(interactions, numBins)
        fNetworkBins = networkCreation.makeFixedBins(interactions, numBins)
        # make coFunctional random subnetworks
        # need 1000 populations where populations are the 5000 networks
        # do 1000 or 5000 times
        coFPopDensities = []
        for i in range(5000):
            coFSubnetworks = networkCreation.makeCoFSubnetworks(interactions, qNetworkBins, lociSubN)
            popD = 0
            for subCoF in coFSubnetworks:
                popD += statistics.calcEdgeDensity(subCoF)
            coFPopDensities.append(popD/5000)

        # calculate the avg of each population -> makeCoFSubnetorks is one population?
        # calculate the pvalue
        # probability edges using cof distribution is greater than avg of loci edged divided by # of random networks
        pval = statistics.empiricalPVal(lociSubN, coFSubnetworks)

        # make a graph showing the edge density distributions
        coFDensities = []
        for network in coFSubnetworks:
            coFDensities.append(statistics.calcEdgeDensity(network))

        lociDensities = []
        for network in lociSubN:
            lociDensities.append(statistics.calcEdgeDensity(network))

        statistics.overlappingHistogram(coFPopDensities, lociDensities)

        #print('P-val : ', pval)

main()