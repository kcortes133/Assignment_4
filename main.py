# Author: Katherina Cortes
# Date: October 30, 2021
# Purpose: Take a tab delimited .gmz file of gene sets, a given tab-delimited STRING file
#   create a sub network of the gene interactions from the input file using the STRING file
#   get statistical significance

import argparse, random, time
import networkCreation, fileParsing, statistics, geneScoring, \
    networkVisualization, geneticAlgorithm, outputFiles

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
parser.add_argument('--topGenes', type=bool, default=False, help='graph only top n genes regardless of loci')
parser.add_argument('--numGenes', type=int, default=3, help='number of genes from each loci or total genes'
                                                                 ' if topGenes is true')

args = parser.parse_args()

def main():
    start = time.time()
    random.seed(5)
    visualize = True

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


    print(time.time() - start)

    # get top numGenes from each loci
    # make network with genes
    if visualize:
        if not args.topGenes:
            genes = geneScoring.getTopLociGenes(geneAvg, lociLists, args.numGenes)
            visualNetwork = networkVisualization.makeCrossLociNetwork(genes, network, lociLists)

        # get top numGenes regardless of loci
        # make network with genes
        if args.topGenes:
            genes = networkSorted[:args.numGenes]
            visualNetwork = networkVisualization.makeCrossLociNetwork(genes, network, lociLists)

        # write to output file the genes loci and gene score of genes in the network
        networkVisualization.outputGeneScores(geneAvg, genes, args.geneOutFile, lociLists)
        # make graph with specified genes
        # 1. gene score - size of node
        # 2. loci of gene - color of node
        # 3. weight of edge - darkness of edge
        # make network between loci genes no edges between genes in same loci
        graph = networkVisualization.makeGraph(visualNetwork, lociLists, geneAvg)
        networkVisualization.visualizeGraph(graph, args.topGenes)

    if args.calcPVal:
        numBins = args.numBins
        # make bins for coFunctional subnetwork creation
        qNetworkBins = networkCreation.makeQuantileBins(interactions, numBins)
        fNetworkBins = networkCreation.makeFixedBins(interactions, numBins)
        # make coFunctional random subnetworks
        # need 1000 populations where populations are the 5000 networks
        coFPopDensities = []
        for i in range(1000):
            coFSubnetworks = networkCreation.makeCoFSubnetworks(interactions, qNetworkBins, newPop)
            popD = 0

            for subCoF in coFSubnetworks:
                popD += statistics.calcEdgeDensityW(subCoF)
            coFPopDensities.append(popD/len(coFSubnetworks))

        # calculate the avg of each population -> makeCoFSubnetorks is one population?
        # calculate the pvalue
        # probability edges using cof distribution is greater than avg of loci edged divided by # of random networks

        pval = statistics.empiricalPVal(newPop, coFPopDensities)

        # make a graph showing the edge density distributions
        coFDensities = []
        for network in coFSubnetworks:
            coFDensities.append(statistics.calcEdgeDensityW(network))

        lociDensities = []
        for network in newPop:
            lociDensities.append(statistics.calcEdgeDensityW(network))

        statistics.overlappingHistogram(lociDensities, coFPopDensities)

        print('P-val : ', pval)
    print(time.time() - start)
    outputFiles.outputNetworks(2, newPop, 10)
    outputFiles.outputGeneScoresinLoci(geneAvg, lociLists)

main()