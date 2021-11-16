# Author: Katherina Cortes
# Date: October 30, 2021
# Purpose: Make network of specified top nodes by gene score with no edges between genes in same locus
#   Make graph of the network where:
#       1. gene score - size of node
#       2. loci of gene - color of node
#       3. weight of edge - darkness of edge


import networkx as nx
import matplotlib.pyplot as plt
import nxviz
from nxviz import annotate
from nxviz.plots import despine, aspect_equal


# @param nodes: list of nodes in network
# @param connections: dictionary of nodes and their edges
# @param lociLists: list of list of genes separated by loci
# @returns network: dictionary of genes of interest and edges, no edges between genes in same loci
def makeCrossLociNetwork(nodes, connections, lociLists):
    network = {}
    for n in nodes:
        network[n] = {}
        # get loci of gene
        for l in lociLists:
            if n in l:
                loci = l
        # make edges
        # no edges allowed between genes in same loci
        for edge in connections[n]:
            if edge in nodes:
                if edge not in loci:
                    network[n][edge] = connections[n][edge]
    return network


# @param network: dictionary of nodes and their edges
# @param lociLists: list of list of genes separated by loci
# @param geneVals: genes and their gene scores
# @returns g: networkx graph object with node, edge, edge weight information
def makeGraph(network, lociLists, geneVals):
    g = nx.Graph(name='Locus Gene Interactions Graph')
    # add nodes and their edges
    for node in network:
        g.add_node(node)
        g.nodes[node]['value'] = geneVals[node]
        g.nodes[node]['label'] = node
        if len(network[node]) > 0:
            for node2 in network[node]:
                weight = float(network[node][node2])
                g.add_weighted_edges_from([(node, node2, weight)])
                g.edges[node, node2]['weight'] = weight

    # add loci information for each node
    for n in network:
        for loci in range(len(lociLists)):
            if n in lociLists[loci]:
                g.nodes[n]['class'] = 'Loci ' + str(loci)
    return g


# make graph with specified genes
# 1. gene score - size of node
# 2. loci of gene - color of node
# 3. weight of edge - darkness of edge
# loci legend if genes are top overall regardless of loci
# make network between loci genes no edges between genes in same loci
# @param G: graph to display
# @param topGenes: bool where True gives loci legend
# @displays: nxviz circos graph
def visualizeGraph(G, topGenes):
    nxviz.circos(G, group_by='class', node_color_by='class', node_size_by='value', edge_alpha_by='weight')

    if topGenes:
        annotate.node_colormapping(G, color_by='class', legend_kwargs={'loc':'best', 'bbox_to_anchor':(0.1, 0.1, 0.1, 0.1)})
    # add legend if genes were picked for overall score instead of loci
    else:
        annotate.circos_group(G, group_by='class')
    annotate.circos_labels(G)
    aspect_equal()
    despine()
    plt.show()
    return


# @param geneAvg: dictionary of gene scores
# @param genes: list of genes of interest
# @param outFile: file name to write gene information to
# @param lociLists: list of list of genes separated by loci
# @outputs: file of gene loci gene Scores
def outputGeneScores(geneAvg, genes, outFile, lociLists):
    with open(outFile, 'w') as f:
        for gene in genes:
            for loci in range(len(lociLists)):
                if gene in lociLists[loci]:
                    break
            f.write(gene + '\t' + str(loci) + '\t' + str(geneAvg[gene]) +'\n')
    return
