import igraph as ig
import operator
from collections import Counter, defaultdict
import sys
import overlap
import show_graph
import numpy as np
import json
import time


def gn(origGraph):
	
	# initialize a list of removed edges that result in a split of the graph
	splits = []
	nClusters = 1
	G = origGraph.copy() 
	allCovers = {nClusters : ig.VertexCover(origGraph)}
	
	while G.es:
		# Calculate all edge betweennesses
		# TODO: only recalculate on relevant portions
		edge_betweennesses = G.edge_betweenness()

		# returns an the first index if there is a tie at max.
		max_index, _ = max(enumerate(edge_betweennesses), key=operator.itemgetter(1))

		# edge with the max betweenness
		edge = G.es[max_index].tuple
		print(edge)

		G.delete_edges(edge)

		if splitGraph(G, edge):

			# edge is a tuple, but we want a list of lists.
			splits += [list(edge)]
			comm = G.components().membership
			cover = get_cover(G, origGraph, comm)
			nClusters += 1
            # short circuit stuff would go here.
			allCovers[nClusters] = cover

	# vd = createDendrogram(origGraph, splits)
	# print(vd)

	# If we don't call this then as_clustering() fails. bugfix in development branch.
	# vd.optimal_count 
	return overlap.CrispOverlap(origGraph, allCovers)
	# return vd

def get_cover(G, OG, comm):
	coverDict = defaultdict(list)
	for i, community in enumerate(comm):
		coverDict[community].append(int(i))
	return ig.clustering.VertexCover(OG, clusters=list(coverDict.values()))

def splitGraph(G, edge):

	# return not G.edge_disjoint_paths(source=edge[0], target=edge[1])
	try:
		return not G.edge_disjoint_paths(source=edge[0], target=edge[1])
		# TODO: 
	except Exception as e:
		return False


def createDendrogram(G, splits):

	# To create a dendrogram, new merges have id of max id + 1
	n = len(splits) + 1
	merges = []

	mergeDict = {}

	while splits:
		# most recent split popped off
		edge = splits.pop()

		# Get the values the dendrogram wants for each vertex by finding
		# where merges have already happened.
		edge = [traverse(vertex, mergeDict) for vertex in edge]

		merges += [edge]

		# Update the dict to reflect a new merge.
		for vertex in edge:
			mergeDict[vertex] = n
		
		n += 1

	return ig.VertexDendrogram(G, merges)


def traverse(vertex, mergeDict):
	
	while vertex in mergeDict:
		vertex = mergeDict[vertex]
	return vertex

def run_demo():
	# f = open("D:\DATA\AA-DATN\DATN_GroupDetection\Source_Code\Data set\Data set\dolphins.txt", "r")
	# data = []
	# data_graph = []
	# for x in f:
	# 	x.replace("\n", " ")
	# 	data.append(x.split(" "))
	# for y in data:
	# 	data_graph.append((int(y[0]),int(y[1][:1])))
	# G = ig.Graph()
	# G.add_vertices(62)
	# G.add_edges(data_graph)
	G = ig.Graph().Famous("Zachary").as_undirected()
	result = gn(G)
	data_graph = result.pretty_print_cover(result.optimal_count, label='CONGO_index')
	show_graph.show_graph(G, data_graph)


def main():
	start_time = time.time()
	run_demo()
	run_time = time.time() - start_time
	print(run_time)

if __name__ == "__main__":
	main()
	