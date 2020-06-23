from collections import Counter, defaultdict
import itertools
import igraph as ig
import numpy as np
import operator
import logging
import argparse
# import overlap1

import overlap
import show_graph
import json

def congo(OG, h=2):

    G = OG.copy()
    # Just in case the original graph is disconnected
    if not G.is_connected():
        raise RuntimeError("Congo only makes sense for connected graphs.")
    
    # initializing attributes of copied graph
    G.vs['CONGO_orig'] = [i.index for i in OG.vs]
    G.es['eb'] = 0 #khởi tạo giá trị edge betweenness cho các cạnh 
    G.vs['pb'] = [{uw : 0 for uw in itertools.combinations(G.neighbors(vertex), 2)} for vertex in G.vs] 
    
    do_initial_betweenness(G, h)
    nClusters = 1

    # The first cover is simply the entire connected graph.
    allCovers = {nClusters : ig.VertexCover(OG)}
    while G.es:
        
        # get the edge with the max edge betweenness, and its betweenness.
        maxEdge, maxEb = max(enumerate(G.es['eb']), key=operator.itemgetter(1))
        G.vs['vb'] = G.betweenness(cutoff=h)

        # TODO check if I need to multiply by 2
        vInteresting = [i for i, b in enumerate(G.vs['vb']) if b > maxEb]

        splitInstr = max_split_betweenness(G, vInteresting)

        # split if max split betweenness > max edge betweenness
        if splitInstr is None or splitInstr[0] <= maxEb:
            split = delete_edge(G, maxEdge, h)
        else:
            split = split_vertex(G, splitInstr[1], splitInstr[2], h)



        if split:
            # there must be a new community
            comm = G.components().membership
            cover = get_cover(G, OG, comm)
            nClusters += 1
            
            # short circuit stuff would go here.
            allCovers[nClusters] = cover
            if nClusters > 2:
                break
    for i in range(len(allCovers)):
        print(allCovers[i+1])
    return overlap.CrispOverlap(OG, allCovers)


def delete_edge(G, edge, h):

    tup = G.es[edge].tuple

    logging.info("Deleted: %s", tup)

    neighborhood = get_neighborhood_edge(G, tup, h)
    
    do_local_betweenness(G, neighborhood, h, operator.neg)
    G.delete_edges(edge)
    # fix_betweennesses(G)
    do_local_betweenness(G, neighborhood, h, operator.pos)
    return check_for_split(G, tup)



def fix_pair_betweennesses(G):
    
    for v in G.vs:
        toDel = []
        neededPairs = {uw for uw in itertools.combinations(G.neighbors(v), 2)}
        for pair in v['pb']:
            if pair not in neededPairs:
                toDel.append(pair)
        for d in toDel:
            del v['pb'][d]
        for pair in neededPairs:
            if pair not in v['pb']:
                v['pb'][pair] = 0


def fix_edge_betweennesses(G):
    
    for e in G.es:
        if e['eb'] is None:
            e['eb'] = 0


def fix_betweennesses(G):
    
    fix_pair_betweennesses(G)
    fix_edge_betweennesses(G)


def split_vertex(G, vToSplit, instr, h):
    
    neighborhood = get_neighborhood_vertex(G, vToSplit, h)
    do_local_betweenness(G, neighborhood, h, operator.neg)
    new_index = G.vcount()
    G.add_vertex()
    G.vs[new_index]['CONGO_orig'] = G.vs[vToSplit]['CONGO_orig']
    G.vs[new_index]['pb'] = {uw : 0 for uw in itertools.combinations(G.neighbors(vToSplit), 2)}

    # adding all relevant edges to new vertex, deleting from old one.
    toAdd = list(zip(itertools.repeat(new_index), instr[0]))
    toDelete = list(zip(itertools.repeat(vToSplit), instr[0]))
    G.add_edges(toAdd)
    G.delete_edges(toDelete)
    neighborhood.append(new_index)
    fix_betweennesses(G)
    logging.info("split: %d, %s", vToSplit, instr)
    do_local_betweenness(G, neighborhood, h, operator.pos)
    # check if the two new vertices are disconnected.
    return check_for_split(G, (vToSplit, new_index))


def max_split_betweenness(G, vInteresting):
    
    maxSplitBetweenness = 0
    vToSplit = None
    
    for v in vInteresting:
        clique = create_clique(G, v, G.vs['pb'][v])
        if clique.size < 4:
            continue

        # initialize a list on how we will map the neighbors to the collapsing matrix
        vMap = [[ve] for ve in G.neighbors(v)]

        while clique.size > 4:
            i,j,clique = reduce_matrix(clique)
            vMap[i] += vMap.pop(j)

        if clique[0,1] > maxSplitBetweenness:
            maxSplitBetweenness = clique[0,1]
            vToSplit = v
            splitInstructions = vMap
    if vToSplit is None:
        return None

    return maxSplitBetweenness, vToSplit, splitInstructions


def do_initial_betweenness(G, h):
    
    # all_pairs_shortest_paths = []
    pathCounts = Counter()
    apsp = []
    for ver in G.vs:
        #tập các đỉnh kề và đỉnh có giá trị đường ngắn nhất đến đỉnh ver.index bằng 2
        neighborhood = get_neighborhood_vertex(G, ver, h)
        neighborhood.remove(ver.index)

        #tập đường đi ngắn nhất từ đỉnh ver đến các đỉnh còn lại
        s_s_shortest_paths = G.get_all_shortest_paths(ver, to=neighborhood)
        # all_pairs_shortest_paths += s_s_shortest_paths

        #tính số đường đi ngắn nhất giữa 2 đỉnh
        for path in s_s_shortest_paths:
            pathCounts[(path[0], path[-1])] += 1
            if len(path) <= h + 1:
                apsp.append(path)

    for path in apsp:
        update_betweenness(G, path, pathCounts[(path[0], path[-1])], operator.pos)


def do_local_betweenness(G, neighborhood, h, op=operator.pos):
    
    # all_pairs_shortest_paths = []
    pathCounts = Counter()
    neighSet = set(neighborhood)
    neighSize = len(neighborhood)
    apsp = []
    for i, v in enumerate(neighborhood):
        s_s_shortest_paths = G.get_all_shortest_paths(v, to=neighborhood)
    
        for path in s_s_shortest_paths:
            if len(neighSet | set(path)) == neighSize:
                pathCounts[(path[0], path[-1])] += 1 
                if len(path) <= h + 1:
                    apsp.append(path)
    for path in apsp:  
        update_betweenness(G, path, pathCounts[(path[0], path[-1])], op)


def update_betweenness(G, path, count, op):
    
    weight = op(1./count)
    pos = 0
    while pos < len(path) - 2:
        G.vs[path[pos + 1]]['pb'][order_tuple((path[pos], path[pos + 2]))] += weight
        G.es[G.get_eid(path[pos], path[pos + 1])]['eb'] += weight
        pos += 1
    if pos < len(path) - 1:
        G.es[G.get_eid(path[pos], path[pos + 1])]['eb'] += weight


def get_cover(G, OG, comm):
    
    coverDict = defaultdict(list)
    for i, community in enumerate(comm):
        vt = int(G.vs[i]['CONGO_orig'])
        if vt not in coverDict[community]:
            coverDict[community].append(vt)
    return ig.clustering.VertexCover(OG, clusters=list(coverDict.values()))


def get_neighborhood_vertex(G, v, h):
    
    return G.neighborhood(v, order=h)


def get_neighborhood_edge(G, e, h):
    
    neigh = set(G.neighborhood(e[0], order=h-1))
    neigh.update(G.neighborhood(e[1], order=h-1))
    return list(neigh)


def order_tuple(toOrder):
    if toOrder[0] <= toOrder[1]:
        return toOrder
    return (toOrder[1], toOrder[0])


def create_clique(G, v, pb):
    
    neighbors = G.neighbors(v)

    # map each neighbor to its index in the adjacency matrix
    mapping = {neigh : i for i, neigh in enumerate(neighbors)}
    n = len(neighbors)

    # Can use ints instead: (dtype=int). Only works if we use matrix_min
    # instead of mat_min.
    clique = np.matrix(np.zeros((n, n)))
    for uw, score in pb.items():
        clique[mapping[uw[0]], mapping[uw[1]]] = score
        clique[mapping[uw[1]], mapping[uw[0]]] = score

    # Ignore any self loops if they're there. If not, this line
    # does nothing and can be removed.
    np.fill_diagonal(clique, 0)
    return clique


def reduce_matrix(M):
    
    i,j = mat_min(M)
    M[i,:] = M[j,:] + M[i,:]

    # delete the jth row
    M = np.delete(M, (j), axis=0)

    # similarly with the columns
    M[:,i] = M[:,j] + M[:,i]
    M = np.delete(M, (j), axis=1)
    np.fill_diagonal(M,0)
    return i,j,M


def check_for_split(G, edge):
    
    try:
        return not G.edge_disjoint_paths(source=edge[0], target=edge[1])
        # TODO: 
    except Exception as e:
        return False


def mat_min(M):
    
    np.fill_diagonal(M, float('inf'))

    i,j = np.unravel_index(M.argmin(), M.shape)
    np.fill_diagonal(M,0)
    return i, j



def run_demo():

    # data = []
    # de = []
    # with open('D:\DATA\AA-Luan_Van_ThS\Thesis\Dong Tam\data_dongtam.txt') as json_file:
    #     database = json.load(json_file)
    # for x in database:
    #     datatest = (x["fbid1"], x["fbid2"])
    #     data.append(datatest)
    # G = ig.Graph()
    # G.add_vertices(37017)
    # G.add_edges(data)
    # comn = G.components()
    # de_cn = []
    # for cn in comn:
    #     if len(cn) < 35000:
    #         for i in cn:
    #             de_cn.append(i)
    # de_cn.sort(reverse=True)
    # for j in de_cn:
    #     G.delete_vertices(j)    
    # vb = G.betweenness()
    # for x, y in enumerate(vb):
    #     if y < 1.0:
    #         de.append(x)
    # de.sort(reverse=True)
    # for i in range(len(de)):
    #     G.delete_vertices(de[i])
    

    # data = []
    # de = []
    # data_fbid = []
    # database = []
    # with open('D:\DATA\AA-Luan_Van_ThS\Thesis\Dong Tam\data_Graph.json') as json_file:
    #     database = json.load(json_file)
    # with open('D:\DATA\AA-Luan_Van_ThS\Thesis\Dong Tam\data_vertex_DongTam.json') as json_file:
    #     data_fbid = json.load(json_file)
    # print(len(data_fbid))
    # for x in database:
    #     datatest = (x["fbid1"], x["fbid2"])
    #     data.append(datatest)
    # G = ig.Graph()
    # G.add_vertices(len(data_fbid))
    # G.add_edges(data)
    # G.vs["fbid"] = data_fbid
    # # for y in G.vs:
    # #     print(y["fbid"])
    # comn = G.components()
    # de_cn = []
    # for cn in comn:
    #     if len(cn) < 18000:
    #         for i in cn:
    #             de_cn.append(i)
    # de_cn.sort(reverse=True)
    # for j in de_cn:
    #     G.delete_vertices(j)  
    # print(G.get_edgelist()) 
    # vb = G.betweenness()
    # for x, y in enumerate(vb):
    #     if y == 0.0:
    #         de.append(x)
    # de.sort(reverse=True)
    # for i in range(len(de)):
    #     G.delete_vertices(de[i])


    
    data = []
    de = []
    database = []
    data_fbid = []
    with open('D:\DATA\AA-Luan_Van_ThS\Thesis\Dong Tam\Covid-19\covid-19.json') as json_file:
        database = json.load(json_file)
    with open('D:\DATA\AA-Luan_Van_ThS\Thesis\Dong Tam\Covid-19\data_vertex_Covid-19.json') as json_file:
        data_fbid = json.load(json_file)
    for x in database:
        datatest = (x["fbid1"], x["fbid2"])
        data.append(datatest)
    G = ig.Graph()
    print(len(data_fbid))
    G.add_vertices(len(data_fbid))
    G.add_edges(data)
    G.vs["fbid"] = data_fbid
    comn = G.components()
    de_cn = []
    for cn in comn:
        if len(cn) < 3400:
            for i in cn:
                de_cn.append(i)
    de_cn.sort(reverse=True)
    for j in de_cn:
        G.delete_vertices(j) 


    # database = []
    # f = open("D:\DATA\AA-DATN\DATN_GroupDetection\Source_Code\Data set\Data set\dolphins.txt", "r")
    # for x in f:
    #     y = x.split('\n')
    #     z = y[0].split(' ')
    #     t = (int(z[0])-1, int(z[1])-1)
    #     database.append(t)
    # G = ig.Graph()
    # G.add_vertices(62)
    # G.add_edges(database)
    # print(G)


    # G = ig.Graph().Famous("Zachary").as_undirected()
    
    result = congo(G, 2)
    data_graph = result.pretty_print_cover(G, result.optimal_count)
    show_graph.show_graph(G, data_graph)

def main():
    run_demo()


if __name__ == "__main__":
    main()

