import networkx as nx
import matplotlib.pyplot as plt
import igraph as ig
import random
import json
import math

# data_dongtam = []
# comt = []
# size = []
# d = {}
# pos = {}
# with open('D:\Thesis\Dong Tam\covid19.json') as file:
#     data_dongtam = json.load(file)
# g = nx.Graph()
# for x in data_dongtam:
#     g.add_edge(x["fbid1"], x["fbid2"])
# print(len(g.nodes))
# pr = nx.pagerank(g)
# for y in g.nodes():
#     size.append(200*pr[y]*len(g.nodes))
#     d[y] = random.randint(1, 3)
# for node in g.nodes():
#     if d[node] % 2 == 0:
#         pos[node] = [random.randrange(2000*d[node], 2000*d[node] + 3000, 50), random.randrange(0, 2000, 50)]
#     else:
#         pos[node] = [random.randrange(-2000*d[node], -2000*d[node] + 1000, 50), random.randrange(0, 2000, 50)]
# plt.xlim(-2500,7500)
# plt.ylim(-500, 2100)
# nx.draw_networkx_edges(g, pos)
# nx.draw_networkx_nodes(g, pos, node_size= size, node_color=list(d.values()), with_labels = True)
# # nx.draw_networkx_labels(g, pos)

# # plt.savefig("simple_path.png") # save as png
# plt.show() # display


def show_graph(G, data):
    list_edge = G.get_edgelist()
    g = nx.Graph()
    for edge in list_edge:
        g.add_edge(edge[0], edge[1])
    # print(g.nodes())
    # print(g.edges())
    # nx.draw(g, with_labels = True) 
    d = {}
    e = {}
    dic = {}
    
    dic_ov = {}
    pos = {}
    label = {}
    lis_noOv = []
    list_Ov = []
    size = []
    pr = nx.pagerank(g)
    commt = len(data)
    row = math.ceil(commt / 3) + 1
    col = math.ceil(commt / (row - 1))
    xlen = 100000
    ylen = 100000
    width = xlen / col
    height = ylen / row
    for count, comm in enumerate(data):
        for x in comm:
            if x.index not in lis_noOv:
                lis_noOv.append(x.index)
                dic[x.index] = count

            else:
                list_Ov.append(x.index)
                dic_ov[x.index] = commt + dic[x.index]
                
    for node in g.nodes:
        size.append(300*pr[node]*len(g.nodes))
        label["node"] = node
        if node in list_Ov:
            d[node] = dic_ov[node]
            minwidth = width * (dic_ov[node] % col)
            maxwidth = minwidth + width
            # minheight = height * (dic_ov[node] / row)
            minheight = height * (row - 1)
            maxheight = minheight + height
            
            pos[node] = [random.uniform(minwidth + 100, maxwidth - 100), random.uniform(minheight + 100, maxheight - 100)]
            # pos[node] = [10*node + (-1)**dic[node]*2000, 10*node + (-1)**dic[node]*2000]
            # pos[node] = [random.randrange((-1)*1000, 1000, 10), random.randrange(0, 4000, 50)]
        else:
            d[node] = dic[node]
            minwidth = width * (dic[node] % col)
            maxwidth = minwidth + width
            minheight = height * math.floor(dic[node]/col)
            maxheight = minheight + height
            pos[node] = [random.uniform(minwidth + 50, maxwidth - 50), random.uniform(minheight + 50, maxheight - 50)]
            # pos[node] = [10*node + (-1)**dic[node]*1000, 10*node + dic[node]*1000]
            # if dic[node] % 2 == 0:
            #     pos[node] = [random.randrange(3000*dic[node] + 1000, 3000*dic[node] + 5000, 50), random.randrange(0, 4000, 50)]
            # else:
            #     pos[node] = [random.randrange(-2000*dic[node] - 5000, -2000*dic[node], 50), random.randrange(0, 4000, 50)]
    # for count, comm in enumerate(data):
    #     for x in comm:
    #         d[x.index] = count
    #         pos[x.index] = [random.uniform((-1)**count, (-1)**count + 1), random.uniform(count, count + 1)]
    
    # plt.xlim(-3000,7500)
    # plt.ylim(-4000, 4000)
    # print("Nodes overlap:")
    # for fbid in list_Ov:
    #     print(G.vs[fbid]["fbid"])
    nx.draw_networkx_edges(g, pos, width = 1.5)
    nx.draw_networkx_nodes(g, pos, node_size= size, node_color=list(d.values()), with_labels = True, linewidths = 1.0, edgecolors ='#008000')
    # nx.draw_networkx_labels(g, pos)

    # plt.savefig("simple_path.png") # save as png
    plt.show()