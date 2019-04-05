


import sys
from itertools import *
from pprint import pprint
import cplex
from math import *
import networkx as nx
from cplex.exceptions import CplexSolverError
import pprint
from collections import deque
import matplotlib.pyplot as plt
import heapq as hq
import time
from random import choice, shuffle
import numpy as np



matfile = "/home/igorm/workspace/recomb_2019_good_simulation/igor_output/matrix_1/nick_matrix_1_case_1.for_ilp.vcf"


matfile = sys.argv[1]
outfile = sys.argv[2]


with open(matfile) as f:
    matdata = f.readlines()
matdata = map(lambda x: x.strip().split(" "), matdata)

mutation_matrix = []
for line in matdata:
    mutation_matrix.append(map(int, line))


#mut_trans = []
#for i in range(len(mutation_matrix[0])):
#    mut_trans.append([x[i] for x in mutation_matrix])
#mutation_matrix = mut_trans

columns = []
nrMutations = len(mutation_matrix) # number of rows - number of mutations
nrCells = len(mutation_matrix[0]) # number of cells

zeroGraph = nx.Graph()


for c2 in range(nrCells):
    c_column = [x[c2] for x in mutation_matrix]
    c_muts = []
    for order, elem in enumerate(c_column):
        if elem == 1:
            c_muts.append(order)
    for i in range(len(c_muts)):
        for j in range(len(c_muts)):
            if i == j:
                continue
            m1_cells = []
            m2_cells = []
            for order1, elem1 in enumerate(mutation_matrix[c_muts[i]]):
                if elem1 == 1 and mutation_matrix[c_muts[j]][order1] == 0:
                    m1_cells.append(order1)
            for order2, elem2 in enumerate(mutation_matrix[c_muts[j]]):
                if elem2 == 1 and mutation_matrix[c_muts[i]][order2] == 0:
                    m2_cells.append(order2)
            m1 = c_muts[i]
            m2 = c_muts[j]
            # enumerate pairs of non-edges - they occur in the same W
            for c1 in m1_cells:
                for c3 in m2_cells:
                    zeroGraph.add_edge(nrCells * m2 + c1, nrCells * m1 + c3)

zeroGraph2 = zeroGraph.copy()
zeroGraph3 = zeroGraph.copy()



##############################

def lexBFS(graph):
    graph_am = nx.adjacency_matrix(graph).todense()
    #print graph_am
    graph_nodes = list(graph.nodes())
    local_to_global = {}
    global_to_local = {}
    for i, node in enumerate(graph_nodes):
        local_to_global[i] = node
        global_to_local[node] = i
    graph = nx.relabel_nodes(graph, global_to_local)
    local_graph_nodes = sorted(graph.nodes())
    lists = map(lambda x: list(nx.neighbors(graph, x)), local_graph_nodes)
    n = len(graph_nodes)
    order = [0 for i in range(n)]
    weights = [0 for i in range(n)]
    for i in range(n - 1, -1, -1):
        M = max(weights)
        j = [xxx for xxx, mmm in enumerate(weights) if mmm == M][0]
        order[i] = j
        weights[j] = -1
        for u in lists[j]:
            if weights[u] != -1:
                weights[u] += 1
    #print order, "MMMMMMMMMMM"
    neworder = []
    #print order, "order"
    for o in order:
        neworder.append(local_to_global[o])
    return neworder
###############################


start = time.time()

vertexCover = set()

connDeque = deque(nx.connected_components(zeroGraph))

while connDeque:
    #print map(len, list(connDeque))
    conn = connDeque.popleft()
    conn = nx.Graph(nx.subgraph(zeroGraph, conn))
    if len(conn.nodes()) <= 1:
        continue
    if len(conn.nodes()) == 2:
        vertexCover.add(list(conn.nodes())[0])
        continue
    #print conn

    while True:
        lex = lexBFS(conn)
        #print lex
        order = np.array(list(reversed(lex)))

        absorbing = set()
        for node1 in order:
            for node2 in nx.neighbors(conn, node1):
                v1 = set([x for x in nx.neighbors(conn, node2) if x not in absorbing]) - set([node1])
                v2 = set([x for x in nx.neighbors(conn, node1) if x not in absorbing]) - set([node2])
                if v1 & v2 == v2:
                    absorbing.add(node2)
        if not absorbing:
            break
        else:
            #print len(absorbing), "absorbing"
            for node in absorbing:
                vertexCover.add(node)
                conn.remove_node(node)
    if len(conn.edges()) == 0:
        pass
    else:
        # find a vertex with maximum degree and remove it
        max_degree_vertex = max(dict(conn.degree()).items(), key=lambda x: x[1])[0]
        vertexCover.add(max_degree_vertex)
        conn.remove_node(max_degree_vertex)
        for conn2 in nx.connected_components(conn):
            if len(conn2) > 1:
                connDeque.append(conn2)



# check
edges1 = set()
for node in vertexCover:
    neigh = nx.neighbors(zeroGraph, node)
    for n in neigh:
        edges1.add(tuple(sorted([n, node])))

edges2 = set()
for x, y in zeroGraph.edges():
    edges2.add(tuple(sorted([x, y])))





for vertex in vertexCover:
    row = vertex / nrCells
    col = vertex - row * nrCells
    mutation_matrix[row][col] = 1 - mutation_matrix[row][col]



with open(outfile, "w") as f:
    for line in mutation_matrix:
        f.write("%s\n" % " ".join(map(str, line)))






