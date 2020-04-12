# -*- coding: utf-8 -*-
import os
import glob
import networkx as nx
import sys
import matplotlib.pyplot as plt
import random
import numpy as np

XX = ["elegans", "coli", "drosop", "arab", "human", "inter_H-Y", "mus_musculus", "rattus_norvegicus", "yeast"]
#XX=["rattus_norvegicus"]

if not os.path.exists('info-networks'):
    os.makedirs('info-networks')

cwd = os.getcwd() 

for i in XX:
    G = nx.read_adjlist(str(cwd) + '/data/edge_lists_ints/'+str(i)+'.txt')
    F = open(str(cwd) + "/info-networks/"+str(i)+"-info.txt", "w")
    F.write(str(i) + " Network" + str("\n"))
    F.write("Nodes = "+str(nx.number_of_nodes(G))+'\n')
    F.write('Edges = '+str(nx.number_of_edges(G))+'\n')
    #G.remove_edges_from(G.selfloop_edges())
    G = nx.to_undirected(G)
    F.write('Density = '+str(nx.density(G))+'\n') 
    F.write("Av. Clustering = "+str(nx.average_clustering(G))+"\n")
    F.write("Num of Connected Components = " + str(nx.number_connected_components(G))+"\n")
    avd = 0.0
    for s in list(G.degree()):
        avd = avd + s[1]
    F.write("Av. Degree = "+str(float(avd/len(list(G.degree()))))+"\n")
    if nx.is_connected(G)==True:
        F.write("Diameter = "+str(nx.diameter(G))+"\n")
        F.write("Av. Shortest Path Length = " + str(nx.average_shortest_path_length(G))+"\n")
        av_cc=0.0
        for s in list(nx.closeness_centrality(G)):
            av_cc+=s[1]
        F.write("Closeness Centrality = "+str(float(av_cc/len(list(nx.closeness_centrality(G)))))+"\n")
        av_cb=0.0
        for s in list(nx.betweenness_centrality(G)):
            av_cb+=s[1]
        F.write("Betweenness Centrality = "+str(float(av_cb/len(list(nx.betweenness_centrality(G)))))+"\n")

    F.close()

    lst = nx.degree_histogram(G)
    #lst.pop(0)

    y=[]
    x=[]
    for p in range(0, len(lst)):
        if lst[p]!=0:
            y.append(lst[p])
            x.append(p)

    fig = plt.figure(figsize=(6, 4)) #ANCHO/LARGO
    plt.subplot(111)
    plt.scatter(x=np.log(np.array(x)), y=np.log(np.array(y)), marker='o', color='red', s=20)
    ##plt.ylim(-0.05, 1.05)
    plt.xlabel("log-Degree")
    plt.ylabel("log-Frequency")
    plt.title('Degree Distribution - '+str(i))
    plt.tight_layout()
    plt.savefig(cwd + '/info-networks/'+str(i)+'-dg.png', format='png')


