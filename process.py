# -*- coding: utf-8 -*-
#!/usr/bin/env python3
import networkx as nx
import pandas as pd
import sys
import os
import numpy as np

from sklearn.manifold import TSNE
from matplotlib import pyplot as plt
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler


############################################################
### COMMUNITIES DETECTION
############################################################

cwd = os.getcwd() 

XX = ["rattus_norvegicus", "elegans", "coli", "drosop", "arab", "human", "inter_H-Y", "mus_musculus", "rattus_norvegicus", "yeast"]

if not os.path.exists('net-communities'):
    os.makedirs('net-communities')

for i in XX:
    G = nx.read_adjlist(str(cwd) + '/data/edge_lists_ints/' + i + '.txt')

    m = nx.community.greedy_modularity_communities(G)
    c = list(m)

    J = open(cwd+'/net-communities/com-fastgreedy-'+i+'_NX.txt', 'w') 
    J.write ('# '+str(len(c))+' communities; '+str(G.number_of_nodes())+' elements \n')
    for q in range(len(c)):
	    J.write ('C'+str(q+1)+'-'+str(len(c[q])))
	    for w in c[q]:
		    J.write(' '+str(w))
	    J.write('\n')
    J.close()

############################################################
### NORMALIZACIÓN DE LOS EMBEDDINGS DE DEEPWALK
############################################################

#X = ["elegans" "coli" "drosop" "arab" "human" "inter_H-Y" "mus_musculus" "rattus_norvegicus" "yeast"]

cwd = os.getcwd() 

if not os.path.exists('dw-emb-norm'):
    os.makedirs('dw-emb-norm')

for i in XX:
	fname= cwd + '/deepwalk-embeddings/w10-80walks-40l-d128-' + str(i) + '.txt'
	data = pd.read_csv(fname, header=0, delimiter=' ')
	data_scale = (data - data.mean()) / data.std()
	name = 'dw-emb-norm/w10-80walks-40l-d128-' + str(i) + '_norm.txt'
	pd.DataFrame.to_csv(data_scale, path_or_buf=name, sep=' ', float_format='%.8f', columns=None, header=False, index=True, index_label=True)

############################################################
### 2D T-SNE SOBRE LOS EMBEDDINGS NORMALIZADOS DE DEEPWALK
############################################################

if not os.path.exists('tsne-2d'):
    os.makedirs('tsne-2d')

#XX = ["elegans" "coli" "drosop" "arab" "human" "inter_H-Y" "mus_musculus" "rattus_norvegicus" "yeast"]

cwd = os.getcwd() 

for i in XX:
	fname = cwd + '/dw-emb-norm/w10-80walks-40l-d128-'+ str(i)+'_norm.txt'
	data = pd.read_csv(fname, delimiter=' ', header=None)
	df = data.values
	X = df[::,1:]
	tsne = TSNE(n_components=2, random_state=0)
	X_2d = tsne.fit_transform(X)
	XD = pd.DataFrame(X_2d)
	result = pd.concat([data[0], XD], axis=1, sort=False)
	np.savetxt(cwd + '/tsne-2d/tsne-w10-80walks-40l-d128-'+str(i), result.values, fmt='%f')
	plt.figure(figsize=(6, 5))
	for m in range(len(X_2d)):
		plt.scatter(X_2d[m, 0], X_2d[m, 1], s=2, c='red', alpha=0.5)
	plt.title("d128 "+str(i))
	plt.savefig(str(cwd) + '/tsne-2d/2d-tsne-' + str(i) + '-d128.png', format='png')

############################################################
### DBSCAN SOBRE LOS EMBEDDINGS 2D DE T-SNE
############################################################

if not os.path.exists('dbscan'):
    os.makedirs('dbscan')

#XX = ["elegans" "coli" "drosop" "arab" "human" "inter_H-Y" "mus_musculus" "rattus_norvegicus" "yeast"]

numC = []

cwd = os.getcwd() 

for i in XX:
    fname = cwd + '/tsne-2d/tsne-w10-80walks-40l-d128-' + str(i)
    data = pd.read_csv(fname, delimiter=' ', header=None)
    df = data.values
    X = df[::, 1:]
    X = StandardScaler().fit_transform(X)
    # Compute DBSCAN
    db = DBSCAN(eps=0.1, min_samples=10).fit(X)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
    dff = pd.DataFrame(labels)
    result = pd.concat([data[0], dff], axis=1, sort=False)
    result.columns = ['a', 'b']
    j = result.groupby(['b']).groups
    J = open(cwd + '/dbscan/c-dbscan-tsne-w10-80walks-40l-d128-' + str(i) + '.txt', 'w') 
    for q in list(j.keys()):
        J.write ('C' + str(q)+ '-' +str(len(list(j[q]))))
        for w in j[q]:
            J.write(' ' + str(w))
        J.write('\n')
    J.close()
# Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)
    numC.append(n_clusters_)    
    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each)
    for each in np.linspace(0, 1, len(unique_labels))]
    plt.figure(figsize=(6, 5))
    for u, col in zip(unique_labels, colors):
		if u == -1:
			# Black used for noise.
			col = [0, 0, 0, 1]
		class_member_mask = (labels == u)
		xy = X[class_member_mask & core_samples_mask]
		plt.scatter(xy[:, 0], xy[:, 1], s=14, c=tuple(col), alpha=0.5)#, edgecolors='k')
		xy = X[class_member_mask & ~core_samples_mask]
		plt.scatter(xy[:, 0], xy[:, 1], s=14, c=tuple(col), alpha=0.5)#, edgecolors='k')
    plt.title('Estimated number of clusters: %d' % n_clusters_)
    plt.savefig(str(cwd) + '/dbscan/c-dbscan-tsne-w10-80walks-40l-d128-' + str(i) + '.png', format='png')

############################################################
### MAP INTEGER NODES TO STRING NODES
############################################################

if not os.path.exists('proteins'):
    os.makedirs('proteins')

cwd = os.getcwd() 

#XX = ["elegans" "coli" "drosop" "arab" "human" "inter_H-Y" "mus_musculus" "rattus_norvegicus" "yeast"]

cn = 0

for i in XX:
    G = nx.read_adjlist(str(cwd) + '/data/edge_lists_ints/' + str(i) + '.txt')
    J = open(str(cwd) + '/data/edge_lists_ints/'+str(i)+'.txt', 'r') 
    K = open(str(cwd) + '/data/edge_lists_strings/'+str(i)+'.txt', 'r') 
    s = []
    l = 1
    count = 0
    while len(s) < G.number_of_nodes() :
        l = J.readline()
        m = l.strip().split(' ')
        k = K.readline()
        n = k.strip().split(' ')
        ee = [m[0], n[0]]
        rr = [m[1], n[1]]
        if count==0:
            s.append(ee)
            s.append(rr)
            count = 1
        if len(s) > 0:
            if ee not in s:
                s.append(ee)

            if rr not in s:
                s.append(rr)
    J.close()
    K.close()

    O = open(str(cwd) + '/dbscan/c-dbscan-tsne-w10-80walks-40l-d128-'+str(i)+'.txt', 'r') 
    T = open(str(cwd) + "/proteins/"+str(i)+"-proteins.txt", "w")

    count=0
    while count < numC[cn]:
        l = O.readline()
        m = l.strip().split(' ')
        for element in m:
            ff = False
            count1=0
            while ff == False and count1<G.number_of_nodes():
                rrr = s[count1]
                if(element==rrr[0]):
                    ff=True
                    T.write(rrr[1]+' ')
                count1+=1
        T.write('\n')
        count+=1 
    O.close()
    T.close()
    cn=cn+1





