#!/usr/bin/env python
# coding: utf-8

import networkx as nx
import matplotlib.pyplot as plt
import json
import heapq
import math

#this function tests to see if all nodes are connected. If they are not, the rest of the 
#program won't work
def conn_graph(g):
	m = nx.is_connected(g)
	if m == False:
		print 'Graph not connected'
		return 0
	else:
		return 1

g = nx.Graph()	#initialize an empty graph for all interactions

in_file = open('rocky_intertidal_interactions.csv', 'r')	#read in the data file
out_file = open('rocky_intertidal_results.txt', 'w')
desc = 'rocky intertidal' 
delim = ','
	 
all_names = []
in_file.next()
for line in in_file:
	line = line.strip()
	row = line.split(delim)
	source_taxon = row[0]
	target_taxon = row[2]
	if source_taxon in all_names:
		pass
	else:
		all_names.append(source_taxon)
	if target_taxon in all_names:
		pass
	else:
		all_names.append(target_taxon)
print all_names	#I have a lot of print statements because I like to see what's going on
print len(all_names)
g.add_nodes_from(all_names, name=all_names)   #add all nodes

for name in all_names:	#add edges to all interaction graph
	in_file.seek(0)
	in_file.next()
	for line in in_file:
		line = line.strip()
		row = line.split(delim)
		source_taxon = row[0]
		interaction_type = row[1]
		target_taxon = row[2]
		if name == source_taxon:
			g.add_edge(name,target_taxon, type='eats')
		else:
			pass
nn = g.number_of_nodes()
ne = g.number_of_edges()
print ne
print nn

#determine if all nodes are connected in one graph for all interactions
j = conn_graph(g)
if j == 1:	#if all nodes are connected, proceed with calculations
	all_nodes = g.nodes()
	all_edges = g.edges()
	deg_c = nx.degree_centrality(g)
	bet_c = nx.betweenness_centrality(g)
	clo_c = nx.closeness_centrality(g)
	eig_c = nx.eigenvector_centrality_numpy(g)
	edg_b = nx.edge_betweenness_centrality(g,normalized=True)
	kat_c = nx.katz_centrality(g,alpha=0.01,beta=1.0,max_iter=1000,normalized=True)
	cent = nx.center(g)
	ecc = nx.eccentricity(g)
	peri = nx.periphery(g)
	dom_s = nx.dominating_set(g)
	clo_v = nx.closeness_vitality(g)
else:	#if all nodes are NOT connected, pick the largest graph and proceed with calculations
	print 'graph not connected'

out_file.write('taxon' + '\t' + 'system' + '\t' + 'degree_centrality' + '\t' + 'betweenness_centrality' + '\t' + 'closeness_centrality' + '\t' + 'eigenvector_centrality' + '\t' + 'katz_centrality' + '\t' + 'eccentricity' + '\t' + 'closeness_vitality' + '\t' + 'center' + '\t' + 'periphery' + '\t' + 'dominating_set' + '\n')
for n in all_nodes:
	if n in cent:
		cn = '1'
	else:
		cn = '0'
	if n in peri:
		pn = '1'
	else:
		pn = '0'
	if n in dom_s:
		dn = '1'
	else:
		dn = '0'

	out_file.write(n + '\t' + desc + '\t' + str(deg_c[n]) + '\t' + str(bet_c[n]) + '\t' + str(clo_c[n]) + '\t' + str(eig_c[n]) + '\t' + str(kat_c[n]) + '\t' + str(ecc[n]) + '\t' + str(clo_v[n]) + '\t' + cn + '\t' + pn + '\t' + dn + '\n')


nx.draw_spring(g)
plt.show()

#Here are some of the not-so-useful calculations that I stopped using.
"""
out_file.write('squares' + '\n')
for a in g.nodes_iter():	#calculate number of squares
	counter = 0
	for b in g.neighbors(a):
		for c in g.neighbors(b):
			for d in g.neighbors(c):
				if g.has_edge(a,d):
					#print a,b,c,d
					counter = counter + 1
				else:
					pass
	out_file.write(a + '\t' + str(counter) + '\n')

out_file.write('\n' + 'pentagons' + '\n')
for a in g.nodes_iter():	#calculate number of pentagons
	counter = 0
	for b in g.neighbors(a):
		for c in g.neighbors(b):
			for d in g.neighbors(c):
				for e in g.neighbors(d):
					if g.has_edge(a,e):
						#print a,b,c,d,e
						counter = counter + 1
					else:
						pass
	out_file.write(a + '\t' + str(counter) + '\n')

out_file.write('\n' + 'hexagons' + '\n')
for a in g.nodes_iter():	#calculate number of hexagons
	counter = 0
	for b in g.neighbors(a):
		for c in g.neighbors(b):
			for d in g.neighbors(c):
				for e in g.neighbors(d):
					for f in g.neighbors(e):
						if g.has_edge(a,f):
							#print a,b,c,d,e,f
							counter = counter + 1
						else:
							pass
	out_file.write(a + '\t' + str(counter) + '\n')

out_file.write('\n' + 'heptagons' + '\n')
for a in g.nodes_iter():	#calculate number of heptagons
	counter = 0
	for b in g.neighbors(a):
		for c in g.neighbors(b):
			for d in g.neighbors(c):
				for e in g.neighbors(d):
					for f in g.neighbors(e):
						for h in g.neighbors(f):
							if g.has_edge(a,h):
								#print a,b,c,d,e,f,h
								counter = counter + 1
							else:
								pass
	out_file.write(a + '\t' + str(counter) + '\n')

out_file.write('\n' + 'octagons' + '\n')
for a in g.nodes_iter():	#calculate number of octagons
	counter = 0
	for b in g.neighbors(a):
		for c in g.neighbors(b):
			for d in g.neighbors(c):
				for e in g.neighbors(d):
					for f in g.neighbors(e):
						for h in g.neighbors(f):
							for i in g.neighbors(h):
								if g.has_edge(a,i):
									#print a,b,c,d,e,f,h,i
									counter = counter + 1
								else:
									pass
	out_file.write(a + '\t' + str(counter) + '\n')
"""
#comm = nx.communicability(g)
#tri = nx.triangles(g)
#trans = nx.transitivity(g)
#clus = nx.clustering(g)
#s_clu = nx.square_clustering(g)
#a_n_c = nx.average_node_connectivity(g)
#e_conn = nx.edge_connectivity(g)
#n_conn = nx.node_connectivity(g)
#f_cli = list(nx.find_cliques(g))
#print f_cli
#a_p_n_c = nx.all_pairs_node_connectivity(g)
#diam = nx.diameter(g)
#radi = nx.radius(g)
#print nx.is_distance_regular(g)
#print nx.is_eulerian(g)
#print nx.is_isolates(g)
#print list(nx.all_shortest_paths(g,source='algae',target='Pisaster'))
#print list(nx.all_simple_paths(g,source='algae',target='Pisaster'))