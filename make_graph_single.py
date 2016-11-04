#!/usr/bin/env python
# coding: utf-8

import networkx as nx
import matplotlib.pyplot as plt
import json
import heapq
import math


g = nx.Graph()	#initialize an empty graph for all interactions

in_file = open('kelp_interactions_hab.csv', 'r')	#read in the data file. This file should contain only the interactions of interest
out_file = open('kelp_hab_results.txt', 'w')
desc = 'kelp forest' #add a notation to remind you what network this is
delim = ',' #the delimiter used in the data file

#This function calculates the 'neighbor of neighbor' parameter that I made up. It adds up 
#the unique number of neighbors and neighbors of neighbors for each node. Input is a graph.
#Output is a dictionary.
def nofn(graph):
	nofn_dict = dict()
	for n in graph:
		#print 'node is ' + n
		dn = []
		non = []
		an = nx.all_neighbors(graph,n)
		for i in an:
			#print 'neighbor is ' + i
			anon = nx.all_neighbors(graph,i)
			if i in dn:
				pass
			else:
				dn.append(i)
			for j in anon:
				#print 'neighbor of neighbor is ' + j
				if j in non:
					pass
				else:
					non.append(j)
		tn = len(dn) + len(non)
		#print ndn
		#print nnon
		#print tn
		nofn_dict[n] = tn
	return nofn_dict

#This function normalizes values against the maximum value by dividing all values by the
#maximum value. Input is a dictionary and output is a dictionary.
def max_normalize(dict):
	values = []
	for k in dict:
		values.append(float(dict[k]))
	denom = max(values)
	for k in dict:
		norm_value = float(dict[k])/denom
		dict[k] = norm_value
	return dict

#This function calculates the sum of the edge betweenness centrality values for every node.
# Input is a dictionary and output is a dictionary
def sum_edge_bet(dict):
	taxa = []
	for k in dict:
		#print k
		name_1, name_2 = k
		if name_1 in taxa:
			pass
		else:
			taxa.append(name_1)
		if name_2 in taxa:
			pass
		else:
			taxa.append(name_2)
	new_dict = {}
	for taxon in taxa:
		total = 0.0
		for k in dict:
			name_1, name_2 = k
			if name_1 == taxon or name_2 == taxon:
				total = total + float(dict[k])
		new_dict[taxon] = total
	return new_dict

#this function tests to see if all nodes are connected. If they are not, the rest of the 
#program won't work
def conn_graph(g):
	m = nx.is_connected(g)
	if m == False:
		print 'Graph not connected'
		return 0
	else:
		return 1

#start program: read in file
source_taxon_ids = []	 
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
print all_names
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
			g.add_edge(name,target_taxon, type=interaction_type)
		else:
			pass
nn = g.number_of_nodes()
ne = g.number_of_edges()
print ne

#determine if all nodes are connected in one graph for all interactions
j = conn_graph(g)
if j == 1:	#if all nodes are connected, proceed with calculations
	all_nodes = g.nodes()
	all_edges = g.edges()
	print len(all_nodes)
	print len(all_edges)
	deg_c = nx.degree_centrality(g)
	bet_c = nx.betweenness_centrality(g)
	clo_c = nx.closeness_centrality(g)
	eig_c = nx.eigenvector_centrality_numpy(g)
	edg_b = nx.edge_betweenness_centrality(g,normalized=True)
	s_edg_b = sum_edge_bet(edg_b)
	kat_c = nx.katz_centrality(g,alpha=0.01,beta=1.0,max_iter=1000,normalized=True)
	neib = nofn(g)
	cent = nx.center(g)
	#print cent
	ecc = nx.eccentricity(g)
	peri = nx.periphery(g)
	#print peri
	dom_s = nx.dominating_set(g)
	#print dom_s
	clo_v = nx.closeness_vitality(g)
else:	#if all nodes are NOT connected, pick the largest graph and proceed with calculations
	max_graph = max(nx.connected_component_subgraphs(g), key=len)
	h = list(nx.connected_component_subgraphs(g))  #This is only necessary to find the taxa that are not connected
	print len(h)									#the first time you look at the graph
	all_nodes = max_graph.nodes()
	all_edges = max_graph.edges()
	print len(all_nodes)
	print len(all_edges)
	deg_c = nx.degree_centrality(max_graph)
	bet_c = nx.betweenness_centrality(max_graph)
	clo_c = nx.closeness_centrality(max_graph)
	eig_c = nx.eigenvector_centrality_numpy(max_graph)
	edg_b = nx.edge_betweenness_centrality(max_graph,normalized=True)
	s_edg_b = sum_edge_bet(edg_b)
	kat_c = nx.katz_centrality(max_graph,alpha=0.01,beta=1.0,max_iter=1000,normalized=True)
	neib = nofn(max_graph)
	cent = nx.center(max_graph)
	#print cent
	ecc = nx.eccentricity(max_graph)
	peri = nx.periphery(max_graph)
	#print peri
	dom_s = nx.dominating_set(max_graph)
	#print dom_s
	clo_v = nx.closeness_vitality(max_graph)

#normalize everything
vars = [deg_c,bet_c,clo_c,eig_c,edg_b,s_edg_b,kat_c,neib,ecc,clo_v]
for var in vars:
	var = max_normalize(var)

out_file.write('taxon' + '\t' + 'system' + '\t' + 'degree_centrality' + '\t' + 'betweenness_centrality' + '\t' + 'closeness_centrality' + '\t' + 'eigenvector_centrality' + '\t' + 'sum_edge_bet' + '\t' + 'katz_centrality' + '\t' + 'eccentricity' + '\t' + 'closeness_vitality' + '\t' + 'nofn' + '\t' + 'center' + '\t' + 'periphery' + '\t' + 'dominating_set' + '\n')
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

	out_file.write(n + '\t' + desc + '\t' + str(deg_c[n]) + '\t' + str(bet_c[n]) + '\t' + str(clo_c[n]) + '\t' + str(eig_c[n]) + '\t' + str(s_edg_b[n]) + '\t' + str(kat_c[n]) + '\t' + str(ecc[n]) + '\t' + str(clo_v[n]) + '\t' + str(neib[n]) + '\t' + cn + '\t' + pn + '\t' + dn + '\n')


nx.draw_spring(g) #this line and the next will make a figure. Comment them out if you do not want a figure
plt.show()