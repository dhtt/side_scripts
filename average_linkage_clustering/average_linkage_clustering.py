#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 11:22:46 2018

@author: dohoangthutrang
"""

import time
import math
import pandas as pd
import itertools
from collections import defaultdict


start_time = time.time()
df = pd.read_csv('methylation.csv', delimiter= ";")
celltype = df.columns.values.tolist()[7:]
for item in celltype:
    df[item] = [x.replace(',','.') for x in df[item]]
    df[item] = df[item].replace('.',0.0)
    df[item] = df[item].astype(float)

average_met_level = {}
for item in celltype: 
    average = df[item].sum()/len(df[item])
    average_met_level[item] = average
print(average_met_level)

def remove_listinlist(A,B):
    #Remove list item in list A if item contains any element of B
    for i in B:
        for item in A:
            if i in item:
                A.remove(item)
    return(A)
    
def get_minkey(d): 
    #Return key of min values in dictionary
    min_key = min(d, key = d.get)
    return(min_key)

def delete_dictkey(d,l):
    for k in d.copy():
        for p in l:
            if p in k:
                d.pop(k, None)
    return(d)
                
#PART C1
def Euclidean_distance(a,b):
    dist = [(a - b)**2 for a, b in zip(df[a], df[b])]
    dist = math.sqrt(sum(dist))
    return dist
        
cells = ['HSC','CD4','EPro']

for pair in (itertools.combinations(cells,2)):
    print('Euclidean distance between ',pair,': ', Euclidean_distance(pair[0],pair[1]))


#PART C2
A = [] #A: List of list for keeping track of new cluster
B = celltype #B: List of unclustered cell type

all_pair_distance = defaultdict(dict)
for pair in (itertools.permutations(B,2)):
    all_pair_distance[pair[0]][pair[1]] = Euclidean_distance(pair[0],pair[1])/2
    
average_linkage_0 = {} 
for pair in (itertools.combinations(B,2)):
    average_linkage_0[pair] = Euclidean_distance(pair[0],pair[1])/2
            
def average_linkage(A,B):
    if len(B) == 0:
        return (A)
    
    else: 
        #Case 1: Cluster between pair of proteins in B:        
        average_linkage_1 = {}
        for pair in (itertools.combinations(B,2)):
            average_linkage_1[pair] = average_linkage_0[pair]
            
        #Case 2: Cluster between formed clusters in A and a cell type in B:
        for b in B:
            l = [b]
            sum = 0
            for i in range(0, len(A)): 
                for j in range(0,len(A[i])):
                    sum += all_pair_distance[b][A[i][j]]
                    l.append(A[i][j])
                average_linkage_1[tuple(l)] = sum/(len(A[i]))
       
        
        #Case 3: Cluster between all formed clusters in A:
        average_linkage_2 = {}
        for pair in (itertools.combinations(A,2)):
            sum = 0
            for i in pair[0]:
                for j in pair[1]:
                    if i != j:
                        l = pair[0]+pair[1]
                        k = len(pair[0])*len(pair[1])
                        sum += all_pair_distance[i][j]
            average_linkage_2[tuple(l)] = sum/k
           
        ####Find smallest distance between all clusters:
        min_al1 = get_minkey(average_linkage_1)    
        
        #If a formed cluster in A has smaller distance than case1 and case2, 
        #append newly joined cluster to A and delete the daughter clusters in A by remove_listinlist()
        if len(average_linkage_2) > 0 and min(average_linkage_2.values()) < min(average_linkage_1.values()):
            min_al2 = get_minkey(average_linkage_2) 
            remove_listinlist(A, min_al2)
            A.append(list(min_al2))
            print(min_al2,"|",min(average_linkage_2.values()))
            
        #else, delete cluster pair or protein in B, append those to A and delete the daughter cluster(if case2)
        #by remove_listinlist()
        else:
            for i in min_al1:
                if i in B: B.remove(i)
            remove_listinlist(A, min_al1)    
            A.append(list(min_al1))
            print(min_al1,"|",min(average_linkage_1.values()))
            
#        print(A)
#        print(B,'\n')
        return(average_linkage(A,B))
        
average_linkage(A,B)
print("--- %s seconds ---" % (time.time() - start_time))