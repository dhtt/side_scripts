#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  6 01:05:20 2019
"""

import sys, time
import itertools

def totalDistance(word,DNA):
    totalDistance = 0 #Initiate a value
    best_s = [] #Initiate a list
    for seq in DNA:
        distance = [0]*(len(seq)-len(word)+1) #Initiate a list of positions that word and seq can align
        for i in range(len(word)):
            for p in range(len(seq)-len(word)+1):
                if word[i] != seq[p+i]:
                    distance[p] += 1 #Compare, everytime they are different, distance at p position increase by 1
        best_s.append(distance.index(min(distance)))
        totalDistance += min(distance) #Add the min distance for each seq into total distance
    return totalDistance, best_s


def BruteForceMedianSearch(DNA, l):
    min_distance = l*len(DNA[0]) #Initiate a value
    best_motif = 'A'*l #Initiate a string
    for combination in itertools.product('ATGC', repeat=l): #Generate combination of l character from A G C T
        motif = ''.join(combination) #Because itertools.product yields a list, we have to concatinate the element inside to make a string
        if totalDistance(motif,DNA)[0] <= min_distance: #everytime calculate distance of motif and DNA, compare to the current smallest distance
            #if the new distance is smaller then update min_distance, best_motif and best_s
            min_distance, best_s = totalDistance(motif,DNA) 
            best_motif = motif
            print(best_motif, min_distance, best_s)
    return(best_motif, min_distance, best_s)
    
#def nextWord(word):

DNA = []
fp = open(sys.argv[1])
for line in fp:
	line = line.strip()
	DNA.append(line)

startTime = time.time()
BruteForceMedianSearch(DNA, 8)
elapsed = time.time() - startTime
print('It spends', elapsed, 'seconds')
