#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  2 21:31:11 2018

@author: dohoangthutrang
"""
import numpy
from itertools import product

number_of_genes = 6
number_of_possible_states = 2**6
'''A: Find all possible initial states:'''
possible_states = []
for i in product([0,1], repeat=number_of_genes):
    possible_states.append(''.join(map(str, i)))

'''B: Regulatory network simulation'''
orbit = []    
for starting_state in possible_states:
    states = [starting_state]
    steps = 20
    nodes = ['A','B','C','D','E','F']
    
    #Give the next 20 states for each initial state found in A
    while len(states)<steps:
        s = states[-1]
        if s[0] == '1':
            A = '1'
        else: 
            A = '0'
        
        if (s[0] == '1' or s[2] == '1') and s[3] =='0':
            B = '1'
        else:
            B= '0'
        
        if s[1] == '1':
            C = '1'
        else: 
            C = '0'
            
        if s[5] == '1':
            D = '1'
        else: 
            D = '0'
            
        if s[2] == '1' and s[3] == '0':
            E = '1'
        else:
            E = '0'
            
        if s[4] =='1' and s[3] == '0':
            F = '1'
        else:
            F = '0'
            
        ss = A+B+C+D+E+F
        states.append(ss)
    #Convert  binary levels of the genes into integer and report 10 next states for each initial state:
    initial_states_set = set()
    initial_states_list = []
    for state in states:
        initial_state = 0
        for index,char in enumerate(state):
            if char == '1': initial_state += 2**index
        initial_states_list.append(initial_state)
    print ('Initial state %s: ' %(starting_state), initial_states_list)
    orbit.append(initial_states_list)
        
'''C: Summary'''
orbit =numpy.asarray(orbit)
attractor = [0,1,4,5]
#Find basin of attractor
attraction_basin_list = []
for att in attractor: 
    attraction_basin = set()
    for row in orbit:
        if att in row:
            for i in row: attraction_basin.add(i)
    attraction_basin_list.append(attraction_basin)
#Calculate coverage
coverage = 0
for item in attraction_basin_list:
    coverage = len(item)/(number_of_possible_states)*100
    print(item, coverage,'%')