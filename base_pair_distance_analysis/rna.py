#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from Bio import SeqIO
from collections import Counter
from collections import defaultdict 
from itertools import combinations 
import matplotlib.pyplot as plt
import random
random.seed(1)
import numpy as np

def find_occurence(char, string):
    return [i for i, letter in enumerate(string) if letter == char]

    
record_id = []
record_seq = []

#Open FASTA file and get all sequence, including pairing mask:
with open("tmrna_ali.fasta.txt", "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        record_id.append(record.id) #Store sequence id
        record_seq.append(record) #Store sequence
        
#Process sequence, only store where there is base pairing:        
processed_seq = []
for seq in record_seq[1:]:
    s, s1, s2 = '', '' , ''
    for i in range(len(record_seq[0])):
        if record_seq[0][i] != '-' and seq[i] != '-' and seq[i] != 'n' and seq[i].isupper():
            s+= record_seq[0][i] 
            s1+= seq[i]
    processed_seq.append([s,s1])
    
# Create a dict(all_basepair) to store all base pairs in each sequence for the whole dataset 
pos = defaultdict(list)
all_basepair = defaultdict(dict)
c = 1 #Just counter to keep track of sequence ID
for item in processed_seq:
    template = item[0]
    seq = item[1]
    counts = Counter(template)
    for key in counts.keys():
        seq_pair = []
        len_seq = counts[key]//2
        pos[key] = find_occurence(key, template)
        seq_pair.append((seq[pos[key][0] : pos[key][0]+ len_seq], seq[pos[key][-1] - len_seq +1 : (pos[key][-1]+1)]))
#        print(template[pos[key][0] : pos[key][0]+ len_seq], template[pos[key][-1] - len_seq +1 : (pos[key][-1]+1)])
        all_basepair[record_id[c]][key] = seq_pair
    c+=1
    
####################################
#Extract_basepair:
def extract_basepair(sequence_id):
    return(all_basepair[sequence_id])
    
    
test_seq1 = record_id[1]
test_seq2 = record_id[2]
test_seq3 = record_id[50]
print('All base pair forming loops in sequence: ', extract_basepair(test_seq1))


#Generate a library of random pairs
all_pairs = []
for seq1, seq2 in combinations(record_seq[1:], 2): 
    all_pairs.append([seq1,seq2])
     
#Hamming distance and histogram
def Hamming_distance(seq1, seq2):
    err = 0
    for i in range(len(seq1)):
        if seq1[i] != '-' and seq2[i] != '-' and seq1[i] != 'n' and seq2[i] != 'n':
            if seq1[i] != seq2[i]:
                err += 1
    return(err)
    
total_Hamming_distance = []
for i in (np.random.choice(len(all_pairs), size = 500, replace = False)):
    seq1 = all_pairs[i][0]
    seq2 = all_pairs[i][1]
    total_Hamming_distance.append(Hamming_distance(seq1, seq2))
figure1 = plt.hist(total_Hamming_distance)
plt.xlabel('Hamming distance')
plt.ylabel('Frequency')
title = 'Histogram of Hamming distance'
plt.savefig("Hamming distance of 500 pairs")
plt.close()

#total_Hamming_distance_all_pairs =[]
#for i in all_pairs:
#    seq1 = i[0]
#    seq2 = i[1]
#    total_Hamming_distance_all_pairs.append(Hamming_distance(seq1,seq2))
#figure2 =plt.hist(total_Hamming_distance_all_pairs)
#plt.xlabel('Hamming distance')
#plt.ylabel('Frequency')
#title = 'Histogram of Hamming distance of all pairs'
#plt.savefig('Hamming distance of all pairs')
#plt.close()

#Basepair distance and histogram
def basepair_distance(seq1, seq2):
    err = 0
    for k in all_basepair[seq1].keys():
        if k in all_basepair[seq2].keys():
            if all_basepair[seq1][k] != all_basepair[seq2][k]:
                err += 1
        else:
            err += 1
    for k in all_basepair[seq2].keys():
        if k not in all_basepair[seq1].keys():
            err += 1
    return(err)
    
print('#'*30)
print('Test base pair distance 1: ', basepair_distance(test_seq1, test_seq2))
print('#'*30)
print('Test base pair distance 2: ', basepair_distance(test_seq2, test_seq3))

total_basepair_distance = []
existing_pair = []
while len(total_basepair_distance) < 500:
    a = random.randint(1, len(record_id)-1)
    b = random.randint(1, len(record_id)-1)
    if a != b and (a,b) not in existing_pair:
        seq1 = record_id[a]
        seq2 = record_id[b]
        total_basepair_distance.append(basepair_distance(seq1,seq2))
        existing_pair.append((a,b))
        existing_pair.append((b,a))
figure3 = plt.hist(total_basepair_distance)
plt.xlabel('Base pair distance')
plt.ylabel('Frequency')
title = 'Histogram of base pair distance'
plt.savefig("Base pair distance of 500 pairs")
plt.close()

