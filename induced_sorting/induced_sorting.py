#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 31 19:29:24 2018

@author: dohoangthutrang
"""
from collections import Counter
s ="ABACABA$"
n = len(s)

# =============================================================================
# Type array
# =============================================================================
def getType(s):
    Type = [0]*n
    for i in range(n-2,-1,-1):
        if s[i] > s[i+1]:
            Type[i] = 1
        elif s[i] == s[i+1]:
            Type[i] = Type[i+1]
    return (Type)
    
def getp(s):
    Type = getType(s)  
    p = []
#    Get LMS position:
    if Type[0] == 0:
            p.append(0)
    for i in range(1,n):
        if Type[i] == 0 and Type[i-1] == 1:
            p.append(i)
    return (p)
    
def LMS(s):            
#    Get LMS substrings:
    LMS_substring = []
    p = getp(s)
    for i in range(len(p)-1):
        lms = s[p[i]:p[i+1]+1]
        LMS_substring.append(lms)
    LMS_substring.append(s[n-1])
    print (LMS_substring)
    return (LMS_substring)   
'''any var inside the function should be define first'''
def getBucket(s):
    char = []
    
    count = 1
    for character in s:
        if character not in char:
            char.append(character)
    char = sorted(char)
    print ('char',char)
    
    count = 0
    bucket_tail = []
    for c in char:
        for character in s:
            if c == character: 
                count +=1
        bucket_tail.append(count)
    print ('bucket_tail',bucket_tail)
    
    p = getp(s)
    size = [1, 4, 2, 1]
    bucket = []
    
    for c in char:
        pos = []
        for i in p:
            if s[i] == c:            
                pos.append(i)
        n = len(pos) 
        if n < size[char.index(c)]:
            pos.insert(0,-1)
        bucket += pos
    print ('Bucket',bucket)

LMS(s)
getBucket(s)
#SAIS(s)