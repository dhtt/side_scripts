#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 21 12:16:45 2018

@author: dohoangthutrang
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 19 01:09:28 2018

@author: dohoangthutrang
"""
import sys
import os
# taking the arguments and reading the files
arguments = sys.argv
if len(arguments) == 1:
    print("Please give Pattern first then the path to text file")
    sys.exit()

pattern_file = arguments[1]
text_file = arguments[2]

if not os.path.exists(text_file):
    print("Path given to text does not exist, try again")
    sys.exit()

g = open(pattern_file, "r")
pattern = g.read().splitlines()[0]
g.close()

f = open(text_file, "r")
T = f.read().splitlines()[0]
if not T.endswith("$"):
    T = T + "$"
f.close()

def suffix(T):
    def suf(i):
        return T[i:]
    return suf


def suffxarray(T):
    pos = list(range(len(T)))
    pos.sort(key = suffix(T))
    return pos


def longest_prefix(p, s):
    x = 0
    if len(p)<len(s):
        for i in range(len(p)):
            if p[i] == s[i]:
                x += 1
            else:
                break
        return x
    else:
        for i in range(len(s)):
            if p[i] == s[i]:
                x += 1
            else:
                break
        return x
        

sa = suffxarray(T)
suf_strings = []
positions = []
for i,p in enumerate(sa):
    suf_strings.append(T[p:])
    positions.append(i)
interval = [] 


"""Find the Lower boundary"""  
L = 0
R = positions[-1]
R - L >= 0 

while True:
    M = (R+L)//2 
    l = longest_prefix(suf_strings[L],pattern)
    r = longest_prefix(suf_strings[R],pattern)
    mlr = min(l,r)    
    if suf_strings[M][mlr:len(pattern)] >= pattern[mlr:]:
        R = M
    else:
        L = M
    if R - L == 1:
        Lp = R
        interval.append(Lp)
        break


"""Find the Upper boundary"""
L = 0
R = positions[-1]
R - L >= 0 
while True:
    M = (R+L)//2 
    l = longest_prefix(suf_strings[L],pattern)
    r = longest_prefix(suf_strings[R],pattern)
    mlr = min(l,r)    
    if suf_strings[M][mlr:len(pattern)] > pattern[mlr:]:
        R = M
    else:
        L = M
    if R - L == 1:
        if R == positions[-1]:
            Rp = R
        else:
            Rp = L
        interval.append(Rp)
        break

"""Return result"""
print ("Text = ",T)
print ("Pattern = ",pattern)
if Rp - Lp < 0:
    print ("NotFound")
else:
    print ("Interval= ",interval)

