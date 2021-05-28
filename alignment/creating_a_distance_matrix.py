#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 16:58:38 2021

@author: chengzev
"""

from Bio import SeqIO
import numpy as np

strings = []
with open('rosalind_pdst.txt', 'r', encoding='utf-8') as handle:
    for r in SeqIO.parse(handle, 'fasta'): strings.append(r.seq)

num = len(strings)
distance = np.zeros((num, num), dtype=np.float32)
for i in range(num):
    for j in range(i, num):
        distance[i, j] = distance[j, i] = \
            round(sum(1 for s, t in zip(strings[i], strings[j]) if s != t) / len(strings[i]), 5)
for d in distance: print(*d)
