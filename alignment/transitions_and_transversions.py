#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 16:17:20 2021

@author: chengzev
"""

from Bio import SeqIO

transitions_map = {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C'}
with open('rosalind_tran.txt', 'r', encoding='utf-8') as handle:
    s1, s2 = map(lambda r: getattr(r, 'seq'), SeqIO.parse(handle, 'fasta'))
    ls = [1 if transitions_map[i] == j else 0 for i, j in zip(s1, s2) if i != j]
print(ls.count(1) / ls.count(0))
