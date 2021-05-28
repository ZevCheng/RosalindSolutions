#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 14:45:35 2021

@author: chengzev
"""

from Bio import SeqIO, pairwise2
from Bio.SubsMat import MatrixInfo

with open('rosalind_glob.txt', 'r', encoding='utf-8') as handle:
    s_, t_ = map(lambda r: getattr(r, 'seq'), SeqIO.parse(handle, 'fasta'))

score = pairwise2.align.globalds(s_, t_, MatrixInfo.blosum62, -5, -5, score_only=True)
print(score)
