#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 16:08:28 2021

@author: chengzev
"""

from Bio import pairwise2, SeqIO
from Bio.SubsMat import MatrixInfo

with open('rosalind_gcon.txt', 'r', encoding='utf-8') as handle:
    s, t = map(lambda r: getattr(r, 'seq'), SeqIO.parse(handle, 'fasta'))

score = pairwise2.align.globalds(s, t, MatrixInfo.blosum62, -5, 0, score_only=True)
print(score)
