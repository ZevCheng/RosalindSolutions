#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 17:29:14 2021

@author: chengzev
"""

from Bio import SeqIO, pairwise2

with open('rosalind_mgap.txt', 'r', encoding='utf-8') as handle:
    s, t = map(lambda x: getattr(x, 'seq'), SeqIO.parse(handle, 'fasta'))

align = pairwise2.align.globalms(s, t, 1, -5000, -1, -1, one_alignment_only=True)[0]
print(align.seqA.count('-') + align.seqB.count('-'))
