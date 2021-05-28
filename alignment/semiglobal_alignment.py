#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 16:35:14 2021

@author: chengzev
"""

from Bio import SeqIO, pairwise2

with open('rosalind_smgb.txt', 'r', encoding='utf-8') as handle:
    s, t = map(lambda x: getattr(x, 'seq'), SeqIO.parse(handle, 'fasta'))
    
align = pairwise2.align.globalms(s, t, 1, -1, -1, -1, penalize_end_gaps=(False, False), one_alignment_only=True)[0]
print('{}\n{}\n{}'.format(int(align.score), align.seqA, align.seqB))
