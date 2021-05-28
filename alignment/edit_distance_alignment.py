#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 00:09:30 2021

@author: chengzev
"""

from Bio import SeqIO, pairwise2

with open('rosalind_edta.txt', 'r', encoding='utf-8') as handle:
    s, t = map(lambda r: getattr(r, 'seq'), SeqIO.parse(handle, 'fasta'))
    
align = pairwise2.align.globalms(s, t, 0, -1, -1, -1, one_alignment_only=True)
print(-int(align[0].score), align[0].seqA, align[0].seqB, sep='\n')
