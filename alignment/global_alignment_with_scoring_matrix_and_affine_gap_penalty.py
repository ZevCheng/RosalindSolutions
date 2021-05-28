#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 00:57:10 2021

@author: chengzev
"""

from Bio import pairwise2, SeqIO
from Bio.SubsMat import MatrixInfo

with open('rosalind_gaff.txt', 'r', encoding='utf-8') as handle:
    s, t = map(lambda x: getattr(x, 'seq'), SeqIO.parse(handle, 'fasta'))

align = pairwise2.align.globalds(s, t, MatrixInfo.blosum62, -11, -1, one_alignment_only=True)[0]
print('{}\n{}\n{}'.format(int(align.score), align.seqA, align.seqB))
