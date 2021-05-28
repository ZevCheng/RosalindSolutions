#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 16:41:02 2021

@author: chengzev
"""

from Bio import SeqIO
from Bio.SubsMat import MatrixInfo
import numpy as np


def lcoal_alignment(s, t):
    n, m = len(s), len(t)
    matrix = MatrixInfo.pam250
    score = np.zeros((n + 1, m + 1), dtype=np.uint32)
    # u: 1, d: 2, l: 3
    path = np.zeros((n + 1, m + 1), dtype=np.uint8)
    
    for i, a in enumerate(s, 1):
        for j, b in enumerate(t, 1):
            score[i, j], path[i, j] = \
                max((0, 0), (score[i-1, j] - 5, 1), \
                    (score[i-1, j-1] + matrix.get((a, b), matrix.get((b, a))), 2), \
                        (score[i, j-1] - 5, 3))
    
    mark = score.max()
    temp = np.where(score == mark)
    row, col = temp[0][0], temp[1][0]
    p, s_, t_ = path[row, col], '', ''
    while p:
        if p == 1:
            row -= 1
            s_ += s[row]
        elif p == 2:
            row -= 1
            col -= 1
            s_ += s[row]
            t_ += t[col]
        else:
            col -= 1
            t_ += t[col]
        p = path[row, col]
    return mark, s_[::-1], t_[::-1]


if __name__ == '__main__':
    with open('rosalind_loca.txt', 'r', encoding='utf-8') as handle:
        s1, t1 = map(lambda r: getattr(r, 'seq'), SeqIO.parse(handle, 'fasta'))
    print(*lcoal_alignment(str(s1), str(t1)), sep='\n')
