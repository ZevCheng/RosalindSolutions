#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 00:34:49 2021

@author: chengzev
"""

import numpy as np
from Bio import SeqIO
import numba

@numba.jit(nopython=True)
def edit_distance_alignment(s, t):
    n, m = len(s), len(t)
    score = np.zeros((n + 1, m + 1), dtype=np.uint32)
    path = np.zeros((n + 1, m + 1), dtype=np.uint32)
    score[0, :] = [i for i in range(m + 1)]
    score[:, 0] = [j for j in range(n + 1)]
    path[0, :] = np.full(m + 1, 1)
    path[:, 0] = np.full(n + 1, 1)

    for i, a in enumerate(s, 1):
        for j, b in enumerate(t, 1):
            x = 0 if a == b else 1
            temp = (score[i-1, j] + 1, score[i-1, j-1] + x, score[i, j-1] + 1)
            score[i, j] = y = min(temp)
            for k, c in enumerate(temp):
                if c != y:
                    continue
                if k == 0: path[i, j] += path[i-1, j] % 134217727
                elif k == 1: path[i, j] += path[i-1, j-1] % 134217727
                else: path[i, j] += path[i, j-1] % 134217727

    return path[n, m] % 134217727

            
if __name__ == '__main__':
    with open('rosalind_ctea.txt', 'r', encoding='utf-8') as handle:
        s_, t_ = map(lambda r: getattr(r, 'seq'), SeqIO.parse(handle, 'fasta'))
    print(edit_distance_alignment(str(s_), str(t_)))
