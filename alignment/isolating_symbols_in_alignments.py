#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 19:25:55 2021

@author: zev
"""

from Bio import SeqIO
import numpy as np
import numba

@numba.jit(nopython=True)
def global_alignment(s, t):
    n, m = len(s), len(t)
    match = np.empty((n+1, m+1), dtype=np.int32)
    insert_x = np.empty_like(match)
    insert_y = np.empty_like(match)
    match[0, :], match[:, 0] = -m-n-1, -m-n-1
    insert_x[0, :], insert_y[:, 0] = -m-n-1, -m-n-1
    match[0, 0] = 0
    insert_x[:, 0] = np.arange(0, -n-1, -1)
    insert_y[0, :] = np.arange(0, -m-1, -1)
    
    for i, a in enumerate(s, 1):
        for j, b in enumerate(t, 1):
            x = 1 if a == b else -1
            match[i, j] = max(match[i-1, j-1]+x, insert_x[i-1, j-1]+x, insert_y[i-1, j-1]+x)
            insert_x[i, j] = max(match[i-1, j]-1, insert_x[i-1, j]-1)
            insert_y[i, j] = max(match[i, j-1]-1, insert_y[i, j-1]-1)
    score = max(match[-1, -1], insert_x[-1, -1], insert_y[-1, -1])
    return match, score

@numba.jit(nopython=True)
def get_matrix_m(s, t):
    n, m, count = len(s), len(t), 0
    match, score = global_alignment(s, t)
    match_r, _ = global_alignment(s[::-1], t[::-1])

    for i, a in enumerate(s):
        for j, b in enumerate(t):
            x = -1 if a == b else 1
            count += match[i+1, j+1] + match_r[n-i, m-j] + x
    return score, count


if __name__ == '__main__':
    with open('rosalind_osym.txt', 'r', encoding='utf-8') as handle:
        s_, t_ = map(lambda x: getattr(x, 'seq'), SeqIO.parse(handle, 'fasta'))
    result = get_matrix_m(str(t_), str(s_))
    print(result)
