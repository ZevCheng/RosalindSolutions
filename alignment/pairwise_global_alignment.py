#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 13:37:11 2021

@author: chengzev
"""

from Bio import Entrez, SeqIO
import numpy as np
import numba

@numba.jit(nopython=True)
def global_alignment(s, t):
    n, m = len(s), len(t)
    bases = ('A', 'T', 'G', 'C')
    score = ((5, -4, -4, -4), (-4, 5, -4, -4), (-4, -4, 5, -4), (-4, -4, -4, 5))
    open_gap, extend_gap, limit = -9, -1, -2147483648 // 2
    match = np.zeros((n + 1, m + 1), dtype=np.int32)
    insert_row = np.zeros_like(match)
    insert_col = np.zeros_like(match)
    
    match[0, 1:], match[1:, 0] = np.full(m, limit), np.full(n, limit)
    insert_row[:, 0] = [open_gap + extend_gap * i for i in range(n + 1)]
    insert_row[0, :] = np.full(m + 1, limit)
    insert_col[0, :] = [open_gap + extend_gap * j for j in range(m + 1)]
    insert_col[:, 0] = np.full(n + 1, limit)

    for i, a in enumerate(s, 1):
        for j, b in enumerate(t, 1):
            match[i, j] = max(match[i-1, j-1] + score[bases.index(a)][bases.index(b)], \
                insert_row[i-1, j-1] + score[bases.index(a)][bases.index(b)], \
                insert_col[i-1, j-1] + score[bases.index(a)][bases.index(b)])
            
            insert_row[i, j] = max(match[i-1, j] + open_gap + extend_gap, \
                insert_row[i-1, j] + extend_gap)
                
            insert_col[i, j] = max(match[i, j-1] + open_gap + extend_gap, \
                insert_col[i, j-1] + extend_gap)

    return max(match[n, m], insert_row[n, m], insert_col[n, m])


if __name__ == '__main__':
    Entrez.email = 'zevcheng@yeah.net'
    
    with open('rosalind_need.txt', 'r', encoding='utf-8') as file:
        ids = file.read().replace(' ', ', ')

    with Entrez.efetch(db='nucleotide', id=[ids], rettype='fasta') as handle:
        record1, record2 = SeqIO.parse(handle, 'fasta')

    print(global_alignment(str(record1.seq), str(record2.seq)))
