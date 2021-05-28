#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 18:23:20 2021

@author: chengzev
"""

from Bio import SeqIO
import numpy as np
import numba

table = []
with open('blosum62.txt', 'r', encoding='utf-8') as file:
    syms = np.array([ord(r) for r in file.readline().split()], dtype=np.int8)
    for line in file:
        table.append([eval(i) for i in line.split()[1:]])
table = np.array(table, dtype=np.int8)

@numba.jit(nopython=True)
def affine_penalty_local_alignment(s, t):
    global syms, table
    n, m, limit = len(s), len(t), -2147483648 // 2
    match, path_match = np.zeros((n + 1, m + 1), dtype=np.int32), np.zeros((n + 1, m + 1), dtype=np.int8)
    insert_row, insert_col = np.zeros_like(match), np.zeros_like(match)
    path_row, path_col = np.zeros_like(path_match), np.zeros_like(path_match)
    insert_row[0, :], insert_row[:, 0] = np.full(m + 1, limit), np.full(n + 1, limit)
    insert_col[0, :], insert_col[:, 0] = np.full(m + 1, limit), np.full(n + 1, limit)
    
    for i, a in enumerate(s, 1):
        for j, b in enumerate(t, 1):
            score = table[np.where(syms == ord(a))[0][0], np.where(syms == ord(b))[0][0]]
            d = np.array([0, match[i-1, j-1] + score, insert_row[i-1, j-1] + score, insert_col[i-1, j-1] + score])
            x = np.array([match[i-1, j] - 11, insert_row[i-1, j] - 1])
            y = np.array([match[i, j-1] - 11, insert_col[i, j-1] - 1])
            match[i, j], path_match[i, j] = d.max(), d.argmax()
            insert_row[i, j], path_row[i, j] = x.max(), x.argmax()
            insert_col[i, j], path_col[i, j] = y.max(), y.argmax()
    
    mark, r, c = match.max(), '', ''
    position = np.where(match == mark)
    row, col = position[0][0], position[1][0]
    while True:
        p = path_match[row, col]
        if p == 0: break
        row -= 1
        col -= 1
        r += s[row]
        c += t[col]
        if p == 1:
            continue
        elif p == 2:
            while True:
                p = path_row[row, col]
                row -= 1
                r += s[row]
                if p == 0: break
        else:
            while True:
                p = path_col[row, col]
                col -= 1
                c += t[col]
                if p == 0: break

    return mark, r[::-1], c[::-1]

if __name__ == '__main__':
    with open('rosalind_laff.txt', 'r', encoding='utf-8') as handle:
        s_, t_ = map(lambda x: getattr(x, 'seq'), SeqIO.parse(handle, 'fasta'))
    print('{}\n{}\n{}'.format(*affine_penalty_local_alignment(str(s_), str(t_))))
