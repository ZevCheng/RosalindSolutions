#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 18:13:35 2021

@author: chengzev
"""

from Bio import SeqIO
from functools import lru_cache

@lru_cache(maxsize=None)
def edit_distance(s, t):
    if not s or not t:
        return len(s) + len(t)
    elif s[-1] == t[-1]:
        return edit_distance(s[:-1], t[:-1])
    else:
        return 1 + min(edit_distance(s[:-1], t[:-1]), edit_distance(s[:-1], t), \
                       edit_distance(s, t[:-1]))

if __name__ == '__main__':
    with open('rosalind_edit.txt', 'r', encoding='utf-8') as handle:
        s_, t_ = map(lambda r: getattr(r, 'seq'), SeqIO.parse(handle, 'fasta'))
    print(edit_distance(s_, t_))
