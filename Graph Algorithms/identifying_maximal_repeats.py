#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 16 18:17:42 2021

@author: zev
"""

from sais import SuffixArray, lcp_array
from suffix_tree import SuffixTree

def collect_maximal_repeats(string, limit=20):
    if string[-1] != '$': string = string + '$'
    suffixarray = SuffixArray(string)
    array, sa = suffixarray.array, suffixarray.sa
    lcp = lcp_array(array, sa)
    suffixtree = SuffixTree(sa, lcp)
    maximal_repeats = set()

    def dfs(node, s):
        split_type = False
        alphabet = set()
        for c in node.child:
            if c.isleaf(): alphabet.add(string[c.label-1])
            else:
                result = dfs(c, s+string[c.location: c.location+c.length])
                if result is True: split_type = True
                else: alphabet.add(result)
        if split_type or len(alphabet) > 1:
            if len(s) >= limit: maximal_repeats.add(s)
            return True
        else: return alphabet.pop()
    
    dfs(suffixtree.root, '')
    return maximal_repeats

if __name__ == '__main__':
    with open('rosalind_mrep.txt') as handle:
        string = handle.readline().rstrip()
    print(*collect_maximal_repeats(string), sep='\n')
