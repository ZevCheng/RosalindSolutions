#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 15 17:29:26 2021

@author: zev
"""

from sais import SuffixArray, lcp_array
from suffix_tree import SuffixTree

if __name__ == '__main__':
    with open('rosalind_suff.txt', encoding='utf-8') as handle:
        string = handle.readline().rstrip()
    suffixarray = SuffixArray(string)
    sa, array = suffixarray.sa, suffixarray.array
    lcp = lcp_array(array, sa)
    suffixtree = SuffixTree(sa, lcp)
    stack = [suffixtree.root]
    while stack:
        node = stack.pop()
        for c in node.child:
            print(string[c.location: c.location+c.length])
            if len(c.child): stack.append(c)
