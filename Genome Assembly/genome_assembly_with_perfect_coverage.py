#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 17:18:34 2021

@author: zev
"""

class DeBruijnGraph:
    class Node:
        def __init__(self, kmer):
            self.kmer = kmer
            self.neighborhood = set()
        
        def __hash__(self):
            return hash(self.kmer)

    def __init__(self, reads):
        self.reads = reads
        self.g = dict()
    
    def _build_graph(self):
        for read in self.reads:
            hash_left, hash_right = hash(read[:-1]), hash(read[1:])
            node_left, node_right = self.g.setdefault(hash_left, self.Node(read[:-1])), \
                self.g.setdefault(hash_right, self.Node(read[1:]))
            node_left.neighborhood.add(node_right)
    
    def cyclic_superstring(self):
        self._build_graph()
        g = self.g.copy()
        _, head = g.popitem()
        string = head.kmer
        k = len(string)

        def recursion(node):
            neighbor = node.neighborhood.pop()
            if neighbor in g.values():
                return node.kmer[-1] + recursion(neighbor)
            else: return ''

        return string[:-1] + recursion(head)[:2-k]

if __name__ == '__main__':
    with open('rosalind_pcov.txt', encoding='utf-8') as handle:
        reads_ = handle.read().split()
    g_ = DeBruijnGraph(reads_)
    print(g_.cyclic_superstring())
