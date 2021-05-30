#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 21:39:27 2021

@author: zev
"""

from Bio import Seq
from itertools import chain

class DeBruijnGraph:
    class Node:
        def __init__(self, kmer):
            self.kmer = kmer
            self.neighborhood, self.inDegree, self.outDegree = set(), 0, 0
        
        def isbalance(self):
            return self.inDegree == self.outDegree
            
        def __hash__(self):
            return hash(self.kmer)

    def __init__(self, reads):
        self.reads, self.reads_rc = reads, [Seq.reverse_complement(read) for read in reads]
        self.g, self.length = dict(), len(reads[0])
        self.k = (self.length-2) if self.length % 2 else (self.length-1)

    def _isEulerian(self):
        return all(node.isbalance() for node in self.g.values())

    def _build_graph(self):
        while True:
            for read in chain(self.reads, self.reads_rc):
                for i in range(self.length - self.k):
                    hash_left, hash_right = hash(read[i:self.k+i]), hash(read[i+1:self.k+1+i])
                    node_left, node_right = self.g.setdefault(hash_left, self.Node(read[i:self.k+i])), \
                        self.g.setdefault(hash_right, self.Node(read[i+1:self.k+1+i]))
                    if node_right not in node_left.neighborhood:
                        node_left.neighborhood.add(node_right)
                        node_left.outDegree += 1
                        node_right.inDegree += 1
            if self._isEulerian(): break
            else: 
                self.k -= 2       
                self.g.clear()

    def cyclic_superstring(self):
        self._build_graph()
        g = self.g.copy()
        _, head = g.popitem()
        string, node = head.kmer, head
        while True:
            neighbor = node.neighborhood.pop()
            string += neighbor.kmer[-1]
            if neighbor not in g.values(): break
            else: node = neighbor
        return string[:-self.k]

if __name__ == '__main__':
    with open('rosalind_gasm.txt', encoding='utf-8') as handle:
        reads_ = handle.read().split()
    g_ = DeBruijnGraph(reads_)
    print(g_.cyclic_superstring())
