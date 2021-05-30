#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 15:56:05 2021

@author: zev
"""

from Bio import Seq
from itertools import chain

class DeBruijnGraph:
    class Node:
        def __init__(self, kmer):
            self.kmer = kmer
            self.neighborhood = set()
        
        def __hash__(self):
            return hash(self.kmer)

    def __init__(self, reads):
        self.reads = reads
        self.reads_rc = [Seq.reverse_complement(read) for read in reads]
        self.g = dict()
        
    def _build_graph(self):
        for read in chain(self.reads, self.reads_rc):
            key_left, key_right = hash(read[:-1]), hash(read[1:])
            node_left, node_right = self.g.setdefault(key_left, self.Node(read[:-1])), \
                self.g.setdefault(key_right, self.Node(read[1:]))
            node_left.neighborhood.add(node_right)
    
    def adjacency_list(self):
        self._build_graph()
        f = open('output_dbru.txt', 'w', encoding='utf-8')
        for node in self.g.values():
            for neighbor in node.neighborhood:
                print('({}, {})'.format(node.kmer, neighbor.kmer), file=f)
        f.close()

if __name__ == '__main__':
    with open('rosalind_dbru.txt', encoding='utf-8') as handle:
        reads_ = [h.rstrip() for h in handle]
    g_ = DeBruijnGraph(reads_)
    g_.adjacency_list()
