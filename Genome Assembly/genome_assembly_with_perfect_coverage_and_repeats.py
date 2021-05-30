#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 19:37:40 2021

@author: zev
"""

class DeBruijnGraph:
    class Node:
        def __init__(self, kmer):
            self.kmer, self.neighborhood = kmer, []
            self.inDegree, self.outDegree = 0, 0

    def __init__(self, reads):
        self.g, self.head = dict(), reads[0]
        for read in reads:
            hash_left, hash_right = hash(read[:-1]), hash(read[1:])
            node_left, node_right = self.g.setdefault(hash_left, self.Node(read[:-1])), \
                self.g.setdefault(hash_right, self.Node(read[1:]))
            node_left.neighborhood.append(node_right)
            node_left.outDegree += 1
            node_right.inDegree += 1

    def circular_string(self):
        heads = [node for node in self.g.values() if node.outDegree > 1]
        assert len(heads)
        contigs, strings, k = [], set(), len(self.head)-1

        def dfs(node, string):
            if node.outDegree > 1: contigs.append(string+node.kmer[-1])
            else: dfs(node.neighborhood[0], string+node.kmer[-1])

        for head in heads:
            for neighbor in head.neighborhood: dfs(neighbor, head.kmer)

        def recursion(current, string, pos):
            if len(pos) == len(contigs): strings.add(string)
            else:
                for p in set(range(len(contigs))).difference(pos):
                    if not current.endswith(contigs[p][:k]): continue
                    recursion(contigs[p], string+contigs[p][:-k], pos+(p,))

        for i, contig in enumerate(contigs):
            if contig.startswith(self.head):
                recursion(contig, contig[:-k], (i,))
                break
        return strings

if __name__ == '__main__':
    with open('rosalind_grep.txt', encoding='utf-8') as handle:
        reads_ = handle.read().split()
    g_ = DeBruijnGraph(reads_)
    with open('output_grep.txt', 'w', encoding='utf-8') as f:
        print(*g_.circular_string(), sep='\n', file=f)
