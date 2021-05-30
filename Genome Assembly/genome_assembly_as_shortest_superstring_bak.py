#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 20 17:11:52 2021

@author: zev
"""

from itertools import permutations
from Bio import SeqIO

class DirectedAcyclicGraph:
    class Vertex:
        def __init__(self, key, read):
            self.id, self.read = key, read
            self.connectedTo = dict()
            self.distance = float('inf')
            self.predecessor = None
            self.inDegree = 0

    def __init__(self, reads):
        self.n, self.time = len(reads), 1
        self.g = {i: self.Vertex(i, read) for i, read in enumerate(reads)}
        source = self.Vertex(-1, '') # 虚拟顶点连通不同单元

        for i, j in permutations(range(self.n), 2):
            left, right = self.g[i], self.g[j]
            limit = max(len(left.read), len(right.read)) // 2
            begin = 0
            while True:
                pos = left.read.find(right.read[:limit], begin)
                if pos == -1: break
                if right.read.startswith(left.read[pos:]):
                    left.connectedTo[j] = len(right.read) - len(left.read) + pos
                    right.inDegree += 1
                    break
                else:
                    begin = pos + 1

        source.connectedTo = {i: 0 for i, v in self.g.items() if not v.inDegree}
        self.g[-1] = source
        self.bellman_Ford_algorithm()

    def get_shortest_superstring(self):
        string = ''
        current = self.g[-1]
        while True:
            for i, weight in current.connectedTo.items():
                neighbor = self.g[i]
                if neighbor.predecessor is current:
                    string += neighbor.read[-weight:]
                    current = neighbor
                    break
            else:
                break
        return string

    def bellman_Ford_algorithm(self): # 最短路径
        self.g[-1].distance = 0
        update, i = True, 0
        while i < self.n and update:
            update = False
            for v in self.g.values():
                for j, weight in v.connectedTo.items():
                    neighbor = self.g[j]
                    if (v.distance + weight) < neighbor.distance:
                        neighbor.distance = v.distance + weight
                        neighbor.predecessor = v
                        update = True
            i += 1

if __name__ == '__main__':
    with open('rosalind_long.txt') as handle:
        record = SeqIO.parse(handle, 'fasta')
        reads_ = [str(r.seq) for r in record]
    dag = DirectedAcyclicGraph(reads_)
    print(dag.get_shortest_superstring())
