#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 14 22:22:14 2021

@author: zev
"""

from collections import namedtuple

class SuffixTree:
    vertex = namedtuple('Vertex', 'pos lcp left right') # suffix index; lcp; left child; right child 
    border = namedtuple('Border', 'begin end child') # left child and right child

    class Node:
        def __init__(self, label=None, location=None, length=None):
            self.label = label          # suffix index
            self.location = location    # string[location: location+length] -> edge
            self.length = length
            self.child = []             # child is empty -> leaf node

        def isleaf(self): return len(self.child) == 0

    def __init__(self, sa, lcp_array):  # SuffixArray; LcpArray
        self.sa, self.m = sa, len(sa)
        root = self.vertex(lcp_array[1:].argmin()+1, lcp_array[1:].min(), [], [])
        stack = [self.border(1, root.pos, root.left), 
                 self.border(root.pos+1, self.m, root.right)]
        while stack:
            b = stack.pop()
            if b.begin < b.end:
                c = self.vertex(lcp_array[b.begin: b.end].argmin()+b.begin, lcp_array[b.begin: b.end].min(), [], [])
                b.child.append(c)
                stack.extend((self.border(b.begin, c.pos, c.left), self.border(c.pos+1, b.end, c.right)))

        self.root = self.Node(length=root.lcp)
        self._dfs(root, self.root)

    def _dfs(self, vertex, node):
        if not len(vertex.right):
            n = self.Node(label=self.sa[vertex.pos])
            n.location = n.label + vertex.lcp
            n.length = self.m - n.label - vertex.lcp
            node.child.append(n)
        else:
            if vertex.lcp != vertex.right[0].lcp:
                n = self.Node(label=self.sa[vertex.right[0].pos], length=vertex.right[0].lcp-vertex.lcp)
                n.location = n.label + vertex.lcp
                node.child.append(n)
                self._dfs(vertex.right[0], n)
            else:
                self._dfs(vertex.right[0], node)

        if not len(vertex.left):
            n = self.Node(label=self.sa[vertex.pos-1])
            n.location = n.label + vertex.lcp
            n.length = self.m - n.label - vertex.lcp
            node.child.append(n)
        else:
            if vertex.lcp != vertex.left[0].lcp:
                n = self.Node(label=self.sa[vertex.left[0].pos-1], length=vertex.left[0].lcp-vertex.lcp)
                n.location = n.label + vertex.lcp
                node.child.append(n)
                self._dfs(vertex.left[0], n)
            else:
                self._dfs(vertex.left[0], node)
