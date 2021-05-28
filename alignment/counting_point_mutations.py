#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 21:47:31 2021

@author: chengzev
"""

def hamming_distance(s, t):
    return sum(1 for i, j in zip(s, t) if i != j)


if __name__ == '__main__':
    with open('rosalind_hamm.txt', 'r', encoding='utf-8') as handle:
        s_, t_ = handle.readlines()
    print(hamming_distance(s_, t_))

