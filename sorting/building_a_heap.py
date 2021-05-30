#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 02:57:28 2021

@author: zev
"""

def heapify(array, n):
    for i in range(2, n+1):
        while True:
            j = i // 2
            if not j: break
            if array[j-1] < array[i-1]:
                array[j-1], array[i-1], i = array[i-1], array[j-1], j
            else: break

if __name__ == '__main__':
    with open('rosalind_hea.txt', encoding='utf-8') as handle:
        n_ = eval(handle.readline())
        array_ = list(map(int, handle.readline().split()))
    heapify(array_, n_)
    with open('output_hea.txt', 'w', encoding='utf-8') as f:
        print(*array_, file=f)
