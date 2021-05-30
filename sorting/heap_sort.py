#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 23:15:39 2021

@author: zev
"""

def heap_sort(array, n):
    for i in range(2, n+1):
        while True:
            j = i // 2
            if not j: break
            elif array[j-1] < array[i-1]:
                array[i-1], array[j-1], i = array[j-1], array[i-1], j
            else: break
    for i in range(n-1, 0, -1):
        array[i], array[0] = array[0], array[i]
        j = 1
        while True:
            if 2*j+1 <= i:
                if array[2*j-1] >= array[2*j]:
                    if array[2*j-1] > array[j-1]:
                        array[2*j-1], array[j-1], j = array[j-1], array[2*j-1], 2*j
                    else: break
                else:
                    if array[2*j] > array[j-1]:
                        array[2*j], array[j-1], j = array[j-1], array[2*j], 2*j+1
                    else: break
            elif 2*j == i:
                if array[2*j-1] > array[j-1]:
                    array[2*j-1], array[j-1] = array[j-1], array[2*j-1]
                break
            else: break

if __name__ == '__main__':
    with open('rosalind_hs.txt', encoding='utf-8') as handle:
        n_ = eval(handle.readline())
        array_ = list(map(int, handle.readline().split()))
    heap_sort(array_, n_)
    with open('output_hs.txt', 'w', encoding='utf-8') as f:
        print(*array_, file=f)
