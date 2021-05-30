#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 05:05:10 2021

@author: zev
"""

from math import ceil

def merge_sort(array, begin, end):
    if end - begin > 1:
        median = ceil((begin+end)/2)

        merge_sort(array, begin, median)
        merge_sort(array, median, end)

        left, right = array[begin: median], array[median: end]

        i, j, k = 0, 0, begin
        while k < end:
            if i == len(left) or j == len(right):
                array[k:end] = left[i:median] + right[j:end]
                break
            elif left[i] <= right[j]:
                array[k] = left[i]
                k += 1
                i += 1
            else:
                array[k] = right[j]
                k += 1
                j += 1

if __name__ == '__main__':
    with open('rosalind_ms.txt', encoding='utf-8') as handle:
        n_ = eval(handle.readline())
        array_ = list(map(int, handle.readline().split()))
    merge_sort(array_, 0, n_)
    with open('output_ms.txt', 'w', encoding='utf-8') as f:
        print(*array_, file=f)

