#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 23:57:22 2021

@author: zev
"""

def quick_sort(array, begin, end):
    if end - begin > 1:
        e, b = three_way_partition(array, begin, end)
        quick_sort(array, begin, e)
        quick_sort(array, b, end)

def three_way_partition(array, begin, end):
    i, k, j = begin+1, begin+1, end-1
    while i <= j:
        if array[i] > array[begin]:
            array[i], array[j] = array[j], array[i]
            j -= 1
        elif array[i] < array[begin]:
            array[i], array[k] = array[k], array[i]
            i += 1
            k += 1
        else: i += 1
    array[begin], array[k-1] = array[k-1], array[begin]
    return k, i

if __name__ == '__main__':
    with open('rosalind_qs.txt', encoding='utf-8') as handle:
        n_ = eval(handle.readline())
        array_ = list(map(int, handle.readline().split()))
    quick_sort(array_, 0, n_)
    with open('output_qs.txt', 'w', encoding='utf-8') as f:
        print(*array_, file=f)
