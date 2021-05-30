#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 22:03:30 2021

@author: zev
"""

def insertion_sort(array, n):
    count = 0
    for i in range(1, n):
        k = i
        while k > 0 and array[k] < array[k-1]:
            array[k], array[k-1] = array[k-1], array[k]
            k -= 1
            count += 1
    return count

if __name__ == '__main__':
    with open('rosalind_ins.txt', encoding='utf-8') as handle:
        n_ = eval(handle.readline())
        array_ = list(map(int, handle.readline().split()))
    print(insertion_sort(array_, n_))
