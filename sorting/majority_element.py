#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 22:08:24 2021

@author: zev
"""

def find_majority_element(array, n):
    "摩尔投票算法"
    maj, count = None, 0
    for a in array:
        if count == 0:
            maj = a
            count = 1
        elif maj == a: count += 1
        else: count -= 1
    if count == 0: return -1
    elif count > 1: return maj
    elif array.count(maj) >= n // 2 + 1: return maj
    else: return -1        

if __name__ == '__main__':
    with open('rosalind_maj.txt', encoding='utf-8') as handle:
        k_, n_ = map(int, handle.readline().split())
        print(*(find_majority_element(i.split(), n_) for i in handle), sep=' ')    
