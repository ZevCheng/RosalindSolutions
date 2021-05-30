#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 14 11:23:11 2021

@author: zev
"""

import numpy as np
from collections import Counter
from numba import njit
from numba.typed import List
import warnings
warnings.filterwarnings('ignore')

class SuffixArray:

    def __init__(self, string):
        if string[-1] != '$': string = string+'$'
        self.array = np.array([ord(s) for s in string], dtype=np.int8)
        self.sa = self._sa_is(self.array, 127)

    def _sa_is(self, array, size): # size -> alphabet size
        char_type = get_char_type(array)
        lms_array = get_lms_array(char_type) # left most substring
        bucket_begin, bucket_end = get_bucket_border(array, size)
        pointer_array = get_pointer_array(array, lms_array, bucket_end.copy())
        induceSortL(array, pointer_array, bucket_begin.copy(), char_type)
        induceSortS(array, pointer_array, bucket_end.copy(), char_type)
        s1, substring_pos, size1 = get_s1(array, pointer_array, lms_array)
        if len(s1) == size1:
            sa1 = np.empty(shape=size1, dtype=np.int32)
            for i, s in enumerate(s1): sa1[s] = i
        else:
            sa1 = self._sa_is(s1, size1)
        sa = get_sa_from_sa1(array, sa1, substring_pos, bucket_end.copy())
        induceSortL(array, sa, bucket_begin, char_type)
        induceSortS(array, sa, bucket_end, char_type)
        return sa

@njit
def get_char_type(array):
    # True -> S; False -> L
    m = len(array)
    char_type = np.empty(shape=m, dtype=np.bool_)
    char_type[-1] = True
    for i, j in zip(range(m-2, -1, -1), range(m-1, 0, -1)):
        if array[i] < array[j]: char_type[i] = True
        elif array[i] > array[j]: char_type[i] = False
        else: char_type[i] = char_type[j]
    return char_type

def get_bucket_border(array, size):
    counter = Counter(array).most_common()
    counter.sort()
    bucket_begin = np.zeros(shape=size, dtype=np.int32)
    bucket_end = np.zeros(shape=size, dtype=np.int32)
    count = 0
    for c, n in counter:
        bucket_begin[c], bucket_end[c] = count, count+n-1
        count += n
    return bucket_begin, bucket_end

@njit
def get_lms_array(char_type):
    m = len(char_type)
    lms_array = np.empty(shape=m, dtype=np.bool_)
    lms_array[0] = False
    for i, j in zip(range(m-2, -1, -1), range(m-1, 0, -1)):
        if char_type[j] and not char_type[i]: lms_array[j] = True
        else: lms_array[j] = False 
    return lms_array

@njit
def get_pointer_array(array, lms_array, bucket_end):
    pointer_array = np.full(shape=len(array), fill_value=-1, dtype=np.int32)
    for i, a in enumerate(lms_array):
        if not a: continue
        pointer_array[bucket_end[array[i]]] = i
        bucket_end[array[i]] -= 1
    return pointer_array

@njit
def induceSortL(array, sa, bucket_begin, char_type):
    for i, a in enumerate(sa):
        if a <= 0 or char_type[a-1]:
            continue
        char = array[a-1]
        sa[bucket_begin[char]] = a-1
        bucket_begin[char] += 1

@njit
def induceSortS(array, sa, bucket_end, char_type):
    for i in range(len(sa)-1, -1, -1):
        a = sa[i]
        if not char_type[a-1] or a <= 0:
            continue
        char = array[a-1]
        sa[bucket_end[char]] = a-1
        bucket_end[char] -= 1

@njit
def get_s1(array, pointer_array, lms_array):
    lms_name = np.full(shape=len(array), fill_value=-1, dtype=np.int32)
    pos = pointer_array[0]
    name = 0
    lms_name[pos] = name
    for p in pointer_array[1:]:
        if not lms_array[p]: continue
        if not are_equal_LMSSubstring(array, pos, p, lms_array):
            name += 1
        pos, lms_name[p] = p, name
    s1, substring_pos = List(), List()
    for i, n in enumerate(lms_name):
        if n == -1: continue
        s1.append(n)
        substring_pos.append(i)
    return s1, substring_pos, name+1

@njit
def are_equal_LMSSubstring(array, pos, p, lms_array):
    if pos == len(array)-1: return False
    elif array[pos] != array[p]: return False
    else:
        k = 1
        while True:
            if lms_array[pos+k] and lms_array[p+k]:
                if array[pos+k] == array[p+k]: return True
                else: return False
            elif lms_array[pos+k] is not lms_array[p+k]:
                return False
            elif array[pos+k] != array[p+k]: return False
            k += 1

@njit
def get_sa_from_sa1(array, sa1, substring_pos, bucket_end):
    sa = np.full(shape=len(array), fill_value=-1, dtype=np.int32)
    for i in range(len(sa1)-1, -1, -1):
        a = sa1[i]
        char = array[substring_pos[a]]
        sa[bucket_end[char]] = substring_pos[a]
        bucket_end[char] -= 1
    return sa

@njit
def lcp_array(array, sa):
    m = len(array)
    lcp = np.zeros(shape=m, dtype=np.int32)
    inverse_array = np.empty_like(sa, dtype=np.int32)
    for i, a in enumerate(sa): inverse_array[a] = i
    n = 0
    for i in range(0, m-1):
        k = inverse_array[i]
        j = sa[k-1]
        while array[i+n] == array[j+n]: 
            n += 1
        lcp[k] = n
        if n > 0: n -= 1
    return lcp

if __name__ == '__main__':
    s_ = 'mmiissiissiippii$'
    suffixarray = SuffixArray(s_)
    print(suffixarray.sa)
