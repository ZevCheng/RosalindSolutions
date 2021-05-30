#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 15:19:47 2021

@author: zev
"""

import numpy as np

def assessing_assembly_quality(reads):
    length = sorted((len(read) for read in reads), reverse=True)
    array = np.array(length, dtype=np.int32)
    cumulative = array.cumsum()
    pos50 = np.where(cumulative >= array.sum()*0.5)[0][0]
    pos75 = np.where(cumulative >= array.sum()*0.75)[0][0]
    return length[pos50], length[pos75]

if __name__ == '__main__':
    with open('rosalind_asmq.txt', encoding='utf-8') as handle:
        reads_ = handle.read().split()
    print(assessing_assembly_quality(reads_))
