#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 13:39:18 2021

@author: zev
"""

from Bio import SeqIO

def hamming_distance(error, correct, correct_rc):
    if sum(e != c for e, c in zip(error, correct)) == 1: return 0
    elif sum(e != c for e, c in zip(error, correct_rc)) == 1: return 1
    else: return

def error_correction(reads):
    reads_rc = [r.reverse_complement() for r in reads]
    hash_values, correct_reads = set(), set()
    for i, read in enumerate(reads):
        if hash(read) not in hash_values and hash(reads_rc[i]) not in hash_values:
            hash_values.add(hash(read))
        else: correct_reads.add(i)

    error_reads = set(range(len(reads))) - correct_reads
    for i in error_reads:
        for j in correct_reads:
            ret = hamming_distance(reads[i], reads[j], reads_rc[j])
            if ret is None: continue
            elif ret == 0: print('{}->{}'.format(reads[i], reads[j]))
            else: print('{}->{}'.format(reads[i], reads_rc[j]))
            break

if __name__ == '__main__':
    with open('rosalind_corr.txt', encoding='utf-8') as handle:
        record = SeqIO.parse(handle, 'fasta')
        reads_ = [r.seq for r in record]
    error_correction(reads_)
