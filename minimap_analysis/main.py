import sys
import threading
from queue import Queue
import copy
import time
import collections
import pprint

def inv_phi(val):
    """ returns the string associated with val"""
    ret = ""
    def num_to_base(base):
        if base == 1:
            return 'A'
        elif base == 2:
            return 'C'
        elif base == 3:
            return 'G'
        else:
            return '4'
    c = val % 4
    ret.append(num_to_base(c))

def phi(kmer):
    """ converts kmer to base 4 representation"""
    def base_to_num(base):
        if base == 'A':
            return 0
        elif base == 'C':
            return 1
        elif base == 'G':
            return 2
        else:
            return 3

    ret = 0
    for i, base in enumerate(kmer):
        ret += base_to_num(base) * (4 ** i)
    #return ret, 2*len(kmer)
    return ret, 30

def invertible_hash(x, p):
    """ x is a p-bit integer"""
    m = 2**p - 1
    x = (~x + (x << 21)) & m
    x = x ^ x >> 24
    x = (x + (x << 3) + (x << 8)) & m
    x = x ^ x >> 14
    x = (x + (x << 2) + (x << 4)) & m
    x = x ^ x >> 28
    x = (x + (x << 31)) & m
    return x

def minimap_hash(x):
    return invertible_hash(*phi(x))

def reverse_complement(s):
    # TODO
    t = list(range(len(s)))
    for i,char in enumerate(s[::-1]):
        if char == 'T':
            t[i] = 'A'
        if char == 'A':
            t[i] = 'T'
        if char == 'G':
            t[i] = 'C'
        if char == 'C':
            t[i] = 'G'
    return ''.join(t)

def compute_minimizers(s, w, k):
    """
    s = string
    w = window
    k = kmer
    """
    M = []
    t = reverse_complement(s)
    for i in range(len(s) - w - k - 1):
        m = 9999999999999
        for j in range(w):
            u = minimap_hash(s[i+j : i+j+k])
            v = minimap_hash(t[i+j : i+j+k])
            if u != v:
                m = min(m, min(u,v))
        for j in range(w):
            u = minimap_hash(s[i+j : i+j+k])
            v = minimap_hash(t[i+j : i+j+k])
            if u < v and u == m:
                M.append((m, i+j, 0))
            elif v < u and v == m:
                M.append((m, i+j, 1))
    return M


def index(T, w, k):
    """
    T is a list of strings
    """
    H = collections.defaultdict(list)
    for t in T:
        M = compute_minimizers(t, w, k)
        for h,i,r in M:
            H[h].append((t,i,r))
    return H

HashObj = collections.namedtuple('HashObj', ['t','r','c','iprime'], verbose=False)

def MinimapMap(H, q, w, k, g):
    A = []
    ret = []
    M = compute_minimizers(q, w, k)
    for h, i1, r1 in M:
        for t, i2, r2 in H[h]:
            if r1 == r2:
                A.append(HashObj(t, 0, i1 - i2, i2))
            else:
                A.append(HashObj(t, 1, i1 + i2, i2))
    A.sort()
    b = 1
    for e in range(len(A)):
        if False or \
                e == len(A)-1 or \
                A[e+1].t != A[e].t or \
                A[e+1].r != A[e].r or \
                A[e+1].c - A[e].c > g:
            # TODO maximal collinear subset
            # add left and right side to ret
            ret.append((b,e))
            b = e + 1
    return ret

H = index(["CAGTCGAAGTA", "AAAAAAAAA"], 3, 3)
out = MinimapMap(H, "AAAAAAAACA", 3, 3, 10)
pprint.pprint(out)




