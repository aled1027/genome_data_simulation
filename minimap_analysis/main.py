import sys
import threading
from queue import Queue
import time


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
    return ret, 2*len(kmer)

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

minimap_hash = lambda x : invertible_hash(*phi(x))

#def go(x, y):
#    for i in range(x, y):
#        print(i, invertible_hash(i, 32))
#
#go(2**21, 2**22)

#for z in range(24,28):
#    print("z = {}".format(z))
#    size = int(((2**(z+1)) - (2**(z))) / 4)
#    #q = Queue()
#    for i in range(4):
#        x = 2**z + (i * size)
#        y = 2**z + (size * (i + 1))
#        go(x,y)
#        #t = threading.Thread(target=go, args=[x, y])
#        #t.start()
# 2 seconds for z = 16
# 4.5 seconds z = 17
# 8.62 seconds for z = 18
# 14.9 for z = 19
