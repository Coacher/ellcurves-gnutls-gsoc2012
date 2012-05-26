#!/usr/bin/env python2

import time
from wMNAF import *

def print_table_of_values(min = 1, max = 64, w = 1, base = 2):
    '''Print out the comparison table of wNAF and wMNAF.'''
    print("{int:>5} {basew:>10} {basewlen:>10} {wnaf:>25} {wnaflen:>8} {wmnaf:>25} {wmnaflen:>8}".format(\
    int = "e", basew = "base**w", basewlen = "baselen", \
    wnaf = "wNAF", wnaflen="wNAFlen", wmnaf="wMNAF", wmnaflen="wMNAFlen"))

    for k in range(min, max):
        k = Integer(k)
        basew = base_representation(k, w, base)
        wnaf = wNAF(k, w, base)
        wmnaf = wMNAF(k, w, base)
        print("{int:>5} {basew:>10} {basewlen:>10} {wnaf:>25} {wnaflen:>8} {wmnaf:>25} {wmnaflen:>8}".format(\
        int = k, basew = basew, basewlen = len(basew), \
        wnaf = wnaf, wnaflen = len(wnaf), wmnaf = wmnaf, wmnaflen = len(wmnaf)))

def generate_data(P, w, min_value = 10, max_value = 1024, step = 11, base = 2, with_timings = False):
    if (with_timings):
        print("{0:>10}\t\t{1:>25}\t\t{2:>25}\t\t{3:>5}".format("n", "wMNAF time", "SAGE time", "Valid"))

    graph = []

    neg, pos = precompute_values(P, w, base)

    for k in [1, 2, 3, 5, 7] + range(min_value, max_value, step):
        n = Integer(k)

        t0 = time.clock()
        wnaf = scalar_multiply(n, neg, pos)
        wnaf_time = time.clock() - t0

        t0 = time.clock()
        default = n*P
        sage_time = time.clock() - t0

        if (sage_time != 0) and (wnaf_time != 0):
            percentage = 100*wnaf_time/sage_time
            graph.append((k, percentage))

        if (with_timings):
            print("{0:>10}\t\t{1:>25}\t\t{2:>25}\t\t{3:>5}".format(k, wnaf_time, sage_time, str(wnaf == default)))

    return graph
