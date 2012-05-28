#!/usr/bin/env python2

# needed for data collecting
import time
from sage import *

# helper functions
def digit_to_char(digit):
    if (digit < 10):
        return str(digit)
    if (digit > 35):
        raise ValueError("Digit must not be greater than 35.")
    return chr(digit - 10 + ord('A'))

def char_to_digit(char):
    if ( (char < '0') or (char > 'Z') or ((char > '9') and (char < 'A')) ):
        raise ValueError("Char must be a digit from 0 to 9 or a letter from A to Z.")
    if (char <= '9'):
        return int(char)
    return int(ord(char) - ord('A') + 10)

def base_representation(x, w = 1, base = 2):
    signed = False
    b = base ** w

    if (x < 0):
        signed = True
        x = abs(x)

    div = x // b
    digits = str(digit_to_char(x % b))

    while (div != 0):
        digits = str(digit_to_char(div % b)) + digits
        div //= b
    if (signed):
        digits = "-" + digits
    return digits

def LSB(x, w = 1, base = 2, bits = 1):
    return (x % (base ** (w*bits)))

def LSB_to_string(x, w = 1, base = 2, bits = 1):
    '''String representation of LSB.

    Returned string will be filled with leading zeros and will have exactly `bits` length.'''
    ret = base_representation(LSB(x, w, base, bits), w, base)
    if (len(ret) < bits):
        ret = "0"*(bits - len(ret)) + ret
    return ret

def get_bit(x, bitnum, w = 1, base = 2):
    '''Returns `bitnum` bit from `x`. Bits counting from 0.''' 
    x //= base ** (w * bitnum)
    return x % (base ** w)

# basic algorithms
def wNAF_openssl(x, w = 1, base = 2):
    '''wNAF implementation used in OpenSSL. Only works with base == 2.'''
    if (x == 0):
        return "0"

    bit = base ** w
    next_bit = bit * base
    mask = next_bit - 1

    sign = 1
    if (x < 0):
        sign = -1
        x = abs(x)

    len = x.nbits()
    ret = ""

    window_val = LSB(x, w + 1, base)
    j = 0

    while (window_val != 0) or (j + w + 1 < len):
        digit = 0
        if (LSB(window_val) & 1):
            if (window_val & bit):
                digit = window_val - next_bit
            else:
                digit = window_val
            window_val -= digit

        ret = [sign*digit] + ret
        j += 1

        window_val //= base
        window_val += get_bit(x, w + j) * bit
    return ret

def wMNAF_openssl(x, w = 1, base = 2):
    '''wMNAF implementation used in OpenSSL. Only works with base == 2.'''
    if (x == 0):
        return "0"

    bit = base ** w
    next_bit = bit * base
    mask = next_bit - 1

    sign = 1
    if (x < 0):
        sign = -1
        x = abs(x)

    len = x.nbits()
    ret = ""

    window_val = LSB(x, w + 1, base)
    j = 0

    while (window_val != 0) or (j + w + 1 < len):
        digit = 0
        if (LSB(window_val) & 1):
            if (window_val & bit):
                digit = window_val - next_bit
                if (j + w + 1 >= len):
                    digit = window_val & (mask // base)
            else:
                digit = window_val
            window_val -= digit

        ret = [sign*digit] + ret
        j += 1

        window_val //= base
        window_val += get_bit(x, w + j) * bit
    return ret

def wNAF(x, w = 1, base = 2):
    '''wNAF implementation.'''
    if (x == 0):
        return "0"

    sign = 1
    if (x < 0):
        sign = -1
        x = abs(x)

    ret = []

    while (x > 0):
        if (LSB(x, 1, base)):
            digit = LSB(x, w + 1, base)
            if (digit >= base**w):
                digit -= base**(w + 1)
            x -= digit
        else:
            digit = 0
        ret = [sign*digit] + ret
        x //= base

    return ret

def wMNAF(x, w = 1, base = 2):
    '''wMNAF implementation.'''
    if (x == 0):
        return "0"

    sign = 1
    if (x < 0):
        sign = -1
        x = abs(x)

    mask = (base ** w) - 1
    ret = []
    j = 0
    len = x.nbits()

    while (x > 0) or (j + w + 1 < len):
        if (LSB(x, 1, base)):
            digit = LSB(x, w + 1, base)
            if (digit >= base**w):
                digit -= base**(w + 1)
                if (j + w + 1 >= len):
                    digit &= mask
            x -= digit
        else:
            digit = 0
        ret = [sign*digit] + ret
        x //= base
        j += 1

    return ret

# multiplication algorithm
def precompute_values(P, w = 1, base = 2):
    '''Precomputes multiples of P.

    We need only values which are not multiples of base.'''
    pos = []
    neg = []
    small_range = range(1, base)
    for i in range(0, base**(w), base):
        for j in small_range:
            pos.append( (i+j)*P)
            neg.append(-(i+j)*P)
    return pos,neg

def scalar_multiply(n, pos, neg, w = 1, base = 2, representation = wMNAF):
    form = representation(n, w, base)
    Q = 0
    for digit in form:
        Q = base*Q
        if (digit != 0):
            if (digit > 0):
                Q += pos[ (digit // base)*(base - 1) + (digit % base) - 1]
            else:
                digit = -digit
                Q += neg[ (digit // base)*(base - 1) + (digit % base) - 1]
    return Q

# data collection tools 
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

def generate_data(P, w, min_value = 100, max_value = 1024, step = 10, base = 2, with_timings = False):
    if (with_timings):
        print("{0:>10}\t\t{1:>25}\t\t{2:>25}\t\t{3:>5}".format("n", "wMNAF time", "SAGE time", "Valid"))

    graph = []
    error = False

    neg, pos = precompute_values(P, w, base)

    for k in [1, 2, 3, 5, 7] + range(min_value, max_value, step):
        n = Integer(k)

        t0 = time.clock()
        wnaf = scalar_multiply(n, neg, pos)
        wnaf_time = time.clock() - t0

        t0 = time.clock()
        default = n*P
        sage_time = time.clock() - t0

        error = not (wnaf == default)
        if error:
            print("wMNAF produced incorrect results!")
            print("multiplier was {k}, window size was {w}, base was {base}".format(k = k, w = w, base = base))
            print("point was {point}".format(point = P))
            raise ArithmeticError()

        if (sage_time != 0) and (wnaf_time != 0):
            percentage = 100*wnaf_time/sage_time
            graph.append((k, percentage))

        if (with_timings):
            print("{0:>10}\t\t{1:>25}\t\t{2:>25}\t\t{3:>5}".format(k, wnaf_time, sage_time, str(not error)))

    return graph
