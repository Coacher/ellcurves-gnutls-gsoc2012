#!/usr/bin/env sage

import sys
from wMNAF import *

l = len(sys.argv)
min_value  = 100
max_value  = 500
step_value = 10
w = 4
base = 2

if ("--help" in sys.argv) or ("-h" in sys.argv):
    print("This program will generate output file with runtime comparison")
    print("of wMNAF-based multiplication algorithm against default SAGE's one.")
    print("\n")
    print("This program requires at least 3 arguments:")
    print("--min x -- minimum value of multiplier")
    print("--max y -- maximum value of multiplier")
    print("--step z -- multiplier will increase with this step")
    print("Also you can specify optional arguments:")
    print("--w x -- window size. defauls is 4.")
    print("--base x -- base to use. defaults is 2. NOT IMPLEMENTED.")
    print("\n")
    exit(EXIT_SUCCESS)

if l < 4:
    print("You must specify at least min, max and step. See --help.")
    exit(EXIT_FAILURE)
else:
    if "--min" in sys.argv:
        pos = sys.argv.index("--min")
        min_value = sys.argv[pos + 1]
        sys.argv.remove("--min")
        sys.argv.remove(min_value)
    else:
        print("You must specify at least min, max and step. See --help.")
        exit(EXIT_FAILURE)
    
    if "--max" in sys.argv:
        pos = sys.argv.index("--max")
        max_value = sys.argv[pos + 1]
        sys.argv.remove("--max")
        sys.argv.remove(max_value)
    else:
        print("You must specify at least min, max and step. See --help.")
        exit(EXIT_FAILURE)
    
    if "--step" in sys.argv:
        pos = sys.argv.index("--step")
        step_value = sys.argv[pos + 1]
        sys.argv.remove("--step")
        sys.argv.remove(step_value)
    else:
        print("You must specify at least min, max and step. See --help.")
        exit(EXIT_FAILURE)
    
    if "--w" in sys.argv:
        pos = sys.argv.index("--w")
        w = sys.argv[pos + 1]
        sys.argv.remove("--w")
        sys.argv.remove(w)
    
    if "--base" in sys.argv:
        print("base specifying is not implemented.")
        sys.argv.remove("--base")

output_name = "wMNAF_test_{min}_{max}_{step}__{w}.pdf".format( \
        min = min_value, max = max_value, step = step_value, w = w)

# elliptic curve to use
E = EllipticCurve([3,5])

# find some points to use
P = E.rational_points(bound = 15)

points = generate_data(P[0], w, min_value, max_value, step_value, base)

E = E.plot()
G = list_plot(points, plotjoined = True)
F = line([(min_value,100), (max_value,100)], rgbcolor = (0,0.5,0), linestyle="--")
Graph = plot(F+G)
(E + Graph).save(output_name)
exit(EXIT_SUCCESS)
