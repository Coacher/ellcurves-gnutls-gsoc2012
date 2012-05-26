from wMNAF import *

def print_table_of_values(w = 1):
    print("{0:>5}\t{1:>8}\t{2:>25}\t{3:>25}".format("e", "binary", "NAF", "MNAF"))
    for k in range(1, 40):
        k = Integer(k)
        print("{0:>5}\t{1:>8}\t{2:>25}\t{3:>25}".format(k, base_representation(k, w), wNAF(k, w), wMNAF(k, w)))

def generate_data(P, w, min_value = 10, max_value = 1024, step = 11, printtimings = False):
    if (printtimings):
        print("{0:>10}\t\t{1:>25}\t\t{2:>25}\t\t{3:>5}".format("n", "wMNAF time", "SAGE time", "Check"))

    graph = []
    for k in [1, 2, 3, 5, 7] + range(min_value, max_value, step):
        n = Integer(k)

        t0 = time.clock()
        wnaf = scalar_multiply(n, P, w)
        wnaf_time = time.clock() - t0

        t0 = time.clock()
        default = n*P
        sage_time = time.clock() - t0

        if (sage_time != 0) and (wnaf_time != 0):
            percentage = 100*wnaf_time/sage_time
            graph.append((k, percentage))

        if (printtimings):
            print("{0:>10}\t\t{1:>25}\t\t{2:>25}\t\t{3:>5}".format(k, wnaf_time, sage_time, wnaf == default))

    return graph

#secp128r1
#p = int("FFFFFFFDFFFFFFFFFFFFFFFFFFFFFFFF", 16)
#F = ZZ.quotient_ring(p*ZZ)
#E = EllipticCurve([int("FFFFFFFDFFFFFFFFFFFFFFFFFFFFFFFC", 16), int("E87579C11079F43DD824993C2CEE5ED3", 16)])
E = EllipticCurve([3,5])
w = 4

E.plot()

P = E.rational_points(bound = 15); P

min_value = 100
max_value = 1000
step = 5

points = generate_data(P[0], w, min_value, max_value, step)

G = list_plot(points, plotjoined = True)
F = line([(min_value,100), (max_value,100)], rgbcolor = (0,0.5,0), linestyle="--")
plot(F+G)

points = generate_data(P[0], 3, min_value, max_value, step)

G = list_plot(points, plotjoined = True)
F = line([(min_value,100), (max_value,100)], rgbcolor = (0,0.5,0), linestyle="--")
plot(F+G)

points = generate_data(P[0], 2, min_value, max_value, step)

G = list_plot(points, plotjoined = True)
F = line([(min_value,100), (max_value,100)], rgbcolor = (0,0.5,0), linestyle="--")
plot(F+G)

min_value = 950
max_value = 2000
step = 5

points = generate_data(P[0], 4, min_value, max_value, step)

G = list_plot(points, plotjoined = True)
F = line([(min_value,100), (max_value,100)], rgbcolor = (0,0.5,0), linestyle="--")
plot(F+G)

min_value = 950
max_value = 1200
step = 5

points = generate_data(P[0], 3, min_value, max_value, step)

G = list_plot(points, plotjoined = True)
F = line([(min_value,100), (max_value,100)], rgbcolor = (0,0.5,0), linestyle="--")
plot(F+G)