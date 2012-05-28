import sys
from wMNAF_data_collect import *

l = len(sys.argv)
out = "./output.pdf"

if (l == 2):
    out = sys.argv[1]

E = EllipticCurve([3,5])
w = 4

E.plot().save(out)

P = E.rational_points(bound = 15)

points = generate_data(P[0], w, 100, 1000, 5)

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
 
