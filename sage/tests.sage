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
 
