import sys

load wMNAF.py

l = len(sys.argv)
min_value  = 500
max_value  = 1000
step_value = 10
w = 4
base = 2
print_curve = False

if ("--help" in sys.argv) or ("-h" in sys.argv):
    print("This program will generate output file with runtime comparison")
    print("of wMNAF-based multiplication algorithm against default SAGE's one.")
    print("\nThis program requires at least 3 arguments:")
    print("--min x -- minimum value of multiplier")
    print("--max y -- maximum value of multiplier")
    print("--step z -- multiplier will increase with this step")
    print("\nAlso you can specify optional arguments:")
    print("--w x -- window size. default is 4.")
    print("--base x -- base to use. default is 2. NOT IMPLEMENTED.")
    print("--print_curve -- will output curve in a separate file")
    sys.exit()

if l < 4:
    sys.exit("You must specify at least min, max and step. See --help.")
else:
    if "--min" in sys.argv:
        pos = sys.argv.index("--min")
        min_value = sys.argv[pos + 1]
        sys.argv.remove("--min")
        sys.argv.remove(min_value)
    else:
        sys.exit("You must specify at least min, max and step. See --help.")
    
    if "--max" in sys.argv:
        pos = sys.argv.index("--max")
        max_value = sys.argv[pos + 1]
        sys.argv.remove("--max")
        sys.argv.remove(max_value)
    else:
        sys.exit("You must specify at least min, max and step. See --help.")
    
    if "--step" in sys.argv:
        pos = sys.argv.index("--step")
        step_value = sys.argv[pos + 1]
        sys.argv.remove("--step")
        sys.argv.remove(step_value)
    else:
        sys.exit("You must specify at least min, max and step. See --help.")
    
    if "--w" in sys.argv:
        pos = sys.argv.index("--w")
        w = sys.argv[pos + 1]
        sys.argv.remove("--w")
        sys.argv.remove(w)
    
    if "--base" in sys.argv:
        print("base specifying is not implemented.")
        sys.argv.remove("--base")
    
    if "--print_curve" in sys.argv:
        print_curve = True
        sys.argv.remove("--print_curve")

output_name = "wMNAF_test_{min}-{max}:{step}_{base}^{w}.pdf".format( \
        min = min_value, max = max_value, step = step_value, base = base, w = w)

# elliptic curve to use
E = EllipticCurve([3,5])

# find some points to use
P = E.rational_points(bound = 15)

start = time.clock()

points = generate_data(P[0], int(w), int(min_value), int(max_value), int(step_value), int(base))

if (print_curve):
    Legend = text("\n\nGraphic of {0}".format(str(E)), (0,0), rgbcolor = (0,0,0), \
                   vertical_alignment = "top", horizontal_alignment = "left", axis_coords = True)
    (Legend.plot() + E.plot()).save("curve.pdf")

G = list_plot(points, plotjoined = True)
F = line([(min_value,100), (max_value,100)], rgbcolor = (0,0.5,0), linestyle="--")
Legend = text("\n\nGraphic of wMNAF run time in percents. Sage's time is 100%", (0,0), rgbcolor = (0,0,0), \
               vertical_alignment = "top", horizontal_alignment = "left", axis_coords = True)
(F+G).plot().save(output_name)

finish = time.clock() - start
print("Finished in {t} secs.".format(t = finish))
sys.exit()
