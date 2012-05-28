import sys

load wMNAF.py

EXIT_SUCCESS = 0
EXIT_FAILURE = 1

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
    print("\n")
    print("This program requires at least 3 arguments:")
    print("--min x -- minimum value of multiplier")
    print("--max y -- maximum value of multiplier")
    print("--step z -- multiplier will increase with this step")
    print("Also you can specify optional arguments:")
    print("--w x -- window size. default is 4.")
    print("--base x -- base to use. default is 2. NOT IMPLEMENTED.")
    print("--print_curve -- will output curve in a separate file")
    print("\n")
    sys.exit(EXIT_SUCCESS)

if l < 4:
    print("You must specify at least min, max and step. See --help.")
    sys.exit(EXIT_FAILURE)
else:
    if "--min" in sys.argv:
        pos = sys.argv.index("--min")
        min_value = sys.argv[pos + 1]
        sys.argv.remove("--min")
        sys.argv.remove(min_value)
    else:
        print("You must specify at least min, max and step. See --help.")
        sys.exit(EXIT_FAILURE)
    
    if "--max" in sys.argv:
        pos = sys.argv.index("--max")
        max_value = sys.argv[pos + 1]
        sys.argv.remove("--max")
        sys.argv.remove(max_value)
    else:
        print("You must specify at least min, max and step. See --help.")
        sys.exit(EXIT_FAILURE)
    
    if "--step" in sys.argv:
        pos = sys.argv.index("--step")
        step_value = sys.argv[pos + 1]
        sys.argv.remove("--step")
        sys.argv.remove(step_value)
    else:
        print("You must specify at least min, max and step. See --help.")
        sys.exit(EXIT_FAILURE)
    
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

output_name = "wMNAF_test_{min}_{max}_{step}__{w}.pdf".format( \
        min = min_value, max = max_value, step = step_value, w = w)

# elliptic curve to use
E = EllipticCurve([3,5])

# find some points to use
P = E.rational_points(bound = 15)

points = generate_data(P[0], int(w), int(min_value), int(max_value), int(step_value), int(base))

if (print_curve):
    E.plot().save("curve.pdf")

G = list_plot(points, plotjoined = True)
F = line([(min_value,100), (max_value,100)], rgbcolor = (0,0.5,0), linestyle="--")
(F+G).plot().save(output_name)
sys.exit(EXIT_SUCCESS)
