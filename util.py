import sys

def error(comment):
    sys.exit("Error: " + comment)

def print_sys_version():
    print("sys.version", sys.version)

def read_string(filename):
    print("read", filename)
    s = []
    f = open(filename, "r")
    for line in f:
       v = line.split()
       s.append(v[0])
    f.close()
    return s

def read_float(filename, i = 0):
    print("read", filename)
    x = []
    f = open(filename, "r")
    for line in f:
       v = line.split()
       x.append(float(v[i]))
    f.close()
    return x

def read_float2(filename, i = 0, j = 1):
    print("read", filename)
    x = []
    y = []
    f = open(filename, "r")
    for line in f:
       v = line.split()
       x.append(float(v[i]))
       y.append(float(v[j]))
    f.close()
    return x, y

def save_vector(x, filename):
    out = open(filename, "w")
    for i in range(len(x)):
        print(x[i], file = out)
    out.close()
    print(filename)

def save_vector2(x, y, filename):
    out = open(filename, "w")
    for i in range(len(x)):
        print(x[i], y[i], file = out)
    out.close()
    print(filename)

