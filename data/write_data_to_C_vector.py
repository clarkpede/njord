import numpy as np
import sys

data = np.genfromtxt("DriverSeegmiller.dat", skip_header=952, usecols = [1,2,3],
                     max_rows=27)
data = np.vstack([[0.0,0.0,0.0],data,[9.0,0,0]])

for j in range(3):
    sys.stdout.write("\n{ ")
    for i in range(29):
        sys.stdout.write("{:.3f}".format(data[i,j]))
        if (i < 28):
            sys.stdout.write(", ")
    sys.stdout.write("}")
sys.stdout.write("\n")
