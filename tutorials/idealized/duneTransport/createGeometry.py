"""
create initial dune geometry

The dune has a Gaussian shape
"""

import numpy as np
import os

path = "./constant/polyMesh"

Hz = 1.  # domain height
Hd = 0.1  # dune height

x0d = 2.  # initial x-coordinates of cone center
sigmaD = 0.1  # width of dune

newlines = []


def getZb(x):
    zb = Hd * np.exp(-((x-x0d)/sigmaD)**2)
    return zb


os.system(f"cp {path}/points {path}/newPoints")

newPoints = open(f"{path}/newPoints", "w")


with open(f"{path}/points", "r") as f:
    for line in f:
        newline = line
        if line[0] == "(" and line[-2:] == ")\n":
            xs, ys, zs = line[:-1].strip("()").split(" ")
            x, z = float(xs), float(zs)
            zb = getZb(x)
            z = zb + z * (1 - zb / Hz)
            zs = str(round(z, 10))
            newline = f"({xs} {ys} {zs})\n"
            print(z)
        newPoints.write(newline)

newPoints.close()

os.system(f"rm {path}/points && mv {path}/newPoints {path}/points")
