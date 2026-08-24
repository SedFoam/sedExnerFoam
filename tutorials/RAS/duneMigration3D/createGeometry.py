"""
create initial dune geometry
"""

import numpy as np
import os

path = "./constant/polyMesh"

Hz = 0.09  # domain height
Hd = 0.0238  # cone height
beta = 28.  # cone slope angle in degrees
Xd = 0.  # initial x-coordinates of cone center
dxd = 0.01  # width on top of dune where to smooth it

newlines = []

betaRad = beta * np.pi / 180
slope = np.tan(betaRad)

# coef for polynomial, shape of top of dune
ap = -np.tan(betaRad) / dxd
bp = (2 * np.tan(betaRad) * Xd) / dxd
cp = Hd - 0.25*dxd*np.tan(betaRad) - np.tan(betaRad) * Xd**2 / dxd


def polCrest(x):
    return ap*x**2 + bp*x + cp


def getZb(x):
    zb = np.where(
        x > Xd, Hd - (x-Xd)*slope, Hd + (x-Xd)*slope)
    zb = np.where(zb > 0, zb, 0)
    zb = np.where(
        x < Xd+0.5*dxd and x > Xd-0.5*dxd,
        polCrest(x), zb
    )
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
        newPoints.write(newline)

newPoints.close()

os.system(f"rm {path}/points && mv {path}/newPoints {path}/points")
