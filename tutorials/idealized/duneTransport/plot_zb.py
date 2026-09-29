"""
Represent the bed elevation at different time
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from fluidfoam import readof as rdf
import os

plt.rcParams["font.size"] = 12
plt.rcParams["lines.linewidth"] = 2.5

pathCase = "./"

timeList = os.popen(f"foamListTimes -case {pathCase}").read().split("\n")[:-1]
ntimes = len(timeList)

colors = ["cornflowerblue", "tomato"]
cm = LinearSegmentedColormap.from_list(
        "Custom", colors, N=ntimes)
colors = cm(np.arange(0, cm.N))

fig, ax = plt.subplots(figsize=(10, 5))

for i, time in enumerate(timeList):
    Xbed, Ybed, Zbed = rdf.readmesh(pathCase, time, boundary="bed")

    ax.plot(Xbed, Zbed, color=colors[i])

ax.set_xlabel(r"$x\,[m]$")
ax.set_ylabel(r"$z_b\,[m^2.s^{-1}]$")

ax.grid()

plt.show()
