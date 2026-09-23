"""
Represent the bed elevation at different time
"""

import numpy as np
import matplotlib.pyplot as plt
from fluidfoam import readof as rdf
import os

pathCase = "./"

timeList = os.popen(f"foamListTimes -case {pathCase}").read().split("\n")[:-1]

fig, ax = plt.subplots()

for time in timeList:
    Xbed, Ybed, Zbed = rdf.readmesh(pathCase, time, boundary="bed")

    ax.plot(Xbed, Zbed)

plt.show()