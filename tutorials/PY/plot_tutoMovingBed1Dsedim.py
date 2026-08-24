"""
plot repartition of the sediment mass between suspended/deposited
color plot, x axis is time, y axis is elevation
color is suspended sediment volume fraction
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from fluidfoam import readof as rdf
import os


save = True

pathCase = "../laminar/movingBed1DSedim"

rhoS = 2650.  # sediment density
CsMax = 0.6  # maximum sediment volume fraction

# colors for cmap colorplot
sedColor = "peru"  # color for deposited sediment
bgColor = "steelblue"  # color for background

# domain and mesh dimensions
domHeight = 1.
domWidth = 0.1
ncells = len(rdf.readmesh(pathCase, time_name="0")[0])

bedArea = domWidth**2

myCmap = LinearSegmentedColormap.from_list(
    "Custom", (bgColor, sedColor), N=20)

foamTimes = os.popen(f"foamListTimes -case {pathCase} -withZero").read()
timeList = foamTimes.split("\n")[:-1]
ntimes = len(timeList)

Zarr = np.zeros((ntimes, ncells))
Tarr = np.zeros((ntimes, ncells))
CsArr = np.zeros((ntimes, ncells))
Zbed = np.zeros(ntimes)

print("-domain dimensions")
print(f"height {domHeight} m")
print(f"area {round(bedArea, 5)} m2")
print(f"number of cells {ncells}")

massSuspension = np.zeros(ntimes)
massBed = np.zeros(ntimes)

for i, t in enumerate(timeList):
    print(f"\ntime: {t} s")
    Zmesh = rdf.readmesh(pathCase, time_name=t)[2]
    Cs = rdf.readscalar(pathCase, time_name=t, name="Cs")
    if Cs.shape == (1,):
        Cs = np.ones_like(Zmesh) * Cs
    Zbed[i] = rdf.readmesh(pathCase, time_name=t, boundary="bed")[2][0]
    print(f"bed elevation: {Zbed[i]} m")
    try:
        cellVolumes = rdf.readscalar(
            pathCase, time_name=t, name="V", verbose=False)
    except FileNotFoundError:
        os.system(f"postProcess -case {pathCase} -func writeCellVolumes")
        cellVolumes = rdf.readscalar(
            pathCase, time_name=t, name="V", verbose=False)
    # bed is only one face, mesh is 1D
    massBed[i] = rhoS * CsMax * Zbed[i] * bedArea
    massSuspension[i] = rhoS * np.sum(Cs * cellVolumes)
    print("bed mass: ", massBed[i])
    print("suspended mass: ", massSuspension[i])

    Zarr[i, :] = Zmesh[:]
    Tarr[i, :] = np.ones(ncells) * float(t)
    CsArr[i, :] = Cs[:]

totMass = massBed + massSuspension
relMassErr = np.abs(100 * (totMass - totMass[0]) / totMass[0])

# total mass of sediments
print(f"\ninitial mass: {round(totMass[0], 5)} kg, "
      + f"final mass: {round(totMass[-1], 5)} kg")


fig = plt.figure(figsize=(7, 7))
gs = fig.add_gridspec(2, 1)
axCs = fig.add_subplot(gs[0, 0])
axMass = fig.add_subplot(gs[1, 0])

# Color plot of suspended sediment over time
axCs.pcolormesh(
    Tarr, Zarr, CsArr, vmin=0, vmax=1.4*np.max(CsArr),
    shading="gouraud", cmap=myCmap)

axCs.fill_between(Tarr[:, 0], Zbed, np.zeros_like(Zbed), color="peru")
axCs.plot(Tarr[:, 0], Zbed, lw=2, color="black")

axCs.set_xlabel("time [s]", fontsize=12)
axCs.set_ylabel("z [m]", fontsize=12)
axCs.set_ylim(0, 1)
axCs.set_xlim(Tarr[0, 0], Tarr[-1, 0])

# mass repartition suspended/deposited
axMass.plot(
    Tarr[:, 0], totMass / bedArea, ls="solid", lw=1.8,
    color="#0072B2", label="total")
axMass.plot(
    Tarr[:, 0], massSuspension / bedArea, ls="dashed", lw=2,
    color="#009E71", label="suspension")
axMass.plot(
    Tarr[:, 0], massBed / bedArea, ls="dashdot",
    lw=2, color="#D55E00", label="bed")
axMass.set_ylabel(r"$[kg.m^{-2}]$", fontsize=12)
axMass.set_xlabel("time [s]", fontsize=12)
axMass.set_xlim(Tarr[0, 0], Tarr[-1, 0])
axMass.legend(fontsize=12)
axMass.grid()

# error on total mass in %
axErr = axMass.twinx()
axErr.scatter(
    Tarr[:, 0], relMassErr,
    marker="x", color="black")
axErr.set_ylim(
    -0.2*np.max(relMassErr),
    1.7 * np.max(relMassErr))
axErr.set_ylabel("relative mass error in %", fontsize=12)

for ax in fig.axes:
    ax.tick_params(
        axis='both', which='major', labelsize=12)

# change zorder so that relative error is behind mass
axMass.set_zorder(axErr.get_zorder()+1)
axMass.set_frame_on(False)

fig.tight_layout()
plt.show()

if save:
    imname = "movingBed1Dsedim_tuto.png"
    fig.savefig("./Figures/" + imname, dpi=300, format="png", transparent=True)
