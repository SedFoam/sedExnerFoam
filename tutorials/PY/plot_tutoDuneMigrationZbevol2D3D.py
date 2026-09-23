"""
compare bed profiles, simulation and experiment
"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from fluidfoam import readof as rdf
import numpy as np
import os

plt.rcParams["font.size"] = 15
plt.rcParams["lines.linewidth"] = 1.5

# color for colorblind
# orange  lightblue green   yellow   darkblue red     pink
# #E69F00 #56B4E9   #009E73 #F0E444 #0072B2   #D55E00 #CC79A7

saveFig = True

dataPath = "../DATA/Sandoungout2019/"

xmin, xmax = -100, 250

plotExp = True
plotInitialCondition = True

t0 = 10

config = 1

if (config == 0):
    caseList = [
        {"path": "../RAS/duneMigration/",
         "color": "#0072B2",
         "ls": "solid",
         "label": r"suspension on"}
    ]
    figName = "zb_evolution"
elif (config == 1):
    caseList = [
        {"path": "../RAS/duneMigration/",
         "color": "#0072B2",
         "ls": "dotted",
         "label": r"2D: suspension on"},
        {"path": "../RAS/duneMigration3D/",
         "color": "#0072B2",
         "ls": "solid",
         "label": r"3D: suspension on"}
    ]
    figName = "zb_evolution_comparison2D3D"

ncases = len(caseList)
colors = ["cornflowerblue", "tomato"]
cm = LinearSegmentedColormap.from_list(
    "Custom", colors, N=ncases)
colors = cm(np.arange(0, cm.N))

# time to represent
timeSet = set(["0", "4", "8", "12"])
# timeSet = set(["0", "4", "8", "14"])

expPrefix = "zb_u0_43_M0_10g_"
expFileList = os.popen(
    f"ls " + dataPath + expPrefix + "*.txt").read().split("\n")[:-1]
ntimeExp = len(expFileList)
expTimeList = [expFile.strip("s.txt").split("_")[-1]
               for expFile in expFileList]
expTimeSet = set(expTimeList).intersection(timeSet)


# parameters for initial dune position
Hd = 0.0238  # cone height
betaDeg = 28.  # cone slope angle in degrees
Xd = 0.  # initial x-coordinates of cone center
dxd = 0.01  # width on top of dune where to smooth it
betaRad = betaDeg * np.pi / 180
slope = np.tan(betaRad)
# coef for polynomial, shape of top of dune
ap = -np.tan(betaRad) / dxd
bp = (2 * np.tan(betaRad) * Xd) / dxd
cp = Hd - 0.25 * dxd * np.tan(betaRad) - np.tan(betaRad) * Xd**2 / dxd


def polCrest(x):
    return ap * x**2 + bp * x + cp


def zb0(x):
    zb = np.where(
        x > Xd, Hd - (x - Xd) * slope, Hd + (x - Xd) * slope)
    zb = np.where(zb > 0, zb, 0)
    zb = np.where(
        np.abs(x) < Xd + 0.5 * dxd, polCrest(x), zb)
    return zb


fig, axZb = plt.subplots(figsize=(14, 4))

for i, case in enumerate(caseList):
    pathNum = case["path"]
    color = case["color"]
    # color = colors[i]
    ls = case.get("ls")
    lw = 1.5
    if ls == "dotted":
        lw = 2.

    labelExp = None
    if i == 0:
        labelExp = "Sandoungout 2019"

    # Numerical data reading
    num_file = Dataset(os.path.join(pathNum, "duneParameters.nc"))
    t0 = num_file.variables["t0"][:]

    timeArr = num_file.variables["time"][:]
    Xb = num_file.variables["Xb"][:]
    Zb = num_file.variables["Zb"][:, :]
    # File closing
    num_file.close()

    numTimeList = list(timeArr)
    dumList = []
    for t in numTimeList:
        dumList.append(str(int(float(t) - t0)))
    numTimeList = dumList
    thisTimeSet = set(numTimeList).intersection(timeSet)

    if i == 0:
        print("exp times:", " ".join(expTimeList))
        # plot initial condition
        if plotInitialCondition:
            X0 = np.linspace(-0.05, 0.05, 200)
            axZb.plot(
                X0 * 1000, zb0(X0) * 1000, color="black",
                label="initial condition")
    print("simulation times:", " ".join(numTimeList))
    print("time set:", thisTimeSet)

    for j, time in enumerate(thisTimeSet):

        # experimental data
        if plotExp:
            if j > 0:
                labelExp = None
            expFilePath = dataPath + expPrefix + time + "s.txt"
            xexp, zbexp = np.loadtxt(
                expFilePath, unpack=True, delimiter=";")
            xexp -= 77.
            axZb.scatter(xexp, zbexp, color="black", label=labelExp)
            if i == 0:
                izbMax = np.argmax(zbexp)
                xtext, ztext = xexp[izbMax] - 5, zbexp[izbMax] + 2
                axZb.text(
                    xtext, ztext,
                    f"t = {time} s")
        # OpenFoam data
        label = case["label"]
        if j == 0:
            label = case.get("label")
        timeNumRead = str(t0 + int(time))
        index = np.where(timeArr[:].data == int(t0) + int(time))[0][0]
        # the next line is for directly reading the openfoam result
        # Xb, Yb, Zb = rdf.readmesh(pathNum, timeNumRead, boundary="bed")
        if time == "0":
            axZb.plot(
                Xb * 1000, Zb[:, 0] * 1000, color=color, ls=ls,
                lw=lw, label=label)
        else:
            axZb.plot(
                Xb * 1000, Zb[:, index] * 1000, color=color, ls=ls, lw=lw)
        if not plotExp and i == 0:
            izbMax = np.argmax(Zb)
            xtext, ztext = 1000 * Xb[izbMax] + 5, 1000 * Zb[izbMax] + 2
            axZb.text(
                xtext, ztext, f"t = {time} s")

handles, labels = axZb.get_legend_handles_labels()

axZb.set_xlabel("x [mm]")
axZb.set_ylabel(r"$z_b$ [mm]")
axZb.set_xlim(xmin, xmax)
axZb.set_ylim(None, 30)
if (config == 0):
    axZb.legend(handles[0:4], labels[0:4], loc="upper left")
else:
    axZb.legend(handles[0:6], labels[0:6], loc="upper left")
axZb.tick_params(
    axis="both", which="major")
axZb.grid()

plt.tight_layout()
plt.show()

if saveFig:
    fig.savefig("Figures/" + figName + ".png", dpi=300, format="png")
    fig.savefig("Figures/" + figName + ".eps", dpi=300, format="eps")
    print("figure saved")
