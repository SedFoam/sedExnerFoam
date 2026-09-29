"""
compare bed profiles, simulation and experiment
"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from fluidfoam import readof as rdf
import os

plt.rcParams["font.size"] = 15


# color for colorblind
# orange  lightblue green   yellow   darkblue red     pink
# #E69F00 #56B4E9 #009E73 #F0E444 #0072B2 #D55E00 #CC79A7

saveFig = True

tmin, tmax = 0, 40
xmax = 400  # 250
expMarkSize = 50

dataPath = "../DATA/Sandoungout2019/"

plotExp = True
labelExp = "Sandoungout 2019, V = 8.61 mm/s"
LogScale = False
plotLinRegress = True  # False  # plot or not linear regression line

config = 1

if (config == 0):
    caseList = [
        {"path": "../RAS/duneMigration/",
         "color": "#0072B2",
         "marker": "o",
         "markersize": 50,
         "ls": "solid",
         "label": r"2D: suspension on"}
    ]
    figName = "morphoParamsRegLin2D"
elif (config == 1):
    caseList = [
        {"path": "../RAS/duneMigration3D/",
         "color": "#E69F00",
         "marker": "s",
         "markersize": 40,
         "ls": "solid",
         "label": r"3D: suspension on"}
    ]
    figName = "morphoParamsRegLin2D3D"

# experimental parameter
Heq = 16.6e-3  # dune height at equilibrium
Vdune = 8.61e-3  # xh position velocity in m/s
betaRep = 28 * np.pi / 180
t0 = 10


def findXh(Xb, Zb):
    """find middle of downstream slope"""
    izbMax = np.argmax(Zb)
    beta = (Zb[1:] - Zb[:-1]) / (Xb[1:] - Xb[:-1])
    zbMax = Zb[izbMax]
    for i in range(izbMax, len(Zb)):
        if Zb[i] < 0.8 * zbMax:
            break
    itop = i
    for i in range(itop, len(Zb)):
        if Zb[i] < 0.2 * zbMax:
            break
    idown = i
    xtop, xdo = Xb[itop], Xb[idown]
    ztop, zdo = Zb[itop], Zb[idown]
    a, b = np.polyfit((xtop, xdo), (ztop, zdo), deg=1)
    zh = 0.5 * (ztop + zdo)
    xh = (zh - b) / a
    return xh, zh, a


def findL(Xb, Zb):
    """find length of dune"""
    xup = findXup(Xb, Zb)
    xdo = findXdo(Xb, Zb)
    return xdo - xup


def findXup(Xb, Zb):
    """find upstream limit of dune"""
    izbMax = np.argmax(Zb)
    zbMax = Zb[izbMax]
    for i in range(izbMax):
        if Zb[i] >= 0.2 * zbMax:
            break
    i1 = i
    for i in range(i1, izbMax):
        if Zb[i] >= 0.8 * zbMax:
            break
    i2 = i
    # upstream slope angle
    betaUp = (Zb[i2] - Zb[i1]) / (Xb[i2] - Xb[i1])
    xMid = 0.5 * (Xb[i1] + Xb[i2])
    return -0.5 * (Zb[i1] + Zb[i2]) / betaUp + xMid


def findXdo(Xb, Zb):
    """find downstream limit of dune"""
    izbMax = np.argmax(Zb)
    zbMax = Zb[izbMax]
    for i in range(izbMax, len(Xb)):
        if Zb[i] <= 0.8 * zbMax:
            break
    i1 = i
    for i in range(i1, len(Xb)):
        if Zb[i] <= 0.2 * zbMax:
            break
    i2 = i
    # downstream slope angle
    betaDo = (Zb[i2] - Zb[i1]) / (Xb[i2] - Xb[i1])
    xMid = 0.5 * (Xb[i1] + Xb[i2])
    return -0.5 * (Zb[i1] + Zb[i2]) / betaDo + xMid


fig = plt.figure(figsize=(10, 12))
gs = fig.add_gridspec(3, 1, height_ratios=(0.6, 0.4, 0.4))

axXh = fig.add_subplot(gs[0])
axH = fig.add_subplot(gs[1])
axL = fig.add_subplot(gs[2])

if plotExp:
    # xh in function of time, from experiment
    texpXh, xhExp = np.loadtxt(
        dataPath + "xht_u0_43_M0_10g.txt", unpack=True, delimiter=";")
    xhExp -= 77
    axXh.scatter(
        texpXh, xhExp, s=expMarkSize, marker="o",
        color="whitesmoke", edgecolors="black",
        label=labelExp)

    # dune height in function of time, from experiment
    texpH, HtExp = np.loadtxt(
        dataPath + "Ht_u0_43_M0_10g.txt", unpack=True, delimiter=";")
    axH.scatter(
        texpH, HtExp, s=expMarkSize, marker="o",
        color="whitesmoke", edgecolors="black", label=labelExp)
    # dune length in function of time, from experiment
    texpL, LtExp = np.loadtxt(
        dataPath + "L_u0_43_M0_10g.txt", unpack=True, delimiter=";")
    axL.scatter(
        texpL, LtExp, s=expMarkSize, marker="o",
        color="whitesmoke", edgecolors="black",
        label=labelExp)

for i, case in enumerate(caseList):
    pathNum = case["path"]
    color = case["color"]
    marker = case["marker"]
    numMarkSize = case["markersize"]
    ls = case.get("ls")
    label = case.get("label")
    edgecolors = None
    if marker == "*":
        edgecolors = "whitesmoke"
    print(f"\ncase {pathNum}")

    foamTimes = os.popen(
        f"foamListTimes -case {pathNum} -withZero").read()
    numTimeList = foamTimes.split("\n")[:-1]
    print(numTimeList)
    timeArr = np.array([float(t) for t in numTimeList])
    ntimes = len(numTimeList)

    XhArr = np.zeros(ntimes)
    Harr = np.zeros(ntimes)
    Larr = np.zeros(ntimes)

    for j, time in enumerate(numTimeList):
        if float(time) > tmax:
            timeArr = timeArr[:j]
            XhArr = XhArr[:j]
            Harr = Harr[:j]
            Larr = Larr[:j]
            break
        print(f"time = {time} s")
        Xb, Yb, Zb = rdf.readmesh(
            pathNum, time, boundary="bed", verbose=False)
        nx = np.size(Xb)
        if j == 0:
            Zbt = np.zeros((nx, ntimes))
        Zbt[:, j] = Zb[:]
        izbMax = np.argmax(Zb)
        Harr[j] = Zb[izbMax]
        xh, zh, slope = findXh(Xb, Zb)
        XhArr[j] = xh
        Ldune = findL(Xb, Zb)
        Larr[j] = Ldune

    aFit, bFit = np.polyfit(timeArr[1:], XhArr[1:], deg=1) * 1000
    print(f"velocity: {round(aFit, 3)} mm/s")
    label += f", V = {round(aFit, 2)} mm/s"

    # NetCDF file creation
    rootgrp = Dataset(pathNum + "duneParameters.nc", "w")

    # Dimensions creation
    rootgrp.createDimension("scalar", 1)
    rootgrp.createDimension("nt", np.size(timeArr))
    rootgrp.createDimension("nx", nx)

    # Variables creation
    V_file = rootgrp.createVariable("Vdune", np.float64, "scalar")
    b_file = rootgrp.createVariable("b", np.float64, "scalar")
    t0_file = rootgrp.createVariable("t0", np.float64, "scalar")
    #
    time_file = rootgrp.createVariable("time", np.float64, "nt")
    Xh_file = rootgrp.createVariable("Xcrest", np.float64, "nt")
    H_file = rootgrp.createVariable("Hdune", np.float64, "nt")
    L_file = rootgrp.createVariable("Ldune", np.float64, "nt")
    Xb_file = rootgrp.createVariable("Xb", np.float64, "nx")
    Zb_file = rootgrp.createVariable("Zb", np.float64, ("nx", "nt"))
    # Writing variables
    V_file[:] = aFit
    b_file[:] = bFit
    t0_file[:] = t0

    time_file[:] = timeArr[:]
    Xh_file[:] = XhArr[:]
    H_file[:] = Harr[:]
    L_file[:] = Larr[:]
    Xb_file[:] = Xb[:]
    Zb_file[:, :] = Zbt[:, :]

    # File closing
    rootgrp.close()

    # axXh.text(
    #    5, aFit*5 + bFit + 5, rotation=slopeDegree,
    #    s=r"$c_h =$" + f"{round(aFit, 2)}" + r" $mm/s$",
    #    fontsize=12)

    axXh.scatter(
        timeArr - t0, XhArr * 1000, s=numMarkSize, marker=marker,
        color=color, zorder=3., edgecolors=edgecolors, label=label)
    if plotLinRegress:
        axXh.plot(
            timeArr - t0, aFit * timeArr + bFit,
            color=edgecolors, linestyle='--', zorder=1.5)

    axH.scatter(
        timeArr - t0, Harr * 1000, s=numMarkSize, marker=marker,
        color=color, edgecolors=edgecolors, label=label)

    axL.scatter(
        timeArr - t0, Larr * 1000, s=numMarkSize, marker=marker,
        color=color, edgecolors=edgecolors, label=label)

axXh.set_ylabel(r"$x_h\,[mm]$")
if (LogScale is False):
    axXh.set_ylim(0, xmax)
else:
    axXh.set_yscale("log")
    axXh.set_ylim(1e1, xmax)
axXh.legend(loc="lower right", fontsize=12)
axXh.tick_params(bottom=False, labelbottom=False)


axH.set_ylabel(r"$h_d\,[mm]$")
axH.set_ylim(0, 30)
axH.tick_params(bottom=False, labelbottom=False)


axL.set_ylabel(r"$L_d\,[mm]$")
axL.set_ylim(0, 150)
axL.set_xlabel("time [s]")

if plotExp:
    axH.axhline(16.6, color="grey", ls="dashed")
    axL.axhline(116, color="grey", ls="dashed")

for ax in fig.axes:
    ax.set_xlim(0, tmax - t0)
    ax.grid()
    ax.tick_params(
        axis="both", which="major")

fig.tight_layout()

plt.show()

if saveFig:
    fig.savefig("./Figures/" + figName + ".eps", format="eps")
    print("figure saved")
