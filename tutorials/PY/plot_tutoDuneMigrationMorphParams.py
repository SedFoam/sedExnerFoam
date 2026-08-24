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

config = 0

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
        {"path": "../RAS/duneMigration/",
         "color": "#0072B2",
         "marker": "o",
         "markersize": 50,
         "ls": "dotted",
         "label": r"2D: suspension on"},
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

    # Numerical data reading
    num_file = Dataset(os.path.join(pathNum, "duneParameters.nc"))
    t0 = num_file.variables["t0"][:]

    timeArr = num_file.variables["time"][:]
    XhArr = num_file.variables["Xcrest"][:]
    Harr = num_file.variables["Hdune"][:]
    Larr = num_file.variables["Ldune"][:]
    # File closing
    num_file.close()

    aFit, bFit = np.polyfit(timeArr[1:], XhArr[1:], deg=1) * 1000
    print(f"velocity: {round(aFit, 3)} mm/s")
    label += f", V = {round(aFit, 2)} mm/s"
    # axXh.text(
    #    5, aFit*5 + bFit + 5, rotation=slopeDegree,
    #    s=r"$c_h =$" + f"{round(aFit, 2)}" + r" $mm/s$",
    #    fontsize=12)

    axXh.scatter(
        timeArr - t0, XhArr * 1000, s=numMarkSize, marker=marker,
        color=color, zorder=3., edgecolors=edgecolors, label=label)
    if plotLinRegress:
        axXh.plot(
            timeArr - t0, aFit * timeArr + bFit, color=edgecolors,
            linestyle='--', zorder=1.5)

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
# axH.legend(fontsize=12)
# axH.set_ylim(15, 25)
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
