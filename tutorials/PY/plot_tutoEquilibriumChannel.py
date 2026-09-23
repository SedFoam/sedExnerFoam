"""
Plot velocity, Reynolds stress and C profiles from available Lyn configurations
"""

import numpy as np
import matplotlib.pyplot as plt
from fluidfoam import readof as rdf
from scipy.optimize import minimize
import sys
import os
from readExpData import readExpParameters


time = "latestTime"

saveFig = False

pathCases = "../RAS/equilibriumChannel/"

colorDict = {
    "1565": "#009E73",
    "1965": "#D55E00",
    "2565": "#0072B2",
    "1957": "#CC79A7"
}

lineStyleDict = {
    "1565": "solid",
    "1965": "dashed",
    "2565": "dashdot",
    "1957": "solid"
}

expMarkerDict = {"1565": "P", "1965": "^", "2565": "x", "1957": "s"}

caseList = os.popen(f"ls {pathCases}").read().split("\n")[:-1]

dirToPop = []

for name in caseList:
    if name[:4] != "case" or len(name) > 10:
        dirToPop.append(name)

for name in dirToPop:
    caseList.pop(caseList.index(name))

pathData = "../DATA/dataLyn1988/"


def getUmean(U, Z):
    """"""
    dz0 = 2*Z[0]
    dZ = [dz0]
    for i in range(len(Z)-1):
        dz = (Z[i+1] - Z[i] - 0.5*dZ[i]) * 2
        dZ.append(dz)
    Umean = sum([u*dz for u, dz in zip(U, dZ)])/Hwater
    return Umean


def RouseProfile(z, z0, c0, Ro):
    """Rouse equilibrium profile for suspended sediments

    Parameters
        z: float
            distance from wall
        z0: float
            reference level close to the bottom
        c0: float
            concentration at level z0
        Ro: float
            Rouse number
    """
    Crouse = ((z0 / (Hwater - z0)) * ((Hwater - z)/z))**Ro
    return c0 * Crouse


fig, (axU, axC) = plt.subplots(
    nrows=1, ncols=2, figsize=(8.3, 6))

for i, caseName in enumerate(caseList):
    expCase = caseName[4:]  # experiment from Lyn 1988
    print("\n\n - case: ", expCase)

    # read experimental data
    case = readExpParameters(
        pathData + "parametersExpLyn.txt", caseName=expCase)
    print(expCase)
    # marker for scatter plot experimental data
    expMarker = expMarkerDict[expCase]
    col = colorDict[expCase]
    ls = lineStyleDict[expCase]

    # simulation parameters
    Hwater = case["Hwater"]  # water depth
    dS = case["dS"]  # sand particle diameter (m)
    ufExp = case["uf"]  # friction velocity from experiments
    UexpMean = case["Umean"]  # experimental mean velocity
    wsExp = case["ws"]  # experimental settling velocity
    nuF = 1e-6  # fluid kinematic viscosity
    rhoS = 2650.  # sediment density
    rhoF = 1000.  # water density
    CsMax = 0.57  # maximum sediment volumic fraction
    g = 9.81  # gravity acceleration
    Sc = 1.  # Schmidt number
    kappa = 0.4  # von Karman constant
    ks = 2.5 * dS  # Nikuradse equivalent roughness
    critShields = 0.047  # critical shields number

    print("- experimental data")
    print(f" sediment: diameter = {dS*100} mm, density = {rhoS} kg.m-3, "
          + f"settling velocity = {wsExp} m/s")

    print(f" hydraulic: water depth = {Hwater} m, mean velocity = {UexpMean}, "
          + f"friction velocity = {ufExp} m/s")
    # load experiment results
    pathExpU = pathData + "Uuf" + expCase + ".txt"
    pathExpC = pathData + "Ceq" + expCase + ".txt"
    ZpExp, Uufexp = np.loadtxt(pathExpU, delimiter=";", unpack=True)
    ZcExp, Cexp = np.loadtxt(pathExpC, delimiter=";", unpack=True)

    # load results from simulation
    pathOF = pathCases + caseName
    print(pathOF)
    Zmesh = rdf.readmesh(pathOF, verbose=False)[2]
    # x component of velocity field
    UxField = rdf.readfield(pathOF, time, "U", verbose=False)[0]

    # concentration Cs
    CsField = rdf.readscalar(pathOF, time, "Cs", verbose=False)
    if CsField.shape == (1,):
        CsField = CsField * np.zeros_like(Zmesh)
    # turbulent eddy viscosity
    nutField = rdf.readscalar(pathOF, time, "nut", verbose=False)

    shields = rdf.readvector(pathOF, time, "shieldsVf",
                             boundary="bed", verbose=False)[0, 0]

    ufOf = np.sqrt(shields * (rhoS/rhoF - 1)*g*dS)
    print("- simulation results")
    print(f" Shields number = {shields}")
    print(f" friction velocity from simulation: {ufOf} m/s")
    kPlus = ks * ufOf / nuF  # dimension less roughness length
    print(f" bed equivalent roughness height, ks = {ks} m")
    print(f" dimension less roughness, k+ = {kPlus}")
    print(f" z+ (z * uf / nu) = {Zmesh[0]*ufOf/nuF}")
    print(f" z1/ks = {Zmesh[0]/ks}")
    print(" simulation mean velocity: ", getUmean(UxField, Zmesh), "m/s")

    # erosion deposition condition
    Dstar = dS * ((rhoS/rhoF - 1) * g / nuF**2)**(1/3)
    Tvr = (shields/critShields - 1)
    nut1 = max(nutField[0], nuF)  # kinematic viscosity in first cell
    suspNum1 = wsExp * Sc / nut1
    zrefdzb = Zmesh[0] - ks
    CbRef = 0.015 * (dS / ks) * Tvr**1.5 * Dstar**(-0.3)
    Cbed = CsField[0] + CbRef * (1 - np.exp(-suspNum1*zrefdzb))
    # CbRef extrapolated at level z1
    CbRef1 = min(CbRef * np.exp(-suspNum1 * zrefdzb), CsMax)
    print(f"- Erosion Deposition boundary condition")
    print(f" D* = {Dstar}")
    print(f" T = {Tvr}")
    print(f" nut first cell = {nut1}")
    print(f" ws*Sc/nut = {suspNum1}")
    print(f" z1 - dzb = {zrefdzb}")
    print(f" c1 = {CsField[0]}")
    print(f" Cb* = {CbRef}")
    print(f" Cb = {Cbed}")
    print(f" CbRef1 = {CbRef1}")

    # plot velocity
    axU.scatter(
        Uufexp, ZpExp * nuF / (ufExp * ks), s=50,
        marker=expMarker, color=col, label=f"{expCase} exp")
    axU.plot(UxField/ufOf, Zmesh/ks, color=col, ls=ls, label=f"{expCase} sim")

    # plot suspension
    axC.scatter(
        Cexp, ZcExp, marker=expMarker, s=50,
        color=col, label=f"{expCase} exp")
    axC.plot(
        CsField, Zmesh/Hwater, ls=ls, color=col, label=f"{expCase} sim")

    axU.set_xlabel(r"$u/u_*$", fontsize=15)
axU.set_ylabel(r"$z/k_s$", fontsize=15)
axU.set_yscale("log")
axU.set_xlim(0, 25)
axU.set_ylim(0.1, None)
axU.tick_params(
    axis='both', which='major', labelsize=15)
axU.grid()

axC.legend(fontsize=14)
axC.set_xlabel(r"$c_s$", fontsize=15)
axC.set_ylabel("z/H", fontsize=15)
axC.set_xscale("log")
axC.set_xlim(1e-7, None)
axC.set_ylim(0, 1)
axC.tick_params(
    axis='both', which='major', labelsize=15)
axC.grid()


fig.tight_layout()

plt.show()

if saveFig:
    imname = "lyn1988_U_Cs_results.eps"
    fig.savefig("./" + imname, format="eps")
