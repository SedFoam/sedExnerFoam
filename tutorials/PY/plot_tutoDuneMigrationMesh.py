"""
colorplot and streamplot of velocity fields with the dune at different times
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.collections import LineCollection
from fluidfoam import readof as rdf
from fluidfoam import MeshVisu

plt.rcParams["font.size"] = 15

case = '2D'

if case == '2D':
    pathCase = "../RAS/duneMigration"
    figName = "Mesh_dune2D"

    xmin, xmax = -0.075, 0.125
    zmin, zmax = -0.005, 0.095

    # box to zoom in on for mesh display:
    mybox = ((xmin, -0.1, 0), (xmax, 0.1, 0.09))
else:
    pathCase = "../RAS/duneMigration3D"
    figName = "Mesh_dune3D"

    xmin, xmax = -0.075, 0.125
    ymin, ymax = -0.001, 0.007
    zmin, zmax = -0.005, 0.095

    # box to zoom in on for mesh display:
    mybox = ((xmin, ymin, 0), (xmax, ymax, 0.09))

saveFig = True

timeList = ["10", "15"]  # at least two times must be given
ntimes = len(timeList)

# box: tuple of box’s dimension: ((xmin, ymin, zmin), (xmax, ymax, zmax))

fig = plt.figure(figsize=(12, 4 * ntimes))
gs = fig.add_gridspec(ntimes, 1)

for i, time in enumerate(timeList):
    ax = fig.add_subplot(gs[i, 0])
    print(f"\n- time: {time} s")
    # readmesh, cell centers and bed faces
    # Xmesh, Ymesh, Zmesh = rdf.readmesh(pathCase, time_name=time)
    # Xbed, Ybed, Zbed = rdf.readmesh(pathCase, boundary="bed", time_name=time)
    # read velocity field
    # Load mesh and create an object called myMesh
    # The box by default is egal to the mesh dimension
    myMesh = MeshVisu(path=pathCase, time_name=time, plane='xz', box=mybox)

    # create a collection with edges and print it
    ln_coll = LineCollection(myMesh.get_all_edgesInBox(),
                             linewidths=0.1, colors='black')
    ax.add_collection(ln_coll, autolim=True)

    # impose the dimensions of the box as the limits of the figure
    # ax.set_xlim(myMesh.get_xlim())
    # ax.set_ylim(myMesh.get_zlim())

    # to avoid distorting the mesh:
    ax.set_aspect('equal')

    # to don't print axis:
    # ax.axis('off')

    # plot text to indicate time
    xtext = xmin + 0.05 * (xmax - xmin)
    if case == '2D':
        ztext = zmax - 0.2 * (zmax - zmin)
    else:
        ztext = ymax - 0.2 * (ymax - ymin)

    timeLegend = str(float(time) - 10)
    ax.text(
        xtext, ztext, f"time = {timeLegend} s",
        fontsize=15,
        bbox=dict(
            boxstyle="round",
            ec="gray",
            fc="white"
        )
    )

    # ax.tick_params(labelsize=12)
    if case == '2D':
        ax.set_ylabel(r"$z\,[m]$")

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(zmin, zmax)
    else:
        ax.set_ylabel(r"$z\,[m]$")

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(zmin, zmax)

    if i < ntimes - 1:
        ax.tick_params(bottom=False, labelbottom=False)
    else:
        ax.set_xlabel(r"$x\,[m]$")

fig.tight_layout()
plt.show()

if saveFig:
    fig.savefig("Figures/" + figName + ".eps", dpi=100, format="eps")
    fig.savefig("Figures/" + figName + ".png", dpi=300, format="png")
    print("figure saved")
