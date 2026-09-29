"""
colorplot and streamplot of velocity fields with the dune at different times
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import griddata
from fluidfoam import readof as rdf

plt.rcParams["font.size"] = 15

pathCase = "../RAS/duneMigration"

saveFig = True
figName = "Ufields_dune2D"

timeList = ["10", "12", "16"]  # at least two times must be given
ntimes = len(timeList)

xmin, xmax = -0.1, 0.4
zmin, zmax = 0, 0.09

vmin, vmax = 0, 0.7  # limit values for colormap

# Number of division for linear interpolation
ngridx = 400
ngridz = 300

# Interpolation grid dimensions
xinterpmin = xmin
xinterpmax = xmax
zinterpmin = zmin
zinterpmax = zmax

# Interpolation grid
xi = np.linspace(xinterpmin, xinterpmax, ngridx)
yi = np.linspace(zinterpmin, zinterpmax, ngridz)

# Structured grid creation
xinterp, zinterp = np.meshgrid(xi, yi)


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        "trunc({n},{a:.2f},{b:.2f})".format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


cmap = plt.get_cmap("Blues")
myCmap = truncate_colormap(cmap, 0, 0.85)

fig = plt.figure(figsize=(12, 4 * ntimes))
gs = fig.add_gridspec(
    ntimes + 1, 1, height_ratios=[0.1] + [1 for i in range(ntimes)])

# ax for colorbar
cax = fig.add_subplot(gs[0, 0])
# make colorbar ax occupy less space on figure
l, b, w, h = cax.get_position().bounds
cax.set_position(pos=[l + 0.1 * w, b, 0.8 * w, h * 0.8])

for i, time in enumerate(timeList):
    ax = fig.add_subplot(gs[i + 1, 0])
    print(f"\n- time: {time} s")
    # readmesh, cell centers and bed faces
    Xmesh, Ymesh, Zmesh = rdf.readmesh(pathCase, time_name=time)
    Xbed, Ybed, Zbed = rdf.readmesh(pathCase, boundary="bed", time_name=time)
    # read velocity field
    Ux, Uy, Uz = rdf.readvector(pathCase, time, "U")
    magU = np.sqrt(Ux**2 + Uz**2)

    # plot stream lines
    ux_i = griddata((Xmesh, Zmesh), Ux, (xinterp, zinterp), method="linear")
    uz_i = griddata((Xmesh, Zmesh), Uz, (xinterp, zinterp), method="linear")
    magU_i = np.sqrt(ux_i**2 + uz_i**2)
    lwStream = 1 * magU_i**(1 / 4)
    im = ax.tripcolor(
        Xmesh * 1000, Zmesh * 1000, magU, vmin=vmin, vmax=vmax,
        cmap=myCmap, shading="gouraud")
    if i == 0:
        cbar = fig.colorbar(
            im, cax=cax,
            orientation="horizontal", location="top",
            ticklocation="bottom", extend="max")
        cbar.set_label(r"$u\,[m.s^{-1}]$")
        cbar.ax.xaxis.set_label_position("top")
        # cbar.ax.tick_params(labelsize=12)
    ax.streamplot(
        xi * 1000, yi * 1000, ux_i, uz_i, color="#000000", density=[2, 1],
        linewidth=lwStream, arrowsize=0.05)

    Zbotbed = np.zeros_like(Xbed)

    ax.fill_between(Xbed * 1000, Zbed * 1000, Zbotbed, color="peru", zorder=4.)
    ax.plot(Xbed * 1000, Zbed * 1000, color="black", lw=3.)

    # plot text to indicate time
    xtext = xmin + 0.05 * (xmax - xmin)
    ztext = zmax - 0.2 * (zmax - zmin)
    timeLegend = str(float(time) - 10)
    ax.text(
        1000 * xtext, ztext * 1000, f"time = {timeLegend} s",
        fontsize=15,
        bbox=dict(
            boxstyle="round",
            ec="gray",
            fc="white"
        )
    )

    # ax.tick_params(labelsize=12)

    ax.set_ylabel(r"$z\,[mm]$")
    # cbarU.set_label(r"U [$m.s^{-1}$]")
    ax.set_xlim(1000 * xmin, 1000 * xmax)
    ax.set_ylim(1000 * zmin, 1000 * zmax)

    if i < ntimes - 1:
        ax.tick_params(bottom=False, labelbottom=False)
    else:
        ax.set_xlabel(r"$x\,[mm]$")

plt.show()

if saveFig:
    fig.savefig("Figures/" + figName + ".png", dpi=300, format="png")
    fig.savefig("Figures/" + figName + ".eps", dpi=300, format="eps")
    print("figure saved")
