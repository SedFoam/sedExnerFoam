"""
Class to generate 1D dune case to run with exnerFoam solver

duneCase1D

zb0_gaussian: initial gaussian condition for dune
"""

from fluidfoam import readof as rdf
import numpy as np
import os
import shutil
import subprocess
from pathlib import Path


def zb0_gaussian(x, hDune=0.1, dunePos=2., duneW=0.6):
    """
    Initial gaussian dune
    - hDune: float, initial dune height
    - dunePos: float, initial position of dune crest
    - duneW: float, parameter linked to dune width
    """
    return hDune * np.exp(-((x-dunePos)/duneW)**2)


def bedformCelerity(duneCase, zb):
    """
    Return bed form wave celerity
    """
    alpha, beta = duneCase.bedloadParams()
    Q, H = duneCase.hydroParams()
    return alpha * beta * Q**beta / (H - zb)**(beta+1)


def require_not_run(method):

    def wrapper(self, *args, **kwargs):

        if self.has_run:
            raise RuntimeError(
                f"{method.__name__} cannot be called"
                + " after simulation execution"
            )

        return method(self, *args, **kwargs)

    return wrapper


def require_exist(method):

    def wrapper(self, *args, **kwargs):

        if not self.exist:
            raise FileNotFoundError(
                f"{method.__name__} cannot be called"
                + f" as path {self._path} does not exist"
            )

        return method(self, *args, **kwargs)

    return wrapper


class DuneCase1D:

    def __init__(self, pathCase0, pathCase, dt, nx, Lx, Hz, Qwater, bedloadDict):
        """
        - pathCase0: str
        - pathCase: str
        - dt: float, time step value
        - nx: int, number of faces making up the mesh
        - Lx: float, domain length
        - Hz: float, water depth
        - Qwater: float, water discharge per unit width
        """
        self._path0 = Path(pathCase0)
        self._path = Path(pathCase)
        self._dt = dt
        self._nx = nx
        if self._nx <= 0:
            raise ValueError("nx must be > 0")
        self._Lx = Lx
        if self._Lx <= 0:
            raise ValueError("domain length Lx must be positive")
        self._Hz = Hz
        if self._Hz <= 0:
            raise ValueError("domain height Hz must be positive")
        self._Qwater = Qwater
        # bedloadProperties parameters
        self._alphaQb = bedloadDict.get("alpha", 0.05)
        self._betaQb = bedloadDict.get("beta", 1.5)
        self._qbScheme = bedloadDict.get("qbScheme", "linear")
        self._ABorder = bedloadDict.get("ABorder", 1)
        self._create()

    @property
    def exist(self):
        """
        False if the simulation case does not exist
        """
        return (
            self._path.is_dir()
            and
            (self._path / "system" / "controlDict").exists()
        )

    @property
    def is_initialized(self):
        """
        False if the initial condition has not been set
        """
        return (self._path / "0" / "finite-area" / "qb").exists()

    @property
    def is_setup(self):
        """
        False if case as not been set up
        """
        return (self._path / "constant" / "polyMesh").exists()

    @property
    def has_run(self) -> bool:
        """
        False if case has not been run
        """
        return len(self.getTimes()) > 1

    def _changeInitialMesh(self, getZb):
        """
        - getZb: function, take x coordinates as input and return
        """
        if not self.is_setup:
            raise ValueError(
                "Case must be set up before initializing the mesh"
                + ", run setUp method")
        pathPoints = self._path / "constant" / "polyMesh" / "points"
        pathNewPoints = self._path / "constant" / "polyMesh" / "newPoints"
        with open(pathNewPoints, "w") as npF:
            with open(pathPoints, "r") as opF:
                for line in opF:
                    newline = line
                    if line[0] == "(" and line[-2:] == ")\n":
                        xs, ys, zs = line[:-1].strip("()").split(" ")
                        x, z = float(xs), float(zs)
                        zb = getZb(x)
                        z = zb + z * (1 - zb / self._Hz)
                        zs = str(round(z, 10))
                        newline = f"({xs} {ys} {zs})\n"
                    npF.write(newline)
        os.system(f"rm {pathPoints} && mv {pathNewPoints} {pathPoints}")

    def _create(self):
        shutil.copytree(self._path0, self._path)

    @require_exist
    def _editBedloadProperties(self):
        # if self.has_run():
        pathBedloadProperties = self._path / "constant" / "bedloadProperties"
        subprocess.run(
            [
                "sed", "-i", f"s/insert_alpha/{self._alphaQb}/",
                pathBedloadProperties
            ], check=True
        )
        subprocess.run(
            [
                "sed", "-i", f"s/insert_beta/{self._betaQb}/",
                pathBedloadProperties
            ], check=True
        )
        subprocess.run(
            [
                "sed", "-i", f"s/insert_qbScheme/{self._qbScheme}/",
                pathBedloadProperties
            ],
            check=True
        )
        subprocess.run(
            [
                "sed", "-i", f"s/insert_ABorder/{self._ABorder}/",
                pathBedloadProperties
            ],
            check=True
        )

    @require_exist
    def _editMesh(self):
        pathBlockMeshDict = self._path / "system" / "blockMeshDict"
        subprocess.run(
            ["sed", "-i", f"s/insert_Lx/{self._Lx}/", pathBlockMeshDict],
            check=True
        )
        subprocess.run(
            ["sed", "-i", f"s/insert_Hz/{self._Hz}/", pathBlockMeshDict],
            check=True
        )
        subprocess.run(
            ["sed", "-i", f"s/insert_nx/{self._nx}/", pathBlockMeshDict],
            check=True
        )

    @require_exist
    def _editControl(self):
        pathControlDict = self._path / "system" / "controlDict"
        subprocess.run(
            ["sed", "-i", f"s/insert_dt/{self._dt}/", pathControlDict],
            check=True
        )

    @require_exist
    def delete(self):
        shutil.rmtree(self._path)

    def setup(self):
        """
        Edit mesh, bedload properties and controlDict
        Create volumic mesh and finite area mesh.
        At this point, bed elevation has not been initialized.
        """
        self._editMesh()
        self._editControl()
        self._editBedloadProperties()
        print("Create volumic mesh with blockMesh utility")
        with open(self._path / "log.blockMesh", "w") as log:
            subprocess.run(
                ["blockMesh"], cwd=self._path, stdout=log,
                stderr=subprocess.STDOUT,  check=True
            )
        print("Create finite-area mesh with makeFaMesh utility")
        with open(self._path / "log.makeFaMesh", "w") as log:
            subprocess.run(
                ["makeFaMesh"], cwd=self._path, stdout=log,
                stderr=subprocess.STDOUT,  check=True
            )

    def setInitialCondition(self, hDune, dunePos, duneWidth):
        """
        Copy the 0_org folder into a 0 file (boundary conditions)
        Set the initial bed elevation
        - hDune: float, initial dune height
        - dunePos: float, initial position of dune crest
        - duneWidth: float, parameter linked to dune width
        """
        shutil.copytree(self._path / "0_org", self._path / "0")

        def getZb(x):
            return zb0_gaussian(x, hDune, dunePos, duneWidth)

        self._changeInitialMesh(getZb)

    @require_not_run
    def runSim(self):
        if not self.is_initialized:
            raise ValueError("The initial condition has not been set")
        print(f"lLaunch dune transport case: {self._path}")
        with open(self._path / "log.exnerFoam", "w") as log:
            subprocess.run(
                ["projExnerFoam"], cwd=self._path, stdout=log,
                stderr=subprocess.STDOUT,  check=True
            )

    def getFaceCenters(self, time, verbose=False):
        print(f"read face centers coordinates at time {time}")
        return rdf.readmesh(
            str(self._path), time_name=time, boundary="bed", verbose=verbose)

    def getField(self, time, fieldName, verbose=False):
        """
        compatible avec les champs scalaires et vectoriels
        """
        print(f"read field {fieldName} at time {time}")
        try:
            field = rdf.readscalar(
                str(self._path), time, fieldName,
                boundary="bed", verbose=verbose
            )
        except ValueError:
            field = rdf.readvector(
                str(self._path), time, fieldName,
                boundary="bed", verbose=verbose
            )
        return field

    def getTimes(self):
        result = subprocess.run(
            ["foamListTimes", "-withZero",  "-case", self._path],
            capture_output=True,
            text=True,
            check=True
        )
        timeList = result.stdout.splitlines()
        return timeList

    def bedloadParams(self) -> tuple[float, float]:
        """
        Return alpha and beta coefficients for bedload transport
        """
        return self._alphaQb, self._betaQb

    def hydroParams(self) -> tuple[float, float]:
        """
        Return Q and H parameters, discharge per unit with and water depth
        """
        return self._Qwater, self._Hz
