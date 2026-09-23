/*---------------------------------------------------------------------------*\
Copyright (C) 2022 Matthias Renaud, Cyrille Bonamy, Julien Chauchat
                   and contributors

License
    This file is part of ScourFOAM.

    ScourFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ScourFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with ScourFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

const faMesh& aMesh = bed.aMesh();

areaScalarField& dH = dHptr.ref();

// access to bedload field
areaVectorField& qb = qbPtr.ref();
areaVectorField& qav = qavPtr.ref();

// access to bedload flux through edges
edgeScalarField& phiqb = phiqbPtr.ref();
edgeScalarField& phiqav = phiqavPtr.ref();

phiqb = (fac::interpolate(qb) & aMesh.Le());
phiqav = (fac::interpolate(qav) & aMesh.Le());

scalar deltaT = runTime.time().deltaT().value();  // current time step

Info << "create exnerEqn matrix" << endl;
faScalarMatrix exnerEqn
(
    fam::ddt(dH)
  ==
  - fac::div(qb)
  - fac::div(qav)
);

exnerEqn.solve();

if (bed.bedMotion())
{
    #include "meshMove.H"
}
