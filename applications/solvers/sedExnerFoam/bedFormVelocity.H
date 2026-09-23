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

vectorField& VelBfE = VelBfEdgesPtr.ref();
areaVectorField& VelBf = VelBfPtr.ref();

// face centers elevation, zb on faces
scalarField zbFaces = -bed.aMesh().areaCentres() & eg;

scalarField normWeight(bed.aMesh().nFaces());
normWeight = Zero;

// reset bedform velocities to zero
VelBf.primitiveFieldRef() = Zero;

for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
{
    // compute bedload flux on internal edges
    VelBfE[edgei] = qbProj[owner[edgei]] - qbProj[neighbour[edgei]];
    VelBfE[edgei] += qavProj[owner[edgei]] - qavProj[neighbour[edgei]];
    scalar dzb = zbFaces[owner[edgei]] - zbFaces[neighbour[edgei]];
    VelBfE[edgei] *= 1 / (dzb + SMALL);

    // interpolate VelBf from edges to face centers, inverse distance
    // distance between edgei center and owner face center
    scalar distOwner =
        Foam::mag(xFacesProj[owner[edgei]] - xEdgesProj[edgei]);
    // distance between edgei center and neighbour center
    scalar distNeighbour =
        Foam::mag(xFacesProj[neighbour[edgei]] - xEdgesProj[edgei]);

    VelBf[owner[edgei]] += (1/distOwner) * VelBfE[edgei];
    normWeight[owner[edgei]] += (1/distOwner);

    VelBf[neighbour[edgei]] += (1/distNeighbour) * VelBfE[edgei];
    normWeight[neighbour[edgei]] += (1/distNeighbour);
}

// normalize after interpolation
VelBf.primitiveFieldRef() /= normWeight;
