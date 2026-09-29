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

// access to projected faMesh
const projectedFaMesh& meshProj = bed.aProjMesh();

// projected edge centers coordinates
const vectorField& xEdgesProj = meshProj.edgeCentres();
// projected face centers coordinates
const vectorField& xFacesProj = meshProj.areaCentres();
// face owner list
const labelList& owner = meshProj.edgeOwner();
// face neighbour list
const labelList& neighbour = meshProj.edgeNeighbour();


// filtered Shields values
vectorField dsh(bed.aMesh().nFaces());

// normalise interpolation weight
scalarField normWeight(bed.aMesh().nFaces());
normWeight = Zero;

// loop on nternal edges
for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
{
    // distance between edgei center and owner face center
    scalar distOwner =
        Foam::mag(xFacesProj[owner[edgei]] - xEdgesProj[edgei]);
    // distance between edgei center and neighbour center
    scalar distNeighbour =
        Foam::mag(xFacesProj[neighbour[edgei]] - xEdgesProj[edgei]);

    // Shields value on edge, linear interpolation from faces
    vector shEi =
        (
            distOwner * ssp[neighbour[edgei]]
            + distNeighbour * ssp[owner[edgei]]
        ) / (distOwner + distNeighbour);

    dsh[owner[edgei]] += (1/distOwner) * shEi;
    normWeight[owner[edgei]] += (1/distOwner);
    
    dsh[neighbour[edgei]] += (1/distNeighbour) * shEi;
    normWeight[neighbour[edgei]] += (1/distNeighbour);
}

dsh /= normWeight;

ssp = dsh;
