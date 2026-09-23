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

// access bedloadModel
const bedloadModels::bedloadModel& bedloadMod = bed.bedloadModel();

// access to avalanche bedload
areaVectorField& qav = qavPtr.ref();

if (not onEdgeAvalanche)
{
    // compute bed slope angle
    const vectorField& nFaces =
        bed.aMesh().faceAreaNormals().internalField();

    // slope direction, align with bed
    vectorField slopeDir = eg - nFaces * (eg & nFaces);
    slopeDir /= (Foam::mag(slopeDir) + SMALL);

    scalarField& betaIntern = betaPtr.ref().primitiveFieldRef();

    // Vinent law (2019), additional bedload corresponding to avalanche
    vectorField& qavIntern = qav.primitiveFieldRef();
    qavIntern = bedloadMod.qbAvalanche(betaIntern, slopeDir, betaRep);

    //qb.primitiveFieldRef() += qavIntern;
    Info<<"Max diffusion due to avalanche:"
        <<gMax(mag(qavIntern)/(betaIntern+SMALL))<<endl;

    qav.correctBoundaryConditions();
}
else
{
    // projected bedload flux through projected faMesh edges
    scalarField& phiqbProj = phiqbProjPtr.ref();

    // faMesh edge normals
    const vectorField& nEdges = bed.aMesh().edgeAreaNormals();

    // bed edges angle
    scalarField slopeAngle = Foam::acos(nEdges & eg);

    // slope direction
    vectorField slopeDir = eg - nEdges * (eg & nEdges);
    slopeDir /= (Foam::mag(slopeDir) + SMALL);

    // edge length vectors
    const vectorField& Le = bed.aMesh().Le();

    vectorField qbAvEdge =
        bedloadMod.qbAvalanche(slopeAngle, slopeDir, betaRep);

    for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
    {
        phiqbProj[edgei] += qbAvEdge[edgei] & Le[edgei];
    }
}
