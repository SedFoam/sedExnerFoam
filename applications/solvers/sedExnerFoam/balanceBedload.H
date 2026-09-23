/*---------------------------------------------------------------------------* \
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


// Access to shear stress induced bedload and avalanche related bedload
areaVectorField& qb = qbPtr.ref();
areaVectorField& qav = qavPtr.ref();

// project qb/qav on horizontal plane
vectorField qbProj = meshProj.project(qb);
vectorField qavProj = meshProj.project(qav);

// projected bedload flux through projected faMesh edges
scalarField& phiqbProj = phiqbProjPtr.ref();

// projected edge centers coordinates
const vectorField& xEdgesProj = meshProj.edgeCentres();
// projected face centers coordinates
const vectorField& xFacesProj = meshProj.areaCentres();
// projected edge vectors, oriented toward owner face
const vectorField& LeProj = meshProj.Le();
// face owner list
const labelList& owner = meshProj.edgeOwner();
// face neighbour list
const labelList& neighbour = meshProj.edgeNeighbour();

if (qbScheme == "centered")
{
    // compute bedload flux on through internal edges
    for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
    {
        // distance between edgei center and owner face center
        scalar distOwner =
            Foam::mag(xFacesProj[owner[edgei]] - xEdgesProj[edgei]);
        // distance between edgei center and neighbour center
        scalar distNeighbour =
            Foam::mag(xFacesProj[neighbour[edgei]] - xEdgesProj[edgei]);
        // linear interpolate qb on edge i
        vector qbEdgei =
            (
                distNeighbour * qbProj[owner[edgei]]
                + distOwner * qbProj[neighbour[edgei]]
            )
            * (1 / (distOwner + distNeighbour));

        //
        vector qavEdgei = qavProj[owner[edgei]];
        if (qavScheme == "centered")
        {
            // linear interpolation for avalanche bedload
            qavEdgei =
                (
                    distNeighbour * qavProj[owner[edgei]]
                    + distOwner * qavProj[neighbour[edgei]]
                )
                * (1 / (distOwner + distNeighbour));
        }
        else if (qavScheme == "upwind")
        {
            // upwind interpolation for avalanche bedload
            if ((qavEdgei & LeProj[edgei]) >= scalar(0.0))
            {
                qavEdgei = qavProj[owner[edgei]];
            }
            else
            {
                qavEdgei = qavProj[neighbour[edgei]];
            }
        }
        else
        {
            FatalError << qavScheme
                << " is not a valid entry for qavScheme"
                << endl;
            Info << abort(FatalError) << endl;
        }
        phiqbProj[edgei] += (qbEdgei+qavEdgei) & LeProj[edgei];
    }
}
else if (qbScheme == "upwind")
{
    // compute bedload flux on through internal edges
    for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
    {
        vector qbEdgei = qbProj[owner[edgei]];
        // upwind interpolation for qb
        if ((qbEdgei & LeProj[edgei]) >= scalar(0.0))
        {
            qbEdgei = qbProj[owner[edgei]];
        }
        else
        {
            qbEdgei = qbProj[neighbour[edgei]];
        }
        // Avalanche flux
        vector qavEdgei = qavProj[owner[edgei]];
        if (qavScheme == "centered")
        {
            // distance between edgei center and owner face center
            scalar distOwner =
                Foam::mag(xFacesProj[owner[edgei]] - xEdgesProj[edgei]);
            // distance between edgei center and neighbour center
            scalar distNeighbour =
                Foam::mag(xFacesProj[neighbour[edgei]] - xEdgesProj[edgei]);
            // linear interpolation for avalanche bedload
            qavEdgei =
                (
                    distNeighbour * qavProj[owner[edgei]]
                    + distOwner * qavProj[neighbour[edgei]]
                )
                * (1 / (distOwner + distNeighbour));
        }
        else if (qavScheme == "upwind")
        {
            // upwind interpolation for avalanche bedload
            if ((qavEdgei & LeProj[edgei]) >= scalar(0.0))
            {
                qavEdgei = qavProj[owner[edgei]];
            }
            else
            {
                qavEdgei = qavProj[neighbour[edgei]];
            }
        }
        else
        {
            FatalError << qavScheme
                << " is not a valid entry for qavScheme"
                << endl;
            Info << abort(FatalError) << endl;
        }
        phiqbProj[edgei] += (qbEdgei+qavEdgei) & LeProj[edgei];
    }
}
else
{
    FatalError << qbScheme
        << " is not a valid entry for qbScheme"
        << endl;
    Info << abort(FatalError) << endl;
}

// /* test JC
// compute bedload flux on through boundary edges
forAll(aMesh.boundary(), patchi)
{
    const faPatchField<Vector<double>>& QbBound =
        qb.boundaryField()[patchi];

    if (aMesh.boundary()[patchi].type()!="empty")
    {
        // faMesh edge adressing
        label edgei = aMesh.boundary()[patchi].start();

        forAll(aMesh.boundary()[patchi], patchEdgei)
        {
            const vector edgeLeProj =
                meshProj.project
                (
                    aMesh.Le().boundaryField()[patchi][patchEdgei]
                );

            // linear interpolation from face to edge
            vector qbEdgei = meshProj.project(QbBound[patchEdgei]);

            phiqbProj[edgei] += qbEdgei & edgeLeProj;

            edgei++;
        }
    }
}
//*/

// cancel computed bedload flux where no more sediments are present
if (bed.rigidBed())
{
    //const scalarField& SfProj = meshProj.S();

    areaScalarField& HsedBed = HsedBedPtr.ref();

    //scalar deltaT = runTime.time().deltaT().value();

    for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
    {
        if (phiqbProj[edgei]>0)
        {
            /*scalar phiqbMax =
                HsedBed[owner[edgei]] * SfProj[owner[edgei]] / deltaT;
                phiqbProj[edgei] = min(phiqbProj[edgei], phiqbMax);*/
            if (HsedBed[owner[edgei]]<=1e-4)
            {
                phiqbProj[edgei] = 0;
            }
        }
        else
        {
            /*scalar phiqbMin =
                - HsedBed[neighbour[edgei]]
                * SfProj[neighbour[edgei]] / deltaT;
                phiqbProj[edgei] = max(phiqbProj[edgei], phiqbMin);*/
            if (HsedBed[neighbour[edgei]]<=1e-4)
            {
                phiqbProj[edgei] = 0;
            }
        }
    }
    // rigid bed treatment on boundaries
    for (
    label edgei = meshProj.nInternalEdges();
    edgei < meshProj.nEdges();
    edgei++)
    {
        if (phiqbProj[edgei]>0)
        {
            if (HsedBed[owner[edgei]]<=0)
            {
                phiqbProj[edgei] = 0;
            }
        }
    }
}

for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
{
    dVsed[owner[edgei]] -= phiqbProj[edgei];
    dVsed[neighbour[edgei]] += phiqbProj[edgei];
}

// apply boundary condition to exner equation
for (
    label edgei = meshProj.nInternalEdges();
    edgei < meshProj.nEdges();
    edgei++
    )
{
    dVsed[owner[edgei]] -= phiqbProj[edgei];
}
