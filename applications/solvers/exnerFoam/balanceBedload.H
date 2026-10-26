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

// project qb/qav on horizontal plane
vectorField qbProj = meshProj.project(qb);

// Initialize to zero
phiqbProj = Zero;

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


// compute gradient of qb at face centers
// necessary for linearUpwind scheme
if (qbScheme == "linearUpwind")
{
    // Re-initialize projected gradient of qb at face centers
    gradQbProj = Zero;

    for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
    {
        // linear interpolation of qb from faces to edges
        scalar distOwner =
            Foam::mag(xFacesProj[owner[edgei]] - xEdgesProj[edgei]);
        scalar distNeighbour =
            Foam::mag(xFacesProj[neighbour[edgei]] - xEdgesProj[edgei]);
        vector qbEdgei = (
            distNeighbour * qbProj[owner[edgei]] + distOwner *
            qbProj[neighbour[edgei]]
        ) / (distOwner + distNeighbour);

        // compute gradient of qb
        tensor phiGradQbEdgei = qbEdgei * LeProj[edgei];
        gradQbProj[owner[edgei]] += phiGradQbEdgei / SfProj[owner[edgei]];
        gradQbProj[neighbour[edgei]] -=
            phiGradQbEdgei / SfProj[neighbour[edgei]];
    }
    forAll(aMesh.boundary(), patchi)
    {
        const faPatchField<vector>& QbBound = qb.boundaryField()[patchi];

        if (aMesh.boundary()[patchi].type() != "empty")
        {
            label edgei = aMesh.boundary()[patchi].start();
            forAll(aMesh.boundary()[patchi], patchEdgei)
            {
                const vector edgeLeProj =
                    meshProj.project
                    (
                        aMesh.Le().boundaryField()[patchi][patchEdgei]
                    );

                vector qbEdgei = meshProj.project(QbBound[patchEdgei]);
                tensor phiGradQbEdgei = qbEdgei * edgeLeProj;
                gradQbProj[owner[edgei]] +=
                    phiGradQbEdgei / SfProj[owner[edgei]];
                edgei++;
            }
        }
    }
}

// === INTERNAL EDGES ===
if (qbScheme == "linear")
{
    for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
    {
        scalar distOwner =
            Foam::mag(xFacesProj[owner[edgei]] - xEdgesProj[edgei]);
        scalar distNeighbour =
            Foam::mag(xFacesProj[neighbour[edgei]] - xEdgesProj[edgei]);
        vector qbEdgei = (
            distNeighbour * qbProj[owner[edgei]] + distOwner *
            qbProj[neighbour[edgei]]
        ) / (distOwner + distNeighbour);

        phiqbProj[edgei] = qbEdgei & LeProj[edgei];
    }
}
else if (qbScheme == "upwind")
{
    for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
    {
        // linear interpolation to get directionality of upwind scheme
        scalar distOwner =
            Foam::mag(xFacesProj[owner[edgei]] - xEdgesProj[edgei]);
        scalar distNeighbour =
            Foam::mag(xFacesProj[neighbour[edgei]] - xEdgesProj[edgei]);
        vector qbEdgei = (
            distNeighbour * qbProj[owner[edgei]] + distOwner *
            qbProj[neighbour[edgei]]
        ) / (distOwner + distNeighbour);

        // upwind scheme
        if ((qbEdgei & LeProj[edgei]) >= scalar(0.0))
        {
            qbEdgei = qbProj[owner[edgei]];
        }
        else
        {
            qbEdgei = qbProj[neighbour[edgei]];
        }

        phiqbProj[edgei] = qbEdgei & LeProj[edgei];
    }
}
else if (qbScheme == "linearUpwind")
{
    // linear interpolation ofProj from face to edges
    for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
    {
        // linear interpolation to get directionality of upwind scheme
        scalar distOwner =
            Foam::mag(xFacesProj[owner[edgei]] - xEdgesProj[edgei]);
        scalar distNeighbour =
            Foam::mag(xFacesProj[neighbour[edgei]] - xEdgesProj[edgei]);
        vector qbEdgei = (
            distNeighbour * qbProj[owner[edgei]] + distOwner *
            qbProj[neighbour[edgei]]
        ) / (distOwner + distNeighbour);

        // upwind scheme
        if ((qbEdgei & LeProj[edgei]) >= scalar(0.0))
        {
            vector dpUp = xEdgesProj[edgei] - xFacesProj[owner[edgei]];
            qbEdgei = qbProj[owner[edgei]];
            qbEdgei += gradQbProj[owner[edgei]] & dpUp;
        }
        else
        {
            vector dpUp = xEdgesProj[edgei] - xFacesProj[neighbour[edgei]];
            qbEdgei = qbProj[neighbour[edgei]];
            qbEdgei += gradQbProj[neighbour[edgei]] & dpUp;
        }

        phiqbProj[edgei] = qbEdgei & LeProj[edgei];
    }
}
else
{
    FatalError
        << "Invalid option " << qbScheme << " for qbSchem" << endl
        << "valid options are: linear upwind linearUpwind" << endl;
    Info << abort(FatalError) << endl;
}

forAll(aMesh.boundary(), patchi)
{
    const faPatchField<vector>& QbBound = qb.boundaryField()[patchi];

    if (aMesh.boundary()[patchi].type() != "empty")
    {
        label edgei = aMesh.boundary()[patchi].start();
        forAll(aMesh.boundary()[patchi], patchEdgei)
        {
            const vector edgeLeProj =
                meshProj.project
                (
                    aMesh.Le().boundaryField()[patchi][patchEdgei]
                );

            vector qbEdgei = meshProj.project(QbBound[patchEdgei]);
            phiqbProj[edgei] = qbEdgei & edgeLeProj;
            edgei++;
        }
    }
}

// Internal edges: flux from owner to neighbour
for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
{
    dVsed[owner[edgei]] -= phiqbProj[edgei];
    dVsed[neighbour[edgei]] += phiqbProj[edgei];
}

// Boundary edges: flux out of domain (owner is always the internal face)
for (
    label edgei = meshProj.nInternalEdges();
    edgei < meshProj.nEdges();
    edgei++
)
{
    dVsed[owner[edgei]] -= phiqbProj[edgei];
}
