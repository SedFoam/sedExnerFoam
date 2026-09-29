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

// reference to erosion and deposition terms
volScalarField& Ero = erosionPtr.ref();
volScalarField& Dep = depositionPtr.ref();

volScalarField& Cs = CsPtr.ref();  // suspended load
volVectorField& Ws = WsPtr.ref();  // settling velocity

areaScalarField& DE = DEptr.ref();  // deposition - erosion
surfaceScalarField& bedPhi = bedPhiPtr.ref();  // mesh flux through faces

//areaScalarField& dH = dHptr.ref();
//scalar deltaT = runTime.time().deltaT().value();


const faMesh& aMesh = bed.aMesh();

for (const label patchID : bed.bedPatchesID())
{
    // deposition boundary field
    scalarField& DepPatch = Dep.boundaryFieldRef()[patchID];
    // erosion boundary field
    scalarField& EroPatch = Ero.boundaryFieldRef()[patchID];
    // reference to bed fvPatch
    const fvPatch& bedPatch = mesh.boundary()[patchID];
    const scalarField& bedPhiPatch = bedPhi.boundaryFieldRef()[patchID];
    
    //const scalarField& diffSedPatch = diffSed.boundaryFieldRef()[patchID];
    // get gradient of Cs on bed patch faces

    forAll(bedPatch, facei)
    {
        // label of owner cell of face
        label idNCell = bedPatch.faceCells()[facei];
        //scalar ci = Cs[idNCell];
        scalar ci = Cs.boundaryField()[patchID][facei];
        vector wsi = Ws[idNCell];
        scalar Sface = aMesh.S()[facei];
        scalar Vcell = mesh.V()[idNCell];
        scalar bedPhiFace = bedPhiPatch[facei];
        const vector& nFace = aMesh.faceAreaNormals()[facei];
        // deposition due to settling
        // option 1
        DepPatch[facei] = ci * ((wsi & nFace) - bedPhiFace/Sface);
        // option 2
        //DepPatch[facei] = ci * (Foam::mag(wsi) + dH[facei]/deltaT);
        // deposition due to diffusion, not implemented
        // error in case of non zero diffusion on bed
        DE[facei] = DepPatch[facei] - EroPatch[facei] * (Vcell/Sface);
    }
}

// deposition and erosion contribution in Exner equation
forAll(bed.aMesh().areaCentres(), facei)
{
    dVsed[facei] += DE[facei] * aMesh.S()[facei];  // option 1
    //dVsed[facei] += DE[facei] * SfProj[facei];  // option 2
}
