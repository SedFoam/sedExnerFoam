/*---------------------------------------------------------------------------* \
Copyright (C) 2026 Arthur Lopez,Matthias Renaud, Cyrille Bonamy, Julien Chauchat
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

#include "coupledFaPatch.H"

inline label projectedGlobalEdgeCount(const projectedFaMesh& meshProj)
{
    label localCount = meshProj.nInternalEdges();

    forAll(meshProj.mesh().boundary(), patchi)
    {
        const faPatch& patch = meshProj.mesh().boundary()[patchi];
        if (patch.coupled())
        {
            localCount += patch.size();
        }
    }

    return returnReduce(localCount, sumOp<label>());
}


inline scalar projectedMinEdgeDistance(const projectedFaMesh& meshProj)
{
    const vectorField& xFacesProj = meshProj.areaCentres();
    const labelList& owner = meshProj.edgeOwner();
    const labelList& neighbour = meshProj.edgeNeighbour();

    scalar localMin = VGREAT;

    for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
    {
        const scalar distance = Foam::mag
        (
            xFacesProj[neighbour[edgei]] - xFacesProj[owner[edgei]]
        );

        if (distance > VSMALL)
        {
            localMin = min(localMin, distance);
        }
    }

    forAll(meshProj.mesh().boundary(), patchi)
    {
        const faPatch& patch = meshProj.mesh().boundary()[patchi];
        if (not patch.coupled())
        {
            continue;
        }

        const coupledFaPatch& coupledPatch =
            refCast<const coupledFaPatch>(patch);
        const tmp<vectorField> tDelta = coupledPatch.delta();
        const vectorField& delta = tDelta();

        forAll(patch, patchEdgei)
        {
            const scalar distance = Foam::mag
            (
                meshProj.project(delta[patchEdgei])
            );

            if (distance > VSMALL)
            {
                localMin = min(localMin, distance);
            }
        }
    }

    const scalar globalMin = returnReduce(localMin, minOp<scalar>());
    return globalMin == VGREAT ? scalar(0) : globalMin;
}

template<class Field, class ExchangeField>
inline Field filterField(
    const label N,
    const scalar Kex,
    const Field& value,
    const projectedFaMesh& meshProj,
    const scalarField& SfProj,
    const vectorField& xFacesProj,
    const vectorField& LeProj,
    const labelList& owner,
    const labelList& neighbour,
    ExchangeField& exchangeField)
{
    typedef typename Field::value_type valueType;

    Field returnValue(value);
    Field exchangeOriginal(exchangeField.primitiveField());

    for (label j = 0; j < N; j++)
    {
        Field laplValue(value);
        laplValue = Zero;

        exchangeField.primitiveFieldRef() = returnValue;
        exchangeField.correctBoundaryConditions();

        for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
        {
            const vector Xon =
                xFacesProj[neighbour[edgei]] - xFacesProj[owner[edgei]];
            const scalar XonMagSqr = Foam::magSqr(Xon);

            if (XonMagSqr <= VSMALL)
            {
                continue;
            }

            const valueType graddVe =
                (((returnValue[neighbour[edgei]] - returnValue[owner[edgei]])
                  * (Xon & LeProj[edgei]))
                 / XonMagSqr);

            laplValue[owner[edgei]] +=
                (Kex / SfProj[owner[edgei]]) * graddVe;
            laplValue[neighbour[edgei]] -=
                (Kex / SfProj[neighbour[edgei]]) * graddVe;
        }

        forAll(meshProj.mesh().boundary(), patchi)
        {
            const faPatch& patch = meshProj.mesh().boundary()[patchi];
            if (not patch.coupled())
            {
                continue;
            }

            const coupledFaPatch& coupledPatch =
                refCast<const coupledFaPatch>(patch);
            const tmp<vectorField> tDelta = coupledPatch.delta();
            const vectorField& delta = tDelta();
            const auto& valueBound = exchangeField.boundaryField()[patchi];

            label edgei = patch.start();
            forAll(patch, patchEdgei)
            {
                const label own = owner[edgei];
                const vector Xon = meshProj.project(delta[patchEdgei]);
                const scalar XonMagSqr = Foam::magSqr(Xon);

                if (XonMagSqr > VSMALL)
                {
                    const valueType graddVe =
                        ((valueBound[patchEdgei] - returnValue[own])
                         * (Xon & LeProj[edgei]))
                       / XonMagSqr;

                    laplValue[own] +=
                        (Kex / SfProj[own]) * graddVe;
                }

                edgei++;
            }
        }

        returnValue += laplValue;
    }

    exchangeField.primitiveFieldRef() = exchangeOriginal;
    exchangeField.correctBoundaryConditions();

    return returnValue;
}


// Projected finite-volume Laplacian used by the implicit biharmonic filter.
// The edge contribution is shared with opposite signs by its two adjacent
// faces, so the area-integrated result is conservative. Clipping negative
// non-orthogonal weights keeps this graph Laplacian negative semi-definite,
// which is required by the conjugate-gradient solve below.
inline scalarField projectedLaplacian
(
    const scalarField& value,
    const projectedFaMesh& meshProj,
    const scalarField& SfProj,
    const vectorField& xFacesProj,
    const vectorField& LeProj,
    const labelList& owner,
    const labelList& neighbour,
    areaScalarField& exchangeField
)
{
    scalarField laplValue(value.size(), Zero);

    exchangeField.primitiveFieldRef() = value;
    exchangeField.correctBoundaryConditions();

    for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
    {
        const label own = owner[edgei];
        const label nei = neighbour[edgei];
        const vector Xon = xFacesProj[nei] - xFacesProj[own];
        const scalar XonMagSqr = Foam::magSqr(Xon);

        if (XonMagSqr <= VSMALL)
        {
            continue;
        }

        const scalar edgeWeight = max
        (
            (Xon & LeProj[edgei])/XonMagSqr,
            scalar(0)
        );
        const scalar edgeFlux =
            edgeWeight*(value[nei] - value[own]);

        laplValue[own] += edgeFlux/SfProj[own];
        laplValue[nei] -= edgeFlux/SfProj[nei];
    }

    forAll(meshProj.mesh().boundary(), patchi)
    {
        const faPatch& patch = meshProj.mesh().boundary()[patchi];
        if (not patch.coupled())
        {
            continue;
        }

        const coupledFaPatch& coupledPatch =
            refCast<const coupledFaPatch>(patch);
        const tmp<vectorField> tDelta = coupledPatch.delta();
        const vectorField& delta = tDelta();
        const auto& valueBound = exchangeField.boundaryField()[patchi];

        label edgei = patch.start();
        forAll(patch, patchEdgei)
        {
            const label own = owner[edgei];
            const vector Xon = meshProj.project(delta[patchEdgei]);
            const scalar XonMagSqr = Foam::magSqr(Xon);

            if (XonMagSqr > VSMALL)
            {
                const scalar edgeWeight = max
                (
                    (Xon & LeProj[edgei])/XonMagSqr,
                    scalar(0)
                );
                const scalar edgeFlux =
                    edgeWeight*(valueBound[patchEdgei] - value[own]);

                laplValue[own] += edgeFlux/SfProj[own];
            }

            edgei++;
        }
    }

    return laplValue;
}


inline scalar projectedAreaDot
(
    const scalarField& lhs,
    const scalarField& rhs,
    const scalarField& SfProj
)
{
    scalar localSum = 0;

    forAll(lhs, facei)
    {
        localSum += Foam::mag(SfProj[facei])*lhs[facei]*rhs[facei];
    }

    return returnReduce(localSum, sumOp<scalar>());
}


inline scalarField applyImplicitBiharmonicOperator
(
    const scalarField& value,
    const scalar coefficient,
    const projectedFaMesh& meshProj,
    const scalarField& SfProj,
    const vectorField& xFacesProj,
    const vectorField& LeProj,
    const labelList& owner,
    const labelList& neighbour,
    areaScalarField& exchangeField
)
{
    const scalarField laplValue = projectedLaplacian
    (
        value,
        meshProj,
        SfProj,
        xFacesProj,
        LeProj,
        owner,
        neighbour,
        exchangeField
    );
    const scalarField laplLaplValue = projectedLaplacian
    (
        laplValue,
        meshProj,
        SfProj,
        xFacesProj,
        LeProj,
        owner,
        neighbour,
        exchangeField
    );

    scalarField result(value);
    result += coefficient*laplLaplValue;

    return result;
}


inline scalarField implicitBiharmonicFilter
(
    const scalarField& value,
    const scalar alpha,
    const scalar dl,
    const projectedFaMesh& meshProj,
    const scalarField& SfProj,
    const vectorField& xFacesProj,
    const vectorField& LeProj,
    const labelList& owner,
    const labelList& neighbour,
    areaScalarField& exchangeField
)
{
    if
    (
        alpha <= 0
     || projectedGlobalEdgeCount(meshProj) == 0
    )
    {
        return value;
    }

    // Solve (I + alpha*dl^4*L^2) valueFiltered = value with CG.
    const scalar coefficient = alpha*Foam::sqr(Foam::sqr(dl));
    const scalarField input(value);
    scalarField filtered(input);
    scalarField residual(value.size(), Zero);
    scalarField direction(value.size(), Zero);

    const scalarField applied = applyImplicitBiharmonicOperator
    (
        filtered,
        coefficient,
        meshProj,
        SfProj,
        xFacesProj,
        LeProj,
        owner,
        neighbour,
        exchangeField
    );

    residual = input - applied;
    direction = residual;

    scalar residual2 = projectedAreaDot(residual, residual, SfProj);
    const scalar initialResidual2 = residual2;
    const scalar tolerance = 1e-10;
    const label maxIter = 1000;
    label iter = 0;

    while
    (
        iter < maxIter
     && residual2 > VSMALL
     && Foam::sqrt(residual2/initialResidual2) > tolerance
    )
    {
        const scalarField appliedDirection = applyImplicitBiharmonicOperator
        (
            direction,
            coefficient,
            meshProj,
            SfProj,
            xFacesProj,
            LeProj,
            owner,
            neighbour,
            exchangeField
        );

        const scalar denominator = projectedAreaDot
        (
            direction,
            appliedDirection,
            SfProj
        );

        if (denominator <= VSMALL)
        {
            break;
        }

        const scalar step = residual2/denominator;
        filtered += step*direction;
        residual -= step*appliedDirection;

        const scalar newResidual2 =
            projectedAreaDot(residual, residual, SfProj);
        const scalar beta = newResidual2/residual2;
        direction = residual + beta*direction;
        residual2 = newResidual2;
        iter++;
    }

    const scalar relativeResidual =
        initialResidual2 > VSMALL
      ? Foam::sqrt(residual2/initialResidual2)
      : 0;

    if (relativeResidual > tolerance)
    {
        WarningInFunction
            << "Implicit biharmonic Exner filter stopped after " << iter
            << " iterations with relative residual " << relativeResidual
            << endl;
    }

    // Remove the small residual mass error from the iterative solve.
    scalar initialIntegral = 0;
    scalar filteredIntegral = 0;
    scalar projectedArea = 0;

    forAll(input, facei)
    {
        initialIntegral += SfProj[facei]*input[facei];
        filteredIntegral += SfProj[facei]*filtered[facei];
        projectedArea += SfProj[facei];
    }

    initialIntegral = returnReduce(initialIntegral, sumOp<scalar>());
    filteredIntegral = returnReduce(filteredIntegral, sumOp<scalar>());
    projectedArea = returnReduce(projectedArea, sumOp<scalar>());

    if (Foam::mag(projectedArea) > VSMALL)
    {
        filtered +=
            (initialIntegral - filteredIntegral)/projectedArea;
    }

    return filtered;
}
