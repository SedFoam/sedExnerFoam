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

template<class Field>
inline Field filterField(
    const label N,
    const scalar Kex,
    const Field& value,
    const projectedFaMesh& meshProj,
    const scalarField& SfProj,
    const vectorField& xFacesProj,
    const vectorField& LeProj,
    const labelList& owner,
    const labelList& neighbour)
{
    typedef typename Field::value_type valueType;

    Field returnValue(value);

    for (label j = 0; j < N; j++)
    {
        Field laplValue(value);
        laplValue = Zero;

        for (label edgei = 0; edgei < meshProj.nInternalEdges(); edgei++)
        {
            const vector Xon =
                xFacesProj[neighbour[edgei]] - xFacesProj[owner[edgei]];

            const valueType graddVe =
                (((returnValue[neighbour[edgei]] - returnValue[owner[edgei]])
                  * (Xon & LeProj[edgei]))
                 / Foam::magSqr(Xon));

            laplValue[owner[edgei]] +=
                (Kex / SfProj[owner[edgei]]) * graddVe;
            laplValue[neighbour[edgei]] -=
                (Kex / SfProj[neighbour[edgei]]) * graddVe;
        }

        returnValue += laplValue;
    }

    return returnValue;
}
