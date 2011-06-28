/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "cyclicFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "scalarCoeffField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Update result field based on interface functionality
template<>
void cyclicFvPatchField<scalar>::updateInterfaceMatrix
(
    const Field<scalar>& psiInternal,
    Field<scalar>& result,
    const BlockLduMatrix<scalar>& m,
    const CoeffField<scalar>& coeffs,
    const Pstream::commsTypes commsType
) const
{
    Field<scalar> pnf(this->size());

    label sizeby2 = this->size()/2;
    const unallocLabelList& faceCells = cyclicPatch_.faceCells();

    for (label facei=0; facei<sizeby2; facei++)
    {
        pnf[facei] = psiInternal[faceCells[facei + sizeby2]];
        pnf[facei + sizeby2] = psiInternal[faceCells[facei]];
    }
    
    // Multiply the field by coefficients and add into the result
    forAll(faceCells, elemI)
    {
        result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
    }
}


//- Update result field based on interface functionality
template<>
void cyclicFvPatchField<vector>::updateInterfaceMatrix
(
    const Field<vector>& psiInternal,
    Field<vector>& result,
    const BlockLduMatrix<vector>& m,
    const CoeffField<vector>& coeffs,
    const Pstream::commsTypes commsType
) const
{
    Field<vector> pnf(this->size());

    label sizeby2 = this->size()/2;
    const unallocLabelList& faceCells = cyclicPatch_.faceCells();

    for (label facei=0; facei<sizeby2; facei++)
    {
        pnf[facei] = psiInternal[faceCells[facei + sizeby2]];
        pnf[facei + sizeby2] = psiInternal[faceCells[facei]];
    }

// Ivor Clifford
// TODO
    // Transform according to the transformation tensors
//    transformCoupleField(pnf, cmpt);

    if (coeffs.activeType() == blockCoeffBase::SCALAR)                          
    {                                                                           
        pnf = coeffs.asScalar() * pnf;                                          
    }                                                                           
    else if (coeffs.activeType() == blockCoeffBase::LINEAR)                     
    {                                                                           
        pnf = cmptMultiply(coeffs.asLinear(), pnf);                             
    }                                                                           
    else if (coeffs.activeType() == blockCoeffBase::SQUARE)
    {
        pnf = coeffs.asSquare() & pnf;
    }

    // Multiply the field by coefficients and add into the result
    forAll(faceCells, elemI)
    {
        result[faceCells[elemI]] -= pnf[elemI];
    }
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makePatchFields(cyclic);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
