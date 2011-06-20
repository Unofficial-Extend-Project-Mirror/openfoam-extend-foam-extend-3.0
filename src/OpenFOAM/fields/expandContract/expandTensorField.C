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

Description
    Global functions for expansion and contraction of tensor field
    to diagonal type

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "tensorField.H"
#include "diagTensorField.H"
#include "sphericalTensorField.H"
#include "expandTensor.H"
#include "expandDiagTensor.H"
#include "expandSphericalTensor.H"

#define TEMPLATE
#include "FieldFunctionsM.C"

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void contractScalar(Field<scalar>& res, const UList<tensor>& f)
{
    forAll (res, i)
    {
        contractScalar(res[i], f[i]);
    }
}


void contractScalar(Field<scalar>& res, const UList<diagTensor>& f)
{
    forAll (res, i)
    {
        contractScalar(res[i], f[i]);
    }
}


void contractScalar(Field<scalar>& res, const UList<sphericalTensor>& f)
{
    forAll (res, i)
    {
        contractScalar(res[i], f[i]);
    }
}


void contractLinear(Field<vector>& res, const UList<tensor>& f)
{
    forAll (res, i)
    {
        contractLinear(res[i], f[i]);
    }
}


void contractLinear(Field<vector>& res, const UList<diagTensor>& f)
{
    forAll (res, i)
    {
        contractLinear(res[i], f[i]);
    }
}


void expandScalar(Field<vector>& res, const UList<scalar>& f)
{
    forAll (res, i)
    {
        expandScalar(res[i], f[i]);
    }
}


void expandScalar(Field<tensor>& res, const UList<scalar>& f)
{
    forAll (res, i)
    {
        expandScalar(res[i], f[i]);
    }
}


void expandScalar(Field<diagTensor>& res, const UList<scalar>& f)
{
    forAll (res, i)
    {
        expandScalar(res[i], f[i]);
    }
}


void expandScalar(Field<sphericalTensor>& res, const UList<scalar>& f)
{
    forAll (res, i)
    {
        expandScalar(res[i], f[i]);
    }
}


void expandLinear(Field<tensor>& res, const UList<vector>& f)
{
    forAll (res, i)
    {
        expandLinear(res[i], f[i]);
    }
}


void expandLinear(Field<diagTensor>& res, const UList<vector>& f)
{
    forAll (res, i)
    {
        expandLinear(res[i], f[i]);
    }
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
