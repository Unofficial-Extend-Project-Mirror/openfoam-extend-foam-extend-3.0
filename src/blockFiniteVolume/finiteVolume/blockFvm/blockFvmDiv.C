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

#include "blockFvmDiv.H"
#include "fvMesh.H"
#include "blockFvMatrix.H"
#include "blockConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace blockFvm
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<blockFvMatrix<Type> >
div
(
    const surfaceScalarField& flux,
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::blockConvectionScheme<Type>::New
    (
        vf.mesh(),
        flux,
        vf.mesh().divScheme(name)
    )().fvmDiv(flux, vf);
}


template<class Type>
tmp<blockFvMatrix<Type> >
div
(
    const tmp<surfaceScalarField>& tflux,
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    tmp<blockFvMatrix<Type> > Div(blockFvm::div(tflux(), vf, name));
    tflux.clear();
    return Div;
}


template<class Type>
tmp<blockFvMatrix<Type> >
div
(
    const surfaceScalarField& flux,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return blockFvm::div(flux, vf, "div("+flux.name()+','+vf.name()+')');
}

template<class Type>
tmp<blockFvMatrix<Type> >
div
(
    const tmp<surfaceScalarField>& tflux,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<blockFvMatrix<Type> > Div(blockFvm::div(tflux(), vf));
    tflux.clear();
    return Div;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockFvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
