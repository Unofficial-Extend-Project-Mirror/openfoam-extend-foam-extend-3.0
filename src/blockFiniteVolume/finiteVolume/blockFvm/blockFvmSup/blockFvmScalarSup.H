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

Namespace
    blockFvm

Description
    Calculate the matrix for implicit and explicit sources.

Author
    Ivor Clifford <ivor.clifford@psu.edu>

SourceFiles
    blockScalarFvmSup.C

\*---------------------------------------------------------------------------*/

#ifndef blockFvmScalarSup_H
#define blockFvmScalarSup_H

#include "volFields.H"
#include "blockFvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvm functions Declaration
\*---------------------------------------------------------------------------*/

namespace blockFvm
{
   // Explicit source

        tmp<blockFvMatrix<scalar> > Su
        (
            const GeometricField<scalar, fvPatchField, volMesh>&,
            GeometricField<scalar, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<scalar> > Su
        (
            const tmp<GeometricField<scalar, fvPatchField, volMesh> >&,
            GeometricField<scalar, fvPatchField, volMesh>&
        );

        zeroField Su
        (
            const zeroField&,
            GeometricField<scalar, fvPatchField, volMesh>&
        );

    // Implicit source

        tmp<blockFvMatrix<scalar> > Sp
        (
            const GeometricField<scalar, fvPatchField, volMesh>&,
            GeometricField<scalar, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<scalar> > Sp
        (
            const tmp<GeometricField<scalar, fvPatchField, volMesh> >&,
            GeometricField<scalar, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<scalar> > Sp
        (
            const dimensioned<scalar>,
            GeometricField<scalar, fvPatchField, volMesh>&
        );

        zeroField Sp
        (
            const zeroField&,
            GeometricField<scalar, fvPatchField, volMesh>&
        );

    // Implicit/Explicit source depending on sign of coefficient

        tmp<blockFvMatrix<scalar> > SuSp
        (
            const GeometricField<scalar, fvPatchField, volMesh>&,
            GeometricField<scalar, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<scalar> > SuSp
        (
            const tmp<GeometricField<scalar, fvPatchField, volMesh> >&,
            GeometricField<scalar, fvPatchField, volMesh>&
        );
        
        zeroField SuSp
        (
            const zeroField&,
            GeometricField<scalar, fvPatchField, volMesh>&
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "blockFvmScalarSup.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
