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
    blockVectorFvmSup.C

\*---------------------------------------------------------------------------*/

#ifndef blockFvmVectorSup_H
#define blockFvmVectorSup_H

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

        tmp<blockFvMatrix<vector> > Su
        (
            const GeometricField<tensor, fvPatchField, volMesh>&,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<vector> > Su
        (
            const GeometricField<diagTensor, fvPatchField, volMesh>&,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<vector> > Su
        (
            const GeometricField<sphericalTensor, fvPatchField, volMesh>&,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<vector> > Su
        (
            const tmp<GeometricField<tensor, fvPatchField, volMesh> >&,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<vector> > Su
        (
            const tmp<GeometricField<diagTensor, fvPatchField, volMesh> >&,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<vector> > Su
        (
            const tmp<GeometricField<sphericalTensor, fvPatchField, volMesh> >&,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        zeroField Su
        (
            const zeroField&,
            GeometricField<vector, fvPatchField, volMesh>&
        );


    // Implicit source

        tmp<blockFvMatrix<vector> > Sp
        (
            const GeometricField<tensor, fvPatchField, volMesh>&,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<vector> > Sp
        (
            const GeometricField<diagTensor, fvPatchField, volMesh>&,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<vector> > Sp
        (
            const GeometricField<sphericalTensor, fvPatchField, volMesh>&,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<vector> > Sp
        (
            const tmp<GeometricField<tensor, fvPatchField, volMesh> >&,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<vector> > Sp
        (
            const tmp<GeometricField<diagTensor, fvPatchField, volMesh> >&,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<vector> > Sp
        (
            const tmp<GeometricField<sphericalTensor, fvPatchField, volMesh> >&,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<vector> > Sp
        (
            const dimensioned<tensor>,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<vector> > Sp
        (
            const dimensioned<diagTensor>,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<vector> > Sp
        (
            const dimensioned<sphericalTensor>,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        zeroField Sp
        (
            const zeroField&,
            GeometricField<vector, fvPatchField, volMesh>&
        );


    // Implicit/Explicit source depending on sign of coefficient

        tmp<blockFvMatrix<vector> > SuSp
        (
            const GeometricField<tensor, fvPatchField, volMesh>&,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<vector> > SuSp
        (
            const GeometricField<diagTensor, fvPatchField, volMesh>&,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<vector> > SuSp
        (
            const GeometricField<sphericalTensor, fvPatchField, volMesh>&,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<vector> > SuSp
        (
            const tmp<GeometricField<tensor, fvPatchField, volMesh> >&,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<vector> > SuSp
        (
            const tmp<GeometricField<diagTensor, fvPatchField, volMesh> >&,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<vector> > SuSp
        (
            const tmp<GeometricField<sphericalTensor, fvPatchField, volMesh> >&,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        zeroField SuSp
        (
            const zeroField&,
            GeometricField<vector, fvPatchField, volMesh>&
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "blockFvmVectorSup.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
