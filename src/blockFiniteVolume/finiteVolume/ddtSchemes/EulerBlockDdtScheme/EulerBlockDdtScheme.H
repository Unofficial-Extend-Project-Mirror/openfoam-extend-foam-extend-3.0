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

Class
    Foam::fv::EulerBlockDdtScheme

Description
    Basic first-order Euler implicit/explicit ddt using only the current and
    previous time-step values.

SourceFiles
    EulerBlockDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef EulerBlockDdtScheme_H
#define EulerBlockDdtScheme_H

#include "blockDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class EulerBlockDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class EulerBlockDdtScheme
:
    public blockDdtScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        EulerBlockDdtScheme(const EulerBlockDdtScheme&);

        //- Disallow default bitwise assignment
        void operator=(const EulerBlockDdtScheme&);


public:

    //- Runtime type information
    TypeName("Euler");


    // Constructors

        //- Construct from mesh
        EulerBlockDdtScheme(const fvMesh& mesh)
        :
            blockDdtScheme<Type>(mesh)
        {}

        //- Construct from mesh and Istream
        EulerBlockDdtScheme(const fvMesh& mesh, Istream& is)
        :
            blockDdtScheme<Type>(mesh, is)
        {}


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return fv::blockDdtScheme<Type>::mesh();
        }

        tmp<blockFvMatrix<Type> > fvmDdt
        (
            GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<Type> > fvmDdt
        (
            const dimensionedScalar&,
            GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<blockFvMatrix<Type> > fvmDdt
        (
            const volScalarField&,
            GeometricField<Type, fvPatchField, volMesh>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "EulerBlockDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
