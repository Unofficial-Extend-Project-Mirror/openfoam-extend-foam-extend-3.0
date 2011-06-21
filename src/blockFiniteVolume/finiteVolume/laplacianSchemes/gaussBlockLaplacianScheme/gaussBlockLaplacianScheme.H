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
    Foam::fv::gaussBlockLaplacianScheme

Description
    Basic second-order laplacian using face-gradients and Gauss' theorem
    for block-coupled matrices.

SourceFiles
    gaussBlockLaplacianScheme.C

\*---------------------------------------------------------------------------*/

#ifndef gaussBlockLaplacianScheme_H
#define gaussBlockLaplacianScheme_H

#include "blockLaplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class gaussBlockLaplacianScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GType>
class gaussBlockLaplacianScheme
:
    public fv::blockLaplacianScheme<Type, GType>
{
    // Private Member Functions

        tmp<blockFvMatrix<Type> > fvmLaplacianUncorrected
        (
            const surfaceScalarField& gammaMagSf,
            GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > gammaSnGradCorr
        (
            const surfaceVectorField& SfGammaCorr,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        //- Disallow default bitwise copy construct
        gaussBlockLaplacianScheme(const gaussBlockLaplacianScheme&);

        //- Disallow default bitwise assignment
        void operator=(const gaussBlockLaplacianScheme&);


public:

    //- Runtime type information
    TypeName("Gauss");


    // Constructors

        //- Construct null
        gaussBlockLaplacianScheme(const fvMesh& mesh)
        :
            blockLaplacianScheme<Type, GType>(mesh)
        {}

        //- Construct from Istream
        gaussBlockLaplacianScheme(const fvMesh& mesh, Istream& is)
        :
            blockLaplacianScheme<Type, GType>(mesh, is)
        {}

        //- Construct from mesh, interpolation and snGradScheme schemes
        gaussBlockLaplacianScheme
        (
            const fvMesh& mesh,
            const tmp<surfaceInterpolationScheme<GType> >& igs,
            const tmp<snGradScheme<Type> >& sngs
        )
        :
            blockLaplacianScheme<Type, GType>(mesh, igs, sngs)
        {}


    // Destructor

        virtual ~gaussBlockLaplacianScheme()
        {}


    // Member Functions

        tmp<blockFvMatrix<Type> > fvmLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            GeometricField<Type, fvPatchField, volMesh>&
        );
};


// Use macros to emulate partial-specialisation of the the Laplacian functions
// for scalar diffusivity gamma

#define defineBlockFvmLaplacianScalarGamma(Type)                            \
                                                                            \
template<>                                                                  \
tmp<blockFvMatrix<Type> >                                                   \
gaussBlockLaplacianScheme<Type, scalar>::fvmLaplacian                       \
(                                                                           \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>&,              \
    GeometricField<Type, fvPatchField, volMesh>&                            \
);


defineBlockFvmLaplacianScalarGamma(scalar);
defineBlockFvmLaplacianScalarGamma(vector);
defineBlockFvmLaplacianScalarGamma(sphericalTensor);
defineBlockFvmLaplacianScalarGamma(symmTensor);
defineBlockFvmLaplacianScalarGamma(tensor);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "gaussBlockLaplacianScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
