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
    Foam::fv::blockLaplacianScheme

Description
    Abstract base class for laplacian schemes for block-coupled matrices.

SourceFiles
    blockLaplacianScheme.C

\*---------------------------------------------------------------------------*/

#ifndef blockLaplacianScheme_H
#define blockLaplacianScheme_H

#include "tmp.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "linear.H"
#include "correctedSnGrad.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class blockFvMatrix;

class fvMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class blockLaplacianScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GType>
class blockLaplacianScheme
:
    public refCount
{

protected:

    // Protected data

        const fvMesh& mesh_;
        tmp<surfaceInterpolationScheme<GType> > tinterpGammaScheme_;
        tmp<snGradScheme<Type> > tsnGradScheme_;


    // Private Member Functions

        //- Disallow copy construct
        blockLaplacianScheme(const blockLaplacianScheme&);

        //- Disallow default bitwise assignment
        void operator=(const blockLaplacianScheme&);


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            blockLaplacianScheme,
            Istream,
            (const fvMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Constructors

        //- Construct from mesh
        blockLaplacianScheme(const fvMesh& mesh)
        :
            mesh_(mesh),
            tinterpGammaScheme_(new linear<GType>(mesh)),
            tsnGradScheme_(new correctedSnGrad<Type>(mesh))
        {}

        //- Construct from mesh and Istream
        blockLaplacianScheme(const fvMesh& mesh, Istream& is)
        :
            mesh_(mesh),
            tinterpGammaScheme_(NULL),
            tsnGradScheme_(NULL)
        {
            tinterpGammaScheme_ = tmp<surfaceInterpolationScheme<GType> >
            (
                surfaceInterpolationScheme<GType>::New(mesh, is)
            );

            tsnGradScheme_ = tmp<snGradScheme<Type> >
            (
                snGradScheme<Type>::New(mesh, is)
            );
        }

        //- Construct from mesh, interpolation and snGradScheme schemes
        blockLaplacianScheme
        (
            const fvMesh& mesh,
            const tmp<surfaceInterpolationScheme<GType> >& igs,
            const tmp<snGradScheme<Type> >& sngs
        )
        :
            mesh_(mesh),
            tinterpGammaScheme_(igs),
            tsnGradScheme_(sngs)
        {}


    // Selectors

        //- Return a pointer to a new blockLaplacianScheme created on freestore
        static tmp<blockLaplacianScheme<Type, GType> > New
        (
            const fvMesh& mesh,
            Istream& schemeData
        );


    // Destructor

        virtual ~blockLaplacianScheme();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        virtual tmp<blockFvMatrix<Type> > fvmLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            GeometricField<Type, fvPatchField, volMesh>&
        ) = 0;

        virtual tmp<blockFvMatrix<Type> > fvmLaplacian
        (
            const GeometricField<GType, fvPatchField, volMesh>&,
            GeometricField<Type, fvPatchField, volMesh>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeBlockFvLaplacianTypeScheme(SS, Type, GType)                        \
                                                                               \
    typedef SS<Type, GType> SS##Type##GType;                                   \
    defineNamedTemplateTypeNameAndDebug(SS##Type##GType, 0);                   \
                                                                               \
blockLaplacianScheme<Type, GType>::                                            \
        addIstreamConstructorToTable<SS<Type, GType> >                         \
    add##SS##Type##GType##IstreamConstructorToTable_;


#define makeBlockFvLaplacianScheme(SS)                                         \
                                                                               \
makeBlockFvLaplacianTypeScheme(SS, scalar, scalar)                             \
makeBlockFvLaplacianTypeScheme(SS, scalar, symmTensor)                         \
makeBlockFvLaplacianTypeScheme(SS, scalar, tensor)                             \
makeBlockFvLaplacianTypeScheme(SS, vector, scalar)                             \
makeBlockFvLaplacianTypeScheme(SS, sphericalTensor, scalar)                    \
makeBlockFvLaplacianTypeScheme(SS, symmTensor, scalar)                         \
makeBlockFvLaplacianTypeScheme(SS, symmTensor, symmTensor)                     \
makeBlockFvLaplacianTypeScheme(SS, tensor, scalar)                             \


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "blockLaplacianScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
