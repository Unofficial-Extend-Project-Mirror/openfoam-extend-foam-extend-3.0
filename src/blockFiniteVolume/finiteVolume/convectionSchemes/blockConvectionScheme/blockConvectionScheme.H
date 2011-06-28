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
    Foam::fv::blockConvectionScheme

Description
    Abstract base class for convection schemes.

SourceFiles
    blockConvectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef blockConvectionScheme_H
#define blockConvectionScheme_H

#include "tmp.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "multivariateSurfaceInterpolationScheme.H"

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
                           Class blockConvectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class blockConvectionScheme
:
    public refCount
{
    // Private data

        const fvMesh& mesh_;


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            blockConvectionScheme,
            Istream,
            (
                const fvMesh& mesh,
                const surfaceScalarField& faceFlux,
                Istream& schemeData
            ),
            (mesh, faceFlux, schemeData)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            blockConvectionScheme,
            Multivariate,
            (
                const fvMesh& mesh,
                const typename multivariateSurfaceInterpolationScheme<Type>::
                    fieldTable& fields,
                const surfaceScalarField& faceFlux,
                Istream& schemeData
            ),
            (mesh, fields, faceFlux, schemeData)
        );


    // Constructors

        //- Copy construct
        blockConvectionScheme(const blockConvectionScheme&);

        //- Construct from mesh, flux and Istream
        blockConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField&
        )
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return a pointer to a new blockConvectionScheme created on freestore
        static tmp<blockConvectionScheme<Type> > New
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        );


        //- Return a pointer to a new multivariate blockConvectionScheme
        //  created on freestore
        static tmp<blockConvectionScheme<Type> > New
        (
            const fvMesh& mesh,
            const typename multivariateSurfaceInterpolationScheme<Type>::
                fieldTable& fields,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        );


    // Destructor

        virtual ~blockConvectionScheme();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        virtual tmp<blockFvMatrix<Type> > fvmDiv
        (
            const surfaceScalarField&,
            GeometricField<Type, fvPatchField, volMesh>&
        ) const = 0;


    // Member operators

        void operator=(const blockConvectionScheme<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeBlockFvConvectionTypeScheme(SS, Type)                              \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
blockConvectionScheme<Type>::addIstreamConstructorToTable<SS<Type> >           \
    add##SS##Type##IstreamConstructorToTable_;


#define makeBlockFvConvectionScheme(SS)                                        \
                                                                               \
makeBlockFvConvectionTypeScheme(SS, scalar)                                    \
makeBlockFvConvectionTypeScheme(SS, vector)                                    \
makeBlockFvConvectionTypeScheme(SS, sphericalTensor)                           \
makeBlockFvConvectionTypeScheme(SS, symmTensor)                                \
makeBlockFvConvectionTypeScheme(SS, tensor)


#define makeMultivariateBlockFvConvectionTypeScheme(SS, Type)                  \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
blockConvectionScheme<Type>::addMultivariateConstructorToTable<SS<Type> >      \
    add##SS##Type##MultivariateConstructorToTable_;


#define makeMultivariateBlockFvConvectionScheme(SS)                            \
                                                                               \
makeMultivariateBlockFvConvectionTypeScheme(SS, scalar)                        \
makeMultivariateBlockFvConvectionTypeScheme(SS, vector)                        \
makeMultivariateBlockFvConvectionTypeScheme(SS, sphericalTensor)               \
makeMultivariateBlockFvConvectionTypeScheme(SS, symmTensor)                    \
makeMultivariateBlockFvConvectionTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "blockConvectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
