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
    Foam::fv::blockDdtScheme

Description
    Abstract base class for ddt schemes for block-coupled matrices.

SourceFiles
    blockDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef blockDdtScheme_H
#define blockDdtScheme_H

#include "tmp.H"
#include "dimensionedType.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
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
                           Class blodkDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class blockDdtScheme
:
    public refCount
{

protected:

    // Protected data

        const fvMesh& mesh_;


    // Private Member Functions

        //- Disallow copy construct
        blockDdtScheme(const blockDdtScheme&);

        //- Disallow default bitwise assignment
        void operator=(const blockDdtScheme&);


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            blockDdtScheme,
            Istream,
            (const fvMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Constructors

        //- Construct from mesh
        blockDdtScheme(const fvMesh& mesh)
        :
            mesh_(mesh)
        {}

        //- Construct from mesh and Istream
        blockDdtScheme(const fvMesh& mesh, Istream&)
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return a pointer to a new ddtScheme created on freestore
        static tmp<blockDdtScheme<Type> > New
        (
            const fvMesh& mesh,
            Istream& schemeData
        );


    // Destructor

        virtual ~blockDdtScheme();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        virtual tmp<blockFvMatrix<Type> > fvmDdt
        (
            GeometricField<Type, fvPatchField, volMesh>&
        ) = 0;

        virtual tmp<blockFvMatrix<Type> > fvmDdt
        (
            const dimensionedScalar&,
            GeometricField<Type, fvPatchField, volMesh>&
        ) = 0;

        virtual tmp<blockFvMatrix<Type> > fvmDdt
        (
            const volScalarField&,
            GeometricField<Type, fvPatchField, volMesh>&
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeBlockFvDdtTypeScheme(SS, Type)                                     \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
blockDdtScheme<Type>::addIstreamConstructorToTable<SS<Type> >                  \
    add##SS##Type##IstreamConstructorToTable_;


#define makeBlockFvDdtScheme(SS)                                               \
                                                                               \
makeBlockFvDdtTypeScheme(SS, scalar)                                           \
makeBlockFvDdtTypeScheme(SS, vector)                                           \
makeBlockFvDdtTypeScheme(SS, sphericalTensor)                                  \
makeBlockFvDdtTypeScheme(SS, symmTensor)                                       \
makeBlockFvDdtTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "blockDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
