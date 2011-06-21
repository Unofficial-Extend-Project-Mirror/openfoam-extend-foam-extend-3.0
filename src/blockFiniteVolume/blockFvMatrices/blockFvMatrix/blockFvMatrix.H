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
    blockFvMatrix<Type>

Description
    Block Finite-Volume matrix.

Author
    Ivor Clifford <ivor.clifford@psu.edu>
    
SourceFiles
    blockFvMatrix.C
    blockFvMatrixSolve.C

\*---------------------------------------------------------------------------*/

#ifndef blockFvMatrix_H
#define blockFvMatrix_H

#include "volFields.H"
#include "surfaceFields.H"
#include "BlockLduMatrix.H"
#include "BlockSolverPerformance.H"
#include "tmp.H"
#include "autoPtr.H"
#include "dimensionedTypes.H"
#include "zeroField.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class blockFvMatrix;

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > operator&
(
    const blockFvMatrix<Type>&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > operator&
(
    const blockFvMatrix<Type>&,
    const tmp<DimensionedField<Type, volMesh> >&
);

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > operator&
(
    const tmp<blockFvMatrix<Type> >&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > operator&
(
    const tmp<blockFvMatrix<Type> >&,
    const tmp<DimensionedField<Type, volMesh> >&
);

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > operator&
(
    const blockFvMatrix<Type>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > operator&
(
    const tmp<blockFvMatrix<Type> >&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
Ostream& operator<<(Ostream&, const blockFvMatrix<Type>&);


/*---------------------------------------------------------------------------*\
                           Class blockFvMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class blockFvMatrix
:
    public BlockLduMatrix<Type>
{
    // Private data

        // Reference to GeometricField<Type, fvPatchField, volMesh>
        GeometricField<Type, fvPatchField, volMesh>& psi_;

        //- Dimension set
        dimensionSet dimensions_;

        //- Face flux field for non-orthogonal correction
        mutable GeometricField<Type, fvsPatchField, surfaceMesh>
            *faceFluxCorrectionPtr_;


    // Private member functions

        //- Add patch contribution to internal field
        template<class Type2>
        void addToInternalField
        (
            const unallocLabelList& addr,
            const Field<Type2>& pf,
            Field<Type2>& intf
        ) const;

        template<class Type2>
        void addToInternalField
        (
            const unallocLabelList& addr,
            const tmp<Field<Type2> >& tpf,
            Field<Type2>& intf
        ) const;

        template<class Type2>
        void addToInternalField
        (
            const unallocLabelList& addr,
            const CoeffField<Type2>& pf,
            CoeffField<Type2>& intf
        ) const;
        
        //- Subtract patch contribution from internal field
        template<class Type2>
        void subtractFromInternalField
        (
            const unallocLabelList& addr,
            const Field<Type2>& pf,
            Field<Type2>& intf
        ) const;

        template<class Type2>
        void subtractFromInternalField
        (
            const unallocLabelList& addr,
            const tmp<Field<Type2> >& tpf,
            Field<Type2>& intf
        ) const;

        template<class Type2>
        void subtractFromInternalField
        (
            const unallocLabelList& addr,
            const CoeffField<Type2>& pf,
            CoeffField<Type2>& intf
        ) const;

public:

    ClassName("blockFvMatrix");


    // Constructors

        //- Construct given a field to solve for
        blockFvMatrix
        (
            GeometricField<Type, fvPatchField, volMesh>&,
            const dimensionSet&
        );

        //- Construct as copy
        blockFvMatrix(const blockFvMatrix<Type>&);

        //- Construct as copy of tmp<blockFvMatrix<Type> > deleting argument
#       ifdef ConstructFromTmp
        blockFvMatrix(const tmp<blockFvMatrix<Type> >&);
#       endif

        //- Construct from Istream given field to solve for
        blockFvMatrix(GeometricField<Type, fvPatchField, volMesh>&, Istream&);


    // Destructor

        virtual ~blockFvMatrix();


    // Member functions

        // Access

            const GeometricField<Type, fvPatchField, volMesh>& psi() const
            {
                return psi_;
            }

            GeometricField<Type, fvPatchField, volMesh>& psi()
            {
                return psi_;
            }

            const dimensionSet& dimensions() const
            {
                return dimensions_;
            }

            //- Declare return type of the faceFluxCorrectionPtr() function
            typedef GeometricField<Type, fvsPatchField, surfaceMesh>
                *surfaceTypeFieldPtr;

            //- Return pointer to face-flux non-orthogonal correction field
            surfaceTypeFieldPtr& faceFluxCorrectionPtr()
            {
                return faceFluxCorrectionPtr_;
            }

        // Matrix completion functionality

            void addBoundaryDiag
            (
                CoeffField<Type>& diag
            ) const;

            void addBoundarySource
            (
                Field<Type>& source,
                const bool couples=true
            ) const;


        // Operations

            //- Solve returning the solution statistics.
            //  Solver controls read from Istream
            BlockSolverPerformance<Type> solve(const dictionary&);

            //- Solve returning the solution statistics.
            //  Solver controls read from fvSolution
            BlockSolverPerformance<Type> solve();

            //- Return the matrix residual
            tmp<Field<Type> > residual() const;

            //- Return the matrix diagonal
            tmp<CoeffField<Type> > D() const;

            //- Return the H operation source
            tmp<GeometricField<Type, fvPatchField, volMesh> > H() const;

    // Member operators

        void operator=(const blockFvMatrix<Type>&);
        void operator=(const tmp<blockFvMatrix<Type> >&);

        void negate();

        void operator+=(const blockFvMatrix<Type>&);
        void operator+=(const tmp<blockFvMatrix<Type> >&);

        void operator-=(const blockFvMatrix<Type>&);
        void operator-=(const tmp<blockFvMatrix<Type> >&);

        void operator+=
        (
            const GeometricField<Type, fvPatchField,volMesh>&
        );
        void operator+=
        (
            const tmp<GeometricField<Type, fvPatchField,volMesh> >&
        );

        void operator-=
        (
            const GeometricField<Type, fvPatchField,volMesh>&
        );
        void operator-=
        (
            const tmp<GeometricField<Type, fvPatchField,volMesh> >&
        );

        void operator+=(const dimensioned<Type>&);
        void operator-=(const dimensioned<Type>&);

        void operator+=(const zeroField&);
        void operator-=(const zeroField&);

        void operator*=(const DimensionedField<scalar, volMesh>&);
        void operator*=(const tmp<DimensionedField<scalar, volMesh> >&);
        void operator*=(const tmp<volScalarField>&);

        void operator*=(const dimensioned<scalar>&);


    // Friend operators

        friend tmp<GeometricField<Type,fvPatchField,volMesh> > operator& <Type>
        (
            const blockFvMatrix<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        friend tmp<GeometricField<Type,fvPatchField,volMesh> > operator& <Type>
        (
            const blockFvMatrix<Type>&,
            const tmp<DimensionedField<Type, volMesh> >&
        );

        friend tmp<GeometricField<Type,fvPatchField,volMesh> > operator& <Type>
        (
            const tmp<blockFvMatrix<Type> >&,
            const DimensionedField<Type, volMesh>&
        );

        friend tmp<GeometricField<Type,fvPatchField,volMesh> > operator& <Type>
        (
            const tmp<blockFvMatrix<Type> >&,
            const tmp<DimensionedField<Type, volMesh> >&
        );

        friend tmp<GeometricField<Type,fvPatchField,volMesh> > operator& <Type>
        (
            const blockFvMatrix<Type>&,
            const tmp<GeometricField<Type,fvPatchField,volMesh> >&
        );

        friend tmp<GeometricField<Type,fvPatchField,volMesh> > operator& <Type>
        (
            const tmp<blockFvMatrix<Type> >&,
            const tmp<GeometricField<Type, fvPatchField, volMesh> >&
        );


    // Ostream operator

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const blockFvMatrix<Type>&
        );        
};


// * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * //

template<class Type>
void checkMethod
(
    const blockFvMatrix<Type>&,
    const blockFvMatrix<Type>&,
    const char*
);

template<class Type>
void checkMethod
(
    const blockFvMatrix<Type>&,
    const GeometricField<Type, fvPatchField, volMesh>&,
    const char*
);

template<class Type>
void checkMethod
(
    const blockFvMatrix<Type>&,
    const dimensioned<Type>&,
    const char*
);


//- Solve returning the solution statistics given convergence tolerance
//  Solver controls read from Istream
template<class Type>
BlockSolverPerformance<Type> solve(blockFvMatrix<Type>&, const dictionary&);


//- Solve returning the solution statistics given convergence tolerance,
//  deleting temporary matrix after solution.
//  Solver controls read from Istream
template<class Type>
BlockSolverPerformance<Type> solve(const tmp<blockFvMatrix<Type> >&, const dictionary&);


//- Solve returning the solution statistics given convergence tolerance
//  Solver controls read fvSolution
template<class Type>
BlockSolverPerformance<Type> solve(blockFvMatrix<Type>&);


//- Solve returning the solution statistics given convergence tolerance,
//  deleting temporary matrix after solution.
//  Solver controls read fvSolution
template<class Type>
BlockSolverPerformance<Type> solve(const tmp<blockFvMatrix<Type> >&);


// * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * //

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const blockFvMatrix<Type>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const tmp<blockFvMatrix<Type> >&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const blockFvMatrix<Type>&,
    const blockFvMatrix<Type>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const tmp<blockFvMatrix<Type> >&,
    const blockFvMatrix<Type>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const blockFvMatrix<Type>&,
    const tmp<blockFvMatrix<Type> >&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const tmp<blockFvMatrix<Type> >&,
    const tmp<blockFvMatrix<Type> >&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const blockFvMatrix<Type>&,
    const blockFvMatrix<Type>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const tmp<blockFvMatrix<Type> >&,
    const blockFvMatrix<Type>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const blockFvMatrix<Type>&,
    const tmp<blockFvMatrix<Type> >&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const tmp<blockFvMatrix<Type> >&,
    const tmp<blockFvMatrix<Type> >&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const blockFvMatrix<Type>&,
    const blockFvMatrix<Type>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const tmp<blockFvMatrix<Type> >&,
    const blockFvMatrix<Type>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const blockFvMatrix<Type>&,
    const tmp<blockFvMatrix<Type> >&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const tmp<blockFvMatrix<Type> >&,
    const tmp<blockFvMatrix<Type> >&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const blockFvMatrix<Type>&,
    const GeometricField<Type, fvPatchField, volMesh>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const tmp<blockFvMatrix<Type> >&,
    const GeometricField<Type, fvPatchField, volMesh>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const blockFvMatrix<Type>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const tmp<blockFvMatrix<Type> >&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const GeometricField<Type, fvPatchField, volMesh>&,
    const blockFvMatrix<Type>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const GeometricField<Type, fvPatchField, volMesh>&,
    const tmp<blockFvMatrix<Type> >&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&,
    const blockFvMatrix<Type>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&,
    const tmp<blockFvMatrix<Type> >&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const blockFvMatrix<Type>&,
    const GeometricField<Type, fvPatchField, volMesh>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const tmp<blockFvMatrix<Type> >&,
    const GeometricField<Type, fvPatchField, volMesh>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const blockFvMatrix<Type>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const tmp<blockFvMatrix<Type> >&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const GeometricField<Type, fvPatchField, volMesh>&,
    const blockFvMatrix<Type>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const GeometricField<Type, fvPatchField, volMesh>&,
    const tmp<blockFvMatrix<Type> >&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&,
    const blockFvMatrix<Type>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&,
    const tmp<blockFvMatrix<Type> >&
);


template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const blockFvMatrix<Type>&,
    const dimensioned<Type>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const tmp<blockFvMatrix<Type> >&,
    const dimensioned<Type>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const dimensioned<Type>&,
    const blockFvMatrix<Type>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const dimensioned<Type>&,
    const tmp<blockFvMatrix<Type> >&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const blockFvMatrix<Type>&,
    const dimensioned<Type>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const tmp<blockFvMatrix<Type> >&,
    const dimensioned<Type>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const dimensioned<Type>&,
    const blockFvMatrix<Type>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const dimensioned<Type>&,
    const tmp<blockFvMatrix<Type> >&
);


template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const blockFvMatrix<Type>&,
    const GeometricField<Type, fvPatchField, volMesh>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const tmp<blockFvMatrix<Type> >&,
    const GeometricField<Type, fvPatchField, volMesh>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const blockFvMatrix<Type>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const tmp<blockFvMatrix<Type> >&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const blockFvMatrix<Type>&,
    const dimensioned<Type>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const tmp<blockFvMatrix<Type> >&,
    const dimensioned<Type>&
);


template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const blockFvMatrix<Type>&,
    const zeroField&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const tmp<blockFvMatrix<Type> >&,
    const zeroField&
);


template<class Type>
tmp<blockFvMatrix<Type> > operator*
(
    const volScalarField&,
    const blockFvMatrix<Type>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator*
(
    const volScalarField&,
    const tmp<blockFvMatrix<Type> >&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator*
(
    const tmp<volScalarField>&,
    const blockFvMatrix<Type>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator*
(
    const tmp<volScalarField>&,
    const tmp<blockFvMatrix<Type> >&
);


template<class Type>
tmp<blockFvMatrix<Type> > operator*
(
    const dimensioned<scalar>&,
    const blockFvMatrix<Type>&
);

template<class Type>
tmp<blockFvMatrix<Type> > operator*
(
    const dimensioned<scalar>&,
    const tmp<blockFvMatrix<Type> >&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "blockFvMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
