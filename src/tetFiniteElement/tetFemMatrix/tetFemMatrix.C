/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    Tetrahedral Finite Element matrix member functions and operators

\*---------------------------------------------------------------------------*/

#include "PstreamReduceOps.H"

#include "tetFemMatrix.H"
#include "tetPointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const label tetFemMatrix<Type>::fixFillIn = 4;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
tetFemMatrix<Type>::tetFemMatrix
(
    GeometricField<Type, tetPolyPatchField, tetPointMesh>& x,
    const dimensionSet& ds
)
:
    BlockLduMatrix<Type>(x.mesh()),
    x_(x),
    dimensions_(ds),
    b_(x.size(), pTraits<Type>::zero),
    boundaryConditionsSet_(false),
    fixedEqns_(Foam::min(x.mesh().lduAddr().size()/fixFillIn, 100))
{
    if (debug)
    {
        Info<< "tetFemMatrix<Type>(GeometricField<Type, tetPolyPatchField, "
            << "tetPointMesh>&, const dimensionSet&) : "
            << "constructing tetFemMatrix<Type> for field " << x_.name()
            << endl;
    }
}


template<class Type>
tetFemMatrix<Type>::tetFemMatrix(const tetFemMatrix<Type>& tetFem)
:
    BlockLduMatrix<Type>(tetFem),
    x_(tetFem.x_),
    dimensions_(tetFem.dimensions_),
    b_(tetFem.b_),
    boundaryConditionsSet_(false),
    fixedEqns_(Foam::min(x_.mesh().lduAddr().size()/fixFillIn, 100))
{
    if (debug)
    {
        InfoIn("tetFemMatrix<Type>::tetFemMatrix(const tetFemMatrix<Type>&)")
            << "copying tetFemMatrix<Type> for field " << x_.name()
            << endl;
    }
}


template<class Type>
tetFemMatrix<Type>::~tetFemMatrix()
{
    if (debug)
    {
        Info<< "tetFemMatrix<Type>::~tetFemMatrix<Type>() : "
            << "destroying tetFemMatrix<Type> for field " << x_.name()
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void tetFemMatrix<Type>::addConstraint
(
    const label vertex,
    const Type& value
)
{
    ConstraintType cp(vertex, value);

    if (!fixedEqns_.found(vertex))
    {
        fixedEqns_.insert(vertex, cp);
    }
    else
    {
        WarningIn
        (
            "void tetFemMatrix<Type>::addConstraint(const label vertex, "
            "const Type& value)"
        )   << "Adding constraint on an already constrained point."
            << "  Point: " << vertex
            << endl;

        fixedEqns_[vertex].combine(cp);
    }
}


template<class Type>
void tetFemMatrix<Type>::relax(const scalar alpha)
{
    this->relax(x_, b_, alpha);
}


template<class Type>
void tetFemMatrix<Type>::relax()
{
    scalar alpha = 0;

    if (x_.mesh().solutionDict().relax(x_.name()))
    {
        alpha = x_.mesh().solutionDict().relaxationFactor(x_.name());
    }

    if (alpha > 0)
    {
        relax(alpha);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void tetFemMatrix<Type>::operator=(const tetFemMatrix<Type>& tetFem)
{
    if (this == &tetFem)
    {
        FatalErrorIn
        (
            "tetFemMatrix<Type>::operator=(const tetFemMatrix<Type>&)"
        )   << "attempted assignment to self"
            << abort(FatalError);
    }

    if (&x_ != &(tetFem.x_))
    {
        FatalErrorIn
        (
            "tetFemMatrix<Type>::operator=(const tetFemMatrix<Type>&)"
        )   << "different fields"
            << abort(FatalError);
    }

    BlockLduMatrix<Type>::operator=(tetFem);
    b_ = tetFem.b_;
    boundaryConditionsSet_ = false;
    fixedEqns_.clear();
}


template<class Type>
void tetFemMatrix<Type>::operator=(const tmp<tetFemMatrix<Type> >& ttetFem)
{
    operator=(ttetFem());
    ttetFem.clear();
}


template<class Type>
void tetFemMatrix<Type>::negate()
{
    BlockLduMatrix<Type>::negate();
    b_.negate();
}


template<class Type>
void tetFemMatrix<Type>::operator+=(const tetFemMatrix<Type>& tetFem)
{
    checkMethod(*this, tetFem, "+=");

    dimensions_ += tetFem.dimensions_;
    BlockLduMatrix<Type>::operator+=(tetFem);
    b_ += tetFem.b_;
}


template<class Type>
void tetFemMatrix<Type>::operator+=(const tmp<tetFemMatrix<Type> >& ttetFem)
{
    operator+=(ttetFem());
    ttetFem.clear();
}


template<class Type>
void tetFemMatrix<Type>::operator+=
(
    const elemTypeGeoField& su
)
{
    checkMethod(*this, su, "+=");
    b() -= distributeField(su.internalField());
}


template<class Type>
void tetFemMatrix<Type>::operator+=
(
    const tmp<elemTypeGeoField>& tsu
)
{
    operator+=(tsu());
    tsu.clear();
}


template<class Type>
void tetFemMatrix<Type>::operator-=(const tetFemMatrix<Type>& tetFem)
{
    checkMethod(*this, tetFem, "+=");

    dimensions_ -= tetFem.dimensions_;
    BlockLduMatrix<Type>::operator-=(tetFem);
    b_ -= tetFem.b_;
}


template<class Type>
void tetFemMatrix<Type>::operator-=(const tmp<tetFemMatrix<Type> >& ttetFem)
{
    operator-=(ttetFem());
    ttetFem.clear();
}


template<class Type>
void tetFemMatrix<Type>::operator-=
(
    const elemTypeGeoField& su
)
{
    checkMethod(*this, su, "-=");
    b() += distributeField(su.internalField());
}


template<class Type>
void tetFemMatrix<Type>::operator-=
(
    const tmp<elemTypeGeoField>& tsu
)
{
    operator-=(tsu());
    tsu.clear();
}


template<class Type>
void tetFemMatrix<Type>::operator+=
(
    const dimensioned<Type>& su
)
{
    checkMethod(*this, su, "+=");
    b() -= distributeField(Field<Type>(x.mesh().nCells(), su.value()));
}


template<class Type>
void tetFemMatrix<Type>::operator-=
(
    const dimensioned<Type>& su
)
{
    checkMethod(*this, su, "-=");
    b() += distributeField(Field<Type>(x.mesh().nCells(), su.value()));
}


template<class Type>
void tetFemMatrix<Type>::operator*=
(
    const dimensioned<scalar>& ds
)
{
    dimensions_ *= ds.dimensions();
    BlockLduMatrix<Type>::operator*=(ds.value());
    b_ *= ds.value();
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
void checkMethod
(
    const tetFemMatrix<Type>& tetFem1,
    const tetFemMatrix<Type>& tetFem2,
    const char* op
)
{
    if (&tetFem1.x() != &tetFem2.x())
    {
        FatalErrorIn
        (
            "checkMethod(const tetFemMatrix<Type>&, "
            "const tetFemMatrix<Type>&) : "
        )   << "incompatible fields for operation "
            << endl << "    "
            << "[" << tetFem1.x().name() << "] "
            << op
            << " [" << tetFem1.x().name() << "]"
            << abort(FatalError);
    }

    if (dimensionSet::debug && tetFem1.dimensions() != tetFem2.dimensions())
    {
        FatalErrorIn
        (
            "checkMethod(const tetFemMatrix<Type>&, "
            "const tetFemMatrix<Type>&) : "
        )   << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << tetFem1.x().name() << tetFem1.dimensions()/dimVolume
            << " ] "
            << op
            << " [" << tetFem1.x().name() << tetFem2.dimensions()/dimVolume
            << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void checkMethod
(
    const tetFemMatrix<Type>& tetFem,
    const GeometricField<Type, elementPatchField, elementMesh>& vf,
    const char* op
)
{
    if
    (
        dimensionSet::debug
     && tetFem.dimensions() != vf.dimensions()
    )
    {
        FatalErrorIn
        (
            "checkMethod(const tetFemMatrix<Type>&, "
            "const GeometricField<Type, elementPatchField, "
            "elementMesh>&) : "
        )   << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << tetFem.x().name() << tetFem.dimensions()/dimVolume
            << " ] "
            << op
            << " [" << tetFem.x().name() << vf.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void checkMethod
(
    const tetFemMatrix<Type>& tetFem,
    const dimensioned<Type>& dt,
    const char* op
)
{
    if
    (
        dimensionSet::debug
     && tetFem.dimensions() != dt.dimensions()
    )
    {
        FatalErrorIn
        (
            "checkMethod(const tetFemMatrix<Type>&, "
            "const dimensioned<Type>&) : "
        )   << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << tetFem.x().name() << tetFem.dimensions()/dimVolume
            << " ] "
            << op
            << " [" << dt.name() << dt.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
BlockSolverPerformance<Type> solve
(
    tetFemMatrix<Type>& tetFem,
    Istream& solverControls
)
{
    return tetFem.solve(solverControls);
}


template<class Type>
BlockSolverPerformance<Type> solve
(
    const tmp<tetFemMatrix<Type> >& ttetFem,
    Istream& solverControls
)
{
    BlockSolverPerformance<Type> solverPerf =
        const_cast<tetFemMatrix<Type>&>(ttetFem()).solve(solverControls);

    ttetFem.clear();

    return solverPerf;
}


template<class Type>
BlockSolverPerformance<Type> solve(tetFemMatrix<Type>& tetFem)
{
    return tetFem.solve();
}


template<class Type>
BlockSolverPerformance<Type> solve(const tmp<tetFemMatrix<Type> >& ttetFem)
{
    BlockSolverPerformance<Type> solverPerf =
        const_cast<tetFemMatrix<Type>&>(ttetFem()).solve();

    ttetFem.clear();
    return solverPerf;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tetFemMatrix<Type>& A,
    const tetFemMatrix<Type>& B
)
{
    checkMethod(A, B, "+");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() += B;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tmp<tetFemMatrix<Type> >& tA,
    const tetFemMatrix<Type>& B
)
{
    checkMethod(tA(), B, "+");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() += B;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tetFemMatrix<Type>& A,
    const tmp<tetFemMatrix<Type> >& tB
)
{
    checkMethod(A, tB(), "+");
    tmp<tetFemMatrix<Type> > tC(tB.ptr());
    tC() += A;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tmp<tetFemMatrix<Type> >& tA,
    const tmp<tetFemMatrix<Type> >& tB
)
{
    checkMethod(tA(), tB(), "+");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() += tB();
    tB.clear();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tetFemMatrix<Type>& A
)
{
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC().negate();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<tetFemMatrix<Type> >& tA
)
{
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC().negate();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tetFemMatrix<Type>& A,
    const tetFemMatrix<Type>& B
)
{
    checkMethod(A, B, "-");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() -= B;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<tetFemMatrix<Type> >& tA,
    const tetFemMatrix<Type>& B
)
{
    checkMethod(tA(), B, "-");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() -= B;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tetFemMatrix<Type>& A,
    const tmp<tetFemMatrix<Type> >& tB
)
{
    checkMethod(A, tB(), "-");
    tmp<tetFemMatrix<Type> > tC(tB.ptr());
    tC() -= A;
    tC().negate();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<tetFemMatrix<Type> >& tA,
    const tmp<tetFemMatrix<Type> >& tB
)
{
    checkMethod(tA(), tB(), "-");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() -= tB();
    tB.clear();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tetFemMatrix<Type>& A,
    const tetFemMatrix<Type>& B
)
{
    checkMethod(A, B, "==");
    return (A - B);
}


template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tmp<tetFemMatrix<Type> >& tA,
    const tetFemMatrix<Type>& B
)
{
    checkMethod(tA(), B, "==");
    return (tA - B);
}


template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tetFemMatrix<Type>& A,
    const tmp<tetFemMatrix<Type> >& tB
)
{
    checkMethod(A, tB(), "==");
    return (A - tB);
}


template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tmp<tetFemMatrix<Type> >& tA,
    const tmp<tetFemMatrix<Type> >& tB
)
{
    checkMethod(tA(), tB(), "==");
    return (tA - tB);
}


template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tetFemMatrix<Type>& A,
    const GeometricField<Type, elementPatchField, elementMesh>& su
)
{
    checkMethod(A, su, "+");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() -= su;
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tmp<tetFemMatrix<Type> >& tA,
    const GeometricField<Type, elementPatchField, elementMesh>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() -= su;
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tetFemMatrix<Type>& A,
    const tmp<GeometricField<Type, elementPatchField, elementMesh> >& tsu
)
{
    checkMethod(A, tsu(), "+");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() -= tsu();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tmp<tetFemMatrix<Type> >& tA,
    const tmp<GeometricField<Type, elementPatchField, elementMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() -= tsu();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const GeometricField<Type, elementPatchField, elementMesh>& su,
    const tetFemMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() -= su;
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const GeometricField<Type, elementPatchField, elementMesh>& su,
    const tmp<tetFemMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() -= su;
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tmp<GeometricField<Type, elementPatchField, elementMesh> >& tsu,
    const tetFemMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "+");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() -= tsu();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tmp<GeometricField<Type, elementPatchField, elementMesh> >& tsu,
    const tmp<tetFemMatrix<Type> >& tA
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() -= tsu();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tetFemMatrix<Type>& A,
    const GeometricField<Type, elementPatchField, elementMesh>& su
)
{
    checkMethod(A, su, "-");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() += su;
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<tetFemMatrix<Type> >& tA,
    const GeometricField<Type, elementPatchField, elementMesh>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() += su;
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tetFemMatrix<Type>& A,
    const tmp<GeometricField<Type, elementPatchField, elementMesh> >& tsu
)
{
    checkMethod(A, tsu(), "-");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() += tsu();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<tetFemMatrix<Type> >& tA,
    const tmp<GeometricField<Type, elementPatchField, elementMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() += tsu();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const GeometricField<Type, elementPatchField, elementMesh>& su,
    const tetFemMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC().negate();
    tC() -= su;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const GeometricField<Type, elementPatchField, elementMesh>& su,
    const tmp<tetFemMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC() -= su;
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<GeometricField<Type, elementPatchField, elementMesh> >& tsu,
    const tetFemMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "-");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC().negate();
    tC() -= tsu();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<GeometricField<Type, elementPatchField, elementMesh> >& tsu,
    const tmp<tetFemMatrix<Type> >& tA
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC() -= tsu();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tetFemMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "+");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() -= su;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tmp<tetFemMatrix<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() -= su;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const dimensioned<Type>& su,
    const tetFemMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() -= su;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const dimensioned<Type>& su,
    const tmp<tetFemMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() -= su;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tetFemMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "-");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() += su;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<tetFemMatrix<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() += su;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const dimensioned<Type>& su,
    const tetFemMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC().negate();
    tC() -= su;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const dimensioned<Type>& su,
    const tmp<tetFemMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC() -= su;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tetFemMatrix<Type>& A,
    const GeometricField<Type, elementPatchField, elementMesh>& su
)
{
    checkMethod(A, su, "==");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() += su;
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tmp<tetFemMatrix<Type> >& tA,
    const GeometricField<Type, elementPatchField, elementMesh>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() += su;
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tetFemMatrix<Type>& A,
    const tmp<GeometricField<Type, elementPatchField, elementMesh> >& tsu
)
{
    checkMethod(A, tsu(), "==");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() += tsu();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tmp<tetFemMatrix<Type> >& tA,
    const tmp<GeometricField<Type, elementPatchField, elementMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "==");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() += tsu();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tetFemMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "==");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() += su;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tmp<tetFemMatrix<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() += su.value();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator*
(
    const dimensioned<scalar>& ds,
    const tetFemMatrix<Type>& A
)
{
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() *= ds;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator*
(
    const dimensioned<scalar>& ds,
    const tmp<tetFemMatrix<Type> >& tA
)
{
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() *= ds;
    return tC;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Ostream& operator<<(Ostream& os, const tetFemMatrix<Type>& tetFem)
{
    os  << static_cast<const BlockLduMatrix<Type>&>(tetFem) << nl
        << tetFem.dimensions_ << nl
        << tetFem.b_ << endl;

    os.check("Ostream& operator<<(Ostream&, tetFemMatrix<Type>&");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
