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
    Tetrahedral Finite Element matrix boundary treatment tools

\*---------------------------------------------------------------------------*/

#include "tetFemMatrices.H"
#include "tetPointRef.H"
#include "constraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void tetFemMatrix<Type>::addBoundarySourceDiag()
{
    // Loop through all boundaries and fix the condition
    const FieldField<tetPolyPatchField, Type>& patches =
        x_.boundaryField();

    forAll (patches, patchI)
    {
        patches[patchI].addBoundarySourceDiag(*this);
    }
}


template<class Type>
void tetFemMatrix<Type>::setBoundaryConditions()
{
    // Make a map of all possible equations and only fix the ones that are
    // requested. Once the list of requested fixes is complete, collapse the
    // list. Note the use of list of pointers for the first part of operation
    // and creation of collapsed list by copying the pointers into a pointer
    // list.   HJ, 28/Feb/2001

    if (!boundaryConditionsSet_)
    {
        boundaryConditionsSet_ = true;

        // Loop through all boundaries and fix the condition
        const FieldField<tetPolyPatchField, Type>& patches =
            x_.boundaryField();

        this->addBoundarySourceDiag();

        forAll (patches, patchI)
        {
            patches[patchI].setBoundaryCondition(fixedEqns_);
        }

        // Loop through all fixed equations and grab the matrix
        labelList toc = fixedEqns_.toc();

        forAll (toc, eqnI)
        {
            fixedEqns_[toc[eqnI]].setMatrix(*this, x_, b_);
        }

        forAll (toc, eqnI)
        {
            fixedEqns_[toc[eqnI]].eliminateEquation(*this, b_);
        }

        forAll (toc, eqnI)
        {
            fixedEqns_[toc[eqnI]].setSourceDiag(*this, x_, b_);
        }
    }
}


template<class Type>
void tetFemMatrix<Type>::reconstructMatrix()
{
    if (!boundaryConditionsSet_)
    {
        FatalErrorIn
        (
            "void tetFemMatrix<Type>::reconstructMatrix()"
        )   << "cannot reconstruct matrix: boundary conditions not set"
            << abort(FatalError);
    }

    // Loop through all fixed equations and reconstruct the matrix
    labelList toc = fixedEqns_.toc();

    forAll (toc, eqnI)
    {
        fixedEqns_[toc[eqnI]].reconstructMatrix(*this);
    }

    // Matrix is restored
    boundaryConditionsSet_ = true;
}


template<class Type>
void tetFemMatrix<Type>::addCouplingCoeffs()
{
    // Diagonal

    if (this->hasDiag())
    {
        // Initialise diagonal transfer
        forAll (x_.boundaryField(), patchI)
        {
            const tetPolyPatchField<Type>& ptf = x_.boundaryField()[patchI];

            if (ptf.coupled())
            {
                ptf.initAddDiag(this->diag());
            }
        }

        // Add the diagonal
        forAll (x_.boundaryField(), patchI)
        {
            const tetPolyPatchField<Type>& ptf = x_.boundaryField()[patchI];

            if (ptf.coupled())
            {
                ptf.addDiag(this->diag());
            }
        }
    }

    // Upper triangle

    if (this->hasUpper())
    {
        // Initialise upper transfer
        forAll (x_.boundaryField(), patchI)
        {
            const tetPolyPatchField<Type>& ptf = x_.boundaryField()[patchI];

            if (ptf.coupled())
            {
                ptf.initAddUpperLower(this->upper());
            }
        }

        // Add the upper
        forAll (x_.boundaryField(), patchI)
        {
            const tetPolyPatchField<Type>& ptf = x_.boundaryField()[patchI];

            if (ptf.coupled())
            {
                ptf.addUpperLower(this->upper());
            }
        }
    }

    // Lower triangle

    if (this->hasLower())
    {
        // Initialise lower transfer
        forAll (x_.boundaryField(), patchI)
        {
            const tetPolyPatchField<Type>& ptf = x_.boundaryField()[patchI];

            if (ptf.coupled())
            {
                ptf.initAddUpperLower(this->lower());
            }
        }

        // Add the lower
        forAll (x_.boundaryField(), patchI)
        {
            const tetPolyPatchField<Type>& ptf = x_.boundaryField()[patchI];

            if (ptf.coupled())
            {
                ptf.addUpperLower(this->lower());
            }
        }
    }
}


template<class Type>
tmp<Field<Type> > tetFemMatrix<Type>::distributeSource
(
    const Field<Type>& ef
) const
{
    const tetPolyMesh& mesh = x_.mesh();

    // Make a field over all points
    tmp<Field<Type> > tdistSu
    (
        new Field<Type>
        (
            mesh.nPoints(),
            pTraits<Type>::zero
        )
    );
    Field<Type>& distSu = tdistSu();

    // Go through all the tets and add a quarter of the volume times ef
    // into each of the vertices. The cell integrals are prepared by the mesh
    labelList localToGlobalBuffer(mesh.maxNPointsForCell());
    labelList globalToLocalBuffer(this->lduAddr().size(), -1);

    scalarField coeffsBuffer
    (
        mesh.maxNPointsForCell()*tetPointRef::nVertices
    );

    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        label nCellPoints =
            mesh.addressing
            (
                cellI,
                localToGlobalBuffer,
                globalToLocalBuffer
            );

        mesh.volIntegral(cellI, coeffsBuffer, globalToLocalBuffer);

        const Type& curEf = ef[cellI];

        for (label localI = 0; localI < nCellPoints; localI++)
        {
            label globalI = localToGlobalBuffer[localI];

            // Insert the source
            distSu[globalI] += coeffsBuffer[localI]*curEf;
        }

        // Clear addressing for element
        mesh.clearAddressing
        (
            cellI,
            nCellPoints,
            localToGlobalBuffer,
            globalToLocalBuffer
        );
    }

    return tdistSu;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
