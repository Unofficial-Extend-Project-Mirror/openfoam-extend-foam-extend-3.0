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
    Tetrahedral Finite Element matrix basic solvers.

\*---------------------------------------------------------------------------*/

#include "blockLduSolvers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Solvers * * * * * * * * * * * * * * * * * //

template<class Type>
BlockSolverPerformance<Type> tetFemMatrix<Type>::solve
(
     const dictionary& solverControls
)
{
    if (debug)
    {
        InfoIn("tetFemMatrix<Type>::solve(const dictionary&)")
            << "solving tetFemMatrix<Type>" << endl;
    }

    // Set boundary conditions
    this->setBoundaryConditions();

    BlockSolverPerformance<Type> solverPerf =
        BlockLduSolver<Type>::New
        (
            x_.name(),
            *this,
            solverControls
        )->solve(x_, b_);
        
    solverPerf.print();

    this->reconstructMatrix();

    if (debug)
    {
        Info<< "tetFemMatrix<Type>::solve : correcting boundary conditions"
            << endl;
    }

    x_.correctBoundaryConditions();

    return solverPerf;
}


template<class Type>
BlockSolverPerformance<Type> tetFemMatrix<Type>::solve()
{
    return solve(x_.mesh().solutionDict().solver(x_.name()));
}


// Return the matrix residual
template<class Type>
tmp<Field<Type> > tetFemMatrix<Type>::residual()
{
    // Set boundary conditions
    this->setBoundaryConditions();

    tmp<Field<Type> > tres =
        BlockLduMatrix<Type>::residual(x_, b_);

    this->reconstructMatrix();

    return tres;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
