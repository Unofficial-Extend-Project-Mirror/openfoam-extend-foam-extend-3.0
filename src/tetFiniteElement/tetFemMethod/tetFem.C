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
    Class of static functions to calculate implicit finite element derivatives
    returning a matrix.

\*---------------------------------------------------------------------------*/

#include "tetCell.H"
#include "tetCellList.H"
#include "tetPointRef.H"
#include "SquareMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<tetFemMatrix<Type> > tetFem::laplacian
(
    GeometricField<Type, tetPolyPatchField, tetPointMesh>& vf
)
{
    elementScalarField Gamma
    (
        IOobject
        (
            "gamma",
            vf.time().constant(),
            vf.db(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        dimensionedScalar("1", dimless, 1.0)
    );

    return tetFem::laplacian(Gamma, vf);
}


template<class Type>
tmp<tetFemMatrix<Type> > tetFem::laplacian
(
    const elementScalarField& gamma,
    GeometricField<Type, tetPolyPatchField, tetPointMesh>& vf
)
{
    const tetPolyMesh& mesh = vf.mesh();

    tmp<tetFemMatrix<Type> > tfem
    (
        new tetFemMatrix<Type>
        (
            vf,
            gamma.dimensions()*vf.dimensions()/dimLength/dimLength
        )
    );
    tetFemMatrix<Type>& fem = tfem();

    // Get reference to upper and diagonal
    scalarField& u = fem.upper().asScalar();
    scalarField& d = fem.diag().asScalar();

    // In order to avoid excessive memory allocation and copying,
    // matrix coefficients and addressing are retrieved using a buffer.
    // HJ, 12/May/2001

    // Get reference to ldu addressing
    const lduAddressing& lduAddr = fem.lduAddr();
    const unallocLabelList& ownerStart = lduAddr.ownerStartAddr();
    const unallocLabelList& neighbour = lduAddr.upperAddr();

    // Create local-to-global and global-to-local addressing arrays
    labelList localToGlobalBuffer(mesh.maxNPointsForCell());
    labelList globalToLocalBuffer(lduAddr.size(), -1);

    SquareMatrix<scalar> denseMatrix
    (
        mesh.maxNPointsForCell(),
        0
    );

    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        scalar curGamma = gamma[cellI];

        label nCellPoints =
            mesh.addressing
            (
                cellI,
                localToGlobalBuffer,
                globalToLocalBuffer
            );

        mesh.gradNiDotGradNj(cellI, denseMatrix, globalToLocalBuffer);

        // Insertion and clearing of the dense matrix is done together!

        for (label localI = 0; localI < nCellPoints; localI++)
        {
            label globalI = localToGlobalBuffer[localI];

            // Insert the diagonal
            d[globalI] += curGamma*denseMatrix[localI][localI];

            // Zero out the coefficient after insertion
            denseMatrix[localI][localI] = 0;

            // Insert the off-diagonal.
            // Here, look for the known global neighbours of the point
            // using ldu owner start and grab their local index using
            // globalToLocalAddressing.  Then insert at the correct
            // place in the global matrix
            label startLabel = ownerStart[globalI];
            label endLabel = ownerStart[globalI + 1];

            for
            (
                label globalJEdge = startLabel;
                globalJEdge < endLabel;
                globalJEdge++
            )
            {
                label localJ = globalToLocalBuffer[neighbour[globalJEdge]];

                // Not all neighbours of globalI are in the dense
                // matrix; only if the global-to-local coefficient is not -1
                if (localJ > -1)
                {
                    // It is unknown if the localI or localJ is lower:
                    // that depends on how the local points have been
                    // selected.  Therefore, the matrix (which is
                    // symmetric!) needs to be "folded"
                    u[globalJEdge] +=
                        curGamma*
                        denseMatrix[min(localI, localJ)][max(localI, localJ)];

                    // Zero out the coefficient after insertion
                    denseMatrix[min(localI, localJ)][max(localI, localJ)] = 0;
                }
            }
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

    return tfem;
}


template<class Type>
tmp<tetFemMatrix<Type> > tetFem::smoother
(
    GeometricField<Type, tetPolyPatchField, tetPointMesh>& vf
)
{
    tmp<tetFemMatrix<Type> > tfem
    (
        new tetFemMatrix<Type>
        (
            vf,
            vf.dimensions()
        )
    );
    tetFemMatrix<Type>& fem = tfem();

    fem.upper() = 1.0;

    fem.negSumDiag();

    return tfem;
}


template<class Type>
tmp<tetFemMatrix<Type> > tetFem::laplacian
(
    const dimensionedScalar& gamma,
    GeometricField<Type, tetPolyPatchField, tetPointMesh>& vf
)
{
    elementScalarField Gamma
    (
        IOobject
        (
            "gamma",
            vf.instance(),
            vf.db(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        gamma
    );

    return tetFem::laplacian(Gamma, vf);
}


template<class Type>
tmp<tetFemMatrix<Type> > tetFem::laplacianTranspose
(
    GeometricField<Type, tetPolyPatchField, tetPointMesh>& vf
)
{
    elementScalarField Gamma
    (
        IOobject
        (
            "gamma",
            vf.time().constant(),
            vf.db(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        dimensionedScalar("1", dimless, 1.0)
    );

    return tetFem::laplacianTranspose(Gamma, vf);
}


template<class Type>
tmp<tetFemMatrix<Type> > tetFem::laplacianTranspose
(
    const elementScalarField& gamma,
    GeometricField<Type, tetPolyPatchField, tetPointMesh>& vf
)
{
    typedef typename CoeffField<Type>::squareType squareType;
    typedef typename CoeffField<Type>::squareTypeField squareTypeField;

    const tetPolyMesh& mesh = vf.mesh();

    tmp<tetFemMatrix<Type> > tfem
    (
        new tetFemMatrix<Type>
        (
            vf,
            gamma.dimensions()*vf.dimensions()/dimLength/dimLength
        )
    );
    tetFemMatrix<Type>& fem = tfem();

    // Get reference to upper,lower and diagonal
    squareTypeField& u = fem.upper().asSquare();
    squareTypeField& d = fem.diag().asSquare();

    // In order to avoid excessive memory allocation and copying,
    // matrix coefficients and addressing are retrieved using a buffer.
    // HJ, 12/May/2001

    // Get reference to ldu addressing
    const lduAddressing& lduAddr = fem.lduAddr();
    const unallocLabelList& ownerStart = lduAddr.ownerStartAddr();
    const unallocLabelList& neighbour = lduAddr.upperAddr();

    // Create local-to-global and global-to-local addressing arrays
    labelList localToGlobalBuffer(mesh.maxNPointsForCell());
    labelList globalToLocalBuffer(lduAddr.size(), -1);

    SquareMatrix<squareType> denseMatrix
    (
        mesh.maxNPointsForCell(),
        pTraits<squareType>::zero
    );

    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        scalar curGamma = gamma[cellI];

        label nCellPoints =
            mesh.addressing
            (
                cellI,
                localToGlobalBuffer,
                globalToLocalBuffer
            );

        mesh.gradNiGradNj(cellI, denseMatrix, globalToLocalBuffer);

        // Insertion and clearing of the dense matrix is done together!

        for (label localI = 0; localI < nCellPoints; localI++)
        {
            label globalI = localToGlobalBuffer[localI];

            // Insert the diagonal
            d[globalI] += curGamma*denseMatrix[localI][localI];

            denseMatrix[localI][localI] = pTraits<squareType>::zero;

            // Insert the off-diagonal.
            // Here, look for the known global neighbours of the point
            // using ldu owner start and grab their local index using
            // globalToLocalAddressing.  Then insert at the correct
            // place in the global matrix
            label startLabel = ownerStart[globalI];
            label endLabel = ownerStart[globalI + 1];

            for
            (
                label globalJEdge = startLabel;
                globalJEdge < endLabel;
                globalJEdge++
            )
            {
                label globalJ = neighbour[globalJEdge];
                label localJ = globalToLocalBuffer[globalJ];

                // Not all neighbours of globalI are in the dense
                // matrix; only if the global-to-local coefficient is not -1
                if (localJ > -1)
                {
                    squareType& curOffDiagCoeff =
                        denseMatrix[min(localI, localJ)][max(localI, localJ)];

                    // It is unknown if the localI or localJ is lower:
                    // that depends on how the local points have been
                    // selected.  Therefore, the matrix (which is
                    // symmetric!) needs to be "folded"
                    u[globalJEdge] += curGamma*curOffDiagCoeff.T();

                    curOffDiagCoeff = pTraits<squareType>::zero;
                }
            }
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

    return tfem;
}


template<class Type>
tmp<tetFemMatrix<Type> > tetFem::laplacianTranspose
(
    const dimensionedScalar& gamma,
    GeometricField<Type, tetPolyPatchField, tetPointMesh>& vf
)
{
    elementScalarField Gamma
    (
        IOobject
        (
            gamma.name(),
            vf.instance(),
            vf.db(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        gamma
    );

    return tetFem::laplacianTranspose(Gamma, vf);
}


template<class Type>
tmp<tetFemMatrix<Type> > tetFem::laplacianTrace
(
    const elementScalarField& gamma,
    GeometricField<Type, tetPolyPatchField, tetPointMesh>& vf
)
{
    typedef typename CoeffField<Type>::squareType squareType;
    typedef typename CoeffField<Type>::squareTypeField squareTypeField;

    const tetPolyMesh& mesh = vf.mesh();

    tmp<tetFemMatrix<Type> > tfem
    (
        new tetFemMatrix<Type>
        (
            vf,
            gamma.dimensions()*vf.dimensions()/dimLength/dimLength
        )
    );
    tetFemMatrix<Type>& fem = tfem();

    // Get reference to upper and diagonal
    squareTypeField& u = fem.upper().asSquare(); 
    squareTypeField& d = fem.diag().asSquare();

    // In order to avoid excessive memory allocation and copying,
    // matrix coefficients and addressing are retrieved using a buffer.
    // HJ, 12/May/2001

    // Get reference to ldu addressing
    const lduAddressing& lduAddr = fem.lduAddr();
    const unallocLabelList& ownerStart = lduAddr.ownerStartAddr();
    const unallocLabelList& neighbour = lduAddr.upperAddr();

    // Create local-to-global and global-to-local addressing arrays
    labelList localToGlobalBuffer(mesh.maxNPointsForCell());
    labelList globalToLocalBuffer(lduAddr.size(), -1);

    SquareMatrix<squareType> denseMatrix
    (
        mesh.maxNPointsForCell(),
        pTraits<squareType>::zero
    );

    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        scalar curGamma = gamma[cellI];

        label nCellPoints =
            mesh.addressing
            (
                cellI,
                localToGlobalBuffer,
                globalToLocalBuffer
            );

        mesh.gradNiGradNj(cellI, denseMatrix, globalToLocalBuffer);

        // Insertion and clearing of the dense matrix is done together!

        for (label localI = 0; localI < nCellPoints; localI++)
        {
            label globalI = localToGlobalBuffer[localI];

            d[globalI] += curGamma*denseMatrix[localI][localI];

            denseMatrix[localI][localI] = pTraits<squareType>::zero;

            // Insert the off-diagonal.
            // Here, look for the known global neighbours of the point
            // using ldu owner start and grab their local index using
            // globalToLocalAddressing.  Then insert at the correct
            // place in the global matrix
            label startLabel = ownerStart[globalI];
            label endLabel = ownerStart[globalI + 1];

            for
            (
                label globalJEdge = startLabel;
                globalJEdge < endLabel;
                globalJEdge++
            )
            {
                label globalJ = neighbour[globalJEdge];
                label localJ = globalToLocalBuffer[globalJ];

                // Not all neighbours of globalI are in the dense
                // matrix; only if the global-to-local coefficient is not -1
                if (localJ > -1)
                {
                    squareType& curOffDiagCoeff =
                        denseMatrix[min(localI, localJ)][max(localI, localJ)];

                    // It is unknown if the localI or localJ is lower:
                    // that depends on how the local points have been
                    // selected.  Therefore, the matrix (which is
                    // symmetric!) needs to be "folded"
                    u[globalJEdge] += curGamma*curOffDiagCoeff;

                    curOffDiagCoeff = pTraits<squareType>::zero;
                }
            }
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

    return tfem;
}


template<class Type>
tmp<tetFemMatrix<Type> > tetFem::laplacianTrace
(
    const dimensionedScalar& gamma,
    GeometricField<Type, tetPolyPatchField, tetPointMesh>& vf
)
{
    elementScalarField Gamma
    (
        IOobject
        (
            gamma.name(),
            vf.instance(),
            vf.db(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        gamma
    );

    return tetFem::laplacianTrace(Gamma, vf);
}


template<class Type>
tmp<tetFemMatrix<Type> > tetFem::laplacianTrace
(
    GeometricField<Type, tetPolyPatchField, tetPointMesh>& vf
)
{
    elementScalarField Gamma
    (
        IOobject
        (
            "1",
            vf.instance(),
            vf.db(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        dimensionedScalar("1", dimless, 1.0)
    );

    return tetFem::laplacianTrace(Gamma, vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
