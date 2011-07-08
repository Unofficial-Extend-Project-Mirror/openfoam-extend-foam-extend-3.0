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

\*---------------------------------------------------------------------------*/

#include "blockFvScalarMatrix.H"
#include "scalarCoeffField.H"
#include "scalarBlockLduMatrix.H"
// #include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// template<>
// void Foam::fvMatrix<Foam::scalar>::setComponentReference
// (
//     const label patchi,
//     const label facei,
//     const direction,
//     const scalar value
// )
// {
//     if (psi_.needReference())
//     {
//         if (Pstream::master())
//         {
//             internalCoeffs_[patchi][facei] +=
//                 diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]];
// 
//             boundaryCoeffs_[patchi][facei] +=
//                 diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]]
//                *value;
//         }
//     }
// }
// 
// 
// template<>
// Foam::autoPtr<Foam::fvMatrix<Foam::scalar>::fvSolver>
// Foam::fvMatrix<Foam::scalar>::solver
// (
//     const dictionary& solverControls
// )
// {
//     if (debug)
//     {
//         Info<< "fvMatrix<scalar>::solver(const dictionary& solverControls) : "
//                "solver for fvMatrix<scalar>"
//             << endl;
//     }
// 
//     scalarField saveDiag = diag();
//     addBoundaryDiag(diag(), 0);
// 
//     // Make a copy of interfaces: no longer a reference
//     // HJ, 20/Nov/2007
//     lduInterfaceFieldPtrsList interfaces = psi_.boundaryField().interfaces();
// 
//     autoPtr<fvMatrix<scalar>::fvSolver> solverPtr
//     (
//         new fvMatrix<scalar>::fvSolver
//         (
//             *this,
//             lduSolver::New
//             (
//                 psi_.name(),
//                 *this,
//                 boundaryCoeffs_,
//                 internalCoeffs_,
//                 interfaces,
//                 solverControls
//             )
//         )
//     );
// 
//     diag() = saveDiag;
// 
//     return solverPtr;
// }
// 
// 
// template<>
// Foam::lduMatrix::solverPerformance
// Foam::fvMatrix<Foam::scalar>::fvSolver::solve
// (
//     const dictionary& solverControls
// )
// {
//     scalarField saveDiag = fvMat_.diag();
//     fvMat_.addBoundaryDiag(fvMat_.diag(), 0);
// 
//     scalarField totalSource = fvMat_.source();
//     fvMat_.addBoundarySource(totalSource, false);
// 
//     // assign new solver controls
//     solver_->read(solverControls);
//     lduSolverPerformance solverPerf =
//         solver_->solve(fvMat_.psi().internalField(), totalSource);
// 
//     solverPerf.print();
// 
//     fvMat_.diag() = saveDiag;
// 
//     fvMat_.psi().correctBoundaryConditions();
// 
//     return solverPerf;
// }
// 
// 
// template<>
// Foam::lduMatrix::solverPerformance Foam::fvMatrix<Foam::scalar>::solve
// (
//     const dictionary& solverControls
// )
// {
//     if (debug)
//     {
//         Info<< "fvMatrix<scalar>::solve(const dictionary& solverControls) : "
//                "solving fvMatrix<scalar>"
//             << endl;
//     }
// 
//     scalarField saveDiag = diag();
//     addBoundaryDiag(diag(), 0);
// 
//     scalarField totalSource = source_;
//     addBoundarySource(totalSource, false);
// 
//     // Make a copy of interfaces: no longer a reference
//     // HJ, 20/Nov/2007
//     lduInterfaceFieldPtrsList interfaces = psi_.boundaryField().interfaces();
// 
//     // Solver call
//     lduSolverPerformance solverPerf = lduSolver::New
//     (
//         psi_.name(),
//         *this,
//         boundaryCoeffs_,
//         internalCoeffs_,
//         interfaces,
//         solverControls
//     )->solve(psi_.internalField(), totalSource);
// 
//     solverPerf.print();
// 
//     diag() = saveDiag;
// 
//     psi_.correctBoundaryConditions();
// 
//     return solverPerf;
// }
// 
// 
// template<>
// Foam::tmp<Foam::scalarField> Foam::fvMatrix<Foam::scalar>::residual() const
// {
//     scalarField boundaryDiag(psi_.size(), 0.0);
//     addBoundaryDiag(boundaryDiag, 0);
// 
//     // Make a copy of interfaces: no longer a reference
//     // HJ, 20/Nov/2007
//     lduInterfaceFieldPtrsList interfaces = psi_.boundaryField().interfaces();
// 
//     tmp<scalarField> tres
//     (
//         lduMatrix::residual
//         (
//             psi_.internalField(),
//             source_ - boundaryDiag*psi_.internalField(),
//             boundaryCoeffs_,
//             interfaces,
//             0
//         )
//     );
// 
//     addBoundarySource(tres());
// 
//     return tres;
// }
// 
// 
// template<>
// Foam::tmp<Foam::volScalarField> Foam::fvMatrix<Foam::scalar>::H() const
// {
//     tmp<volScalarField> tHphi
//     (
//         new volScalarField
//         (
//             IOobject
//             (
//                 "H("+psi_.name()+')',
//                 psi_.instance(),
//                 psi_.mesh(),
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             psi_.mesh(),
//             dimensions_/dimVol,
//             zeroGradientFvPatchScalarField::typeName
//         )
//     );
//     volScalarField& Hphi = tHphi();
// 
//     Hphi.internalField() = (lduMatrix::H(psi_.internalField()) + source_);
//     addBoundarySource(Hphi.internalField());
// 
//     Hphi.internalField() /= psi_.mesh().V();
//     Hphi.correctBoundaryConditions();
// 
//     return tHphi;
// }
// 
// 
// template<>
// Foam::tmp<Foam::volScalarField> Foam::fvMatrix<Foam::scalar>::H1() const
// {
//     tmp<volScalarField> tH1
//     (
//         new volScalarField
//         (
//             IOobject
//             (
//                 "H(1)",
//                 psi_.instance(),
//                 psi_.mesh(),
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             psi_.mesh(),
//             dimensions_/(dimVol*psi_.dimensions()),
//             zeroGradientFvPatchScalarField::typeName
//         )
//     );
//     volScalarField& H1_ = tH1();
// 
//     H1_.internalField() = lduMatrix::H1();
//     //addBoundarySource(Hphi.internalField());
// 
//     H1_.internalField() /= psi_.mesh().V();
//     H1_.correctBoundaryConditions();
// 
//     return tH1;
// }

template<>
void Foam::blockFvMatrix<Foam::scalar>::addBoundarySource
(
    Field<scalar>& source,
    const bool couples
) const
{
    forAll(psi_.boundaryField(), patchI)
    {
        const fvPatchField<scalar>& ptf = psi_.boundaryField()[patchI];

        const Field<scalar>& interfaceSource
                = BlockLduMatrix<scalar>::interfaceSource()[patchI];

        addToInternalField
        (
            BlockLduMatrix<scalar>::lduAddr().patchAddr(patchI),
            interfaceSource, 
            source
        );
        
        if (ptf.coupled() && couples)
        {
            const CoeffField<scalar>& upperCoeffs
                    = BlockLduMatrix<scalar>::coupleUpper()[patchI];

            if(upperCoeffs.activeType() != blockCoeffBase::UNALLOCATED)
            {
                Field<scalar> pnf = upperCoeffs.asScalar()*ptf.patchNeighbourField();
                
                const unallocLabelList& addr = BlockLduMatrix<scalar>::lduAddr().patchAddr(patchI);

                forAll(addr, facei)
                {
                    source[addr[facei]] += pnf[facei];
                }
            }
        }
    }
}


template<>
void Foam::blockFvMatrix<Foam::scalar>::setReference
(
    const label celli,
    const scalar& value,
    const bool forceReference
)
{
    if (celli >= 0 && (psi_.needReference() || forceReference))
    {
        // Bug fix: force reference only on master for parallel runs
        // HJ, 12/Feb/2010
        if (Pstream::parRun())
        {
            // Parallel run:
            // - only set reference on master processor: one place is enough
            // - make sure that cellI is not out of range
            if (Pstream::master())
            {
                label parCelli = celli;

                while (parCelli >= BlockLduMatrix<scalar>::diag().size())
                {
                    // Out of range, pick a local cell
                    parCelli /= Pstream::nProcs();
                }

                scalar& diagCoeff = BlockLduMatrix<scalar>::diag()[parCelli];
                
                BlockLduMatrix<scalar>::source()[parCelli] += diagCoeff*value;
                
                diagCoeff += diagCoeff;
            }
        }
        else
        {
            // Serial run, standard practice
                scalar& diagCoeff = BlockLduMatrix<scalar>::diag()[celli];
                
                BlockLduMatrix<scalar>::source()[celli] += diagCoeff*value;
                
                diagCoeff += diagCoeff;

        }
    }
}


template<>
Foam::tmp<Foam::scalarField > Foam::blockFvMatrix<Foam::scalar>::D() const
{
    CoeffField<scalar> diag = blockD();
    
    return tmp<scalarField>(new scalarField(diag.asScalar()));
}


template<>
Foam::tmp<Foam::scalarField> Foam::blockFvMatrix<Foam::scalar>::DD() const
{
    tmp<scalarField> tdiag = D();
    scalarField& diag = tdiag();

    forAll(psi_.boundaryField(), patchI)
    {
        const fvPatchField<scalar>& ptf = psi_.boundaryField()[patchI];

        if (!ptf.coupled() && ptf.size())
        {
            addToInternalField
            (
                BlockLduMatrix<scalar>::lduAddr().patchAddr(patchI),
                BlockLduMatrix<scalar>::interfaceDiag()[patchI],
                diag
            );
        }
    }
    
    return tdiag;
}


template<>
Foam::tmp<Foam::volScalarField> Foam::blockFvMatrix<Foam::scalar>::H() const
{
    tmp<volScalarField> tHphi
    (
        new volScalarField
        (
            IOobject
            (
                "H("+psi_.name()+')',
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions_/dimVol,
            zeroGradientFvPatchField<scalar>::typeName
        )
    );
    volScalarField& Hphi = tHphi();
    
    Hphi.internalField() = BlockLduMatrix<scalar>::H(psi_.internalField());    
    Hphi.internalField() += BlockLduMatrix<scalar>::source();
    addBoundarySource(Hphi.internalField());
    
    Hphi.internalField() /= psi_.mesh().V();
    Hphi.correctBoundaryConditions();

    return tHphi;
}


template<>
Foam::tmp<Foam::surfaceScalarField>
Foam::blockFvMatrix<Foam::scalar>::flux() const
{
    if (!psi_.mesh().fluxRequired(psi_.name()))
    {
        FatalErrorIn("blockFvMatrix<Type>::flux()")
            << "flux requested but " << psi_.name()
            << " not specified in the fluxRequired sub-dictionary"
               " of fvSchemes."
            << abort(FatalError);
    }

    // construct GeometricField<Type, fvsPatchField, surfaceMesh>
    tmp<surfaceScalarField> tfieldFlux
    (
        new surfaceScalarField
        (
            IOobject
            (
                "flux("+psi_.name()+')',
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions()
        )
    );
    surfaceScalarField& fieldFlux = tfieldFlux();

    fieldFlux.internalField() = BlockLduMatrix<scalar>::faceH
    (
        psi_.internalField()
    );

    //- Boundary diagonal contribution
    forAll(psi_.boundaryField(), patchI)
    {
        fieldFlux.boundaryField()[patchI] =
        (
            BlockLduMatrix<scalar>::interfaceDiag()[patchI].asScalar()
          * psi_.boundaryField()[patchI].patchInternalField()
        );
    }
        
    //- Boundary diagonal contribution
    forAll(psi_.boundaryField(), patchI)
    {        
        //- Coupled Boundaries off-diagonal contribution
        if (psi_.boundaryField()[patchI].coupled())
        {
            fieldFlux.boundaryField()[patchI] -=
            (
                BlockLduMatrix<scalar>::coupleUpper()[patchI].asScalar()
              * psi_.boundaryField()[patchI].patchNeighbourField()
            );
        }
    }

    if (faceFluxCorrectionPtr_)
    {
        fieldFlux += *faceFluxCorrectionPtr_;
    }

    return tfieldFlux;
}


// ************************************************************************* //
