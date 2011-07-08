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
    
\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "surfaceFields.H"
#include "calculatedFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "coupledFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template<class Type2>
void blockFvMatrix<Type>::addToInternalField
(
    const unallocLabelList& addr,
    const Field<Type2>& pf,
    Field<Type2>& intf
) const
{
    if (addr.size() != pf.size())
    {
        FatalErrorIn
        (
            "blockFvMatrix<Type>::addToInternalField(const unallocLabelList&, "
            "const Field&, Field&)"
        )   << "sizes of addressing and field are different"
            << abort(FatalError);
    }

    forAll(addr, faceI)
    {
        intf[addr[faceI]] += pf[faceI];
    }
}


template<class Type>
template<class Type2>
void blockFvMatrix<Type>::addToInternalField
(
    const unallocLabelList& addr,
    const tmp<Field<Type2> >& tpf,
    Field<Type2>& intf
) const
{
    addToInternalField(addr, tpf(), intf);
    tpf.clear();
}


template<class Type>
template<class Type2>
void blockFvMatrix<Type>::subtractFromInternalField
(
    const unallocLabelList& addr,
    const Field<Type2>& pf,
    Field<Type2>& intf
) const
{
    if (addr.size() != pf.size())
    {
        FatalErrorIn
        (
            "blockFvMatrix<Type>::addToInternalField(const unallocLabelList&, "
            "const Field&, Field&)"
        )   << "sizes of addressing and field are different"
            << abort(FatalError);
    }

    forAll(addr, faceI)
    {
        intf[addr[faceI]] -= pf[faceI];
    }
}


template<class Type>
template<class Type2>
void blockFvMatrix<Type>::subtractFromInternalField
(
    const unallocLabelList& addr,
    const tmp<Field<Type2> >& tpf,
    Field<Type2>& intf
) const
{
    subtractFromInternalField(addr, tpf(), intf);
    tpf.clear();
}


template<class Type>
void blockFvMatrix<Type>::addBoundaryDiag
(
    CoeffField<Type>& diag
) const
{
    forAll(BlockLduMatrix<Type>::interfaceDiag(), patchI)
    {
        diag.addSubset
        (
            BlockLduMatrix<Type>::interfaceDiag()[patchI],
            BlockLduMatrix<Type>::lduAddr().patchAddr(patchI)
        );
    }
}


template<class Type>
void blockFvMatrix<Type>::addBoundarySource
(
    Field<Type>& source,
    const bool couples
) const
{
    forAll(psi_.boundaryField(), patchI)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        const Field<Type>& interfaceSource
                = BlockLduMatrix<Type>::interfaceSource()[patchI];

        addToInternalField
        (
            BlockLduMatrix<Type>::lduAddr().patchAddr(patchI),
            interfaceSource, 
            source
        );
        
        if (ptf.coupled() && couples)
        {
            const CoeffField<Type>& upperCoeffs
                    = BlockLduMatrix<Type>::coupleUpper()[patchI];

            if(upperCoeffs.activeType() != blockCoeffBase::UNALLOCATED)
            {
                tmp<Field<Type> > tpnf = ptf.patchNeighbourField();
                Field<Type>& pnf = tpnf();

                if (upperCoeffs.activeType() == blockCoeffBase::SCALAR)
                {
                    pnf = upperCoeffs.asScalar()*pnf;
                }
                else if (upperCoeffs.activeType() == blockCoeffBase::LINEAR)
                {
                    pnf = cmptMultiply(upperCoeffs.asLinear(), pnf);
                }
                else if (upperCoeffs.activeType() == blockCoeffBase::SQUARE)
                {
                    pnf = upperCoeffs.asSquare() & pnf;
                }
                
                const unallocLabelList& addr = BlockLduMatrix<Type>::lduAddr().patchAddr(patchI);

                forAll(addr, facei)
                {
                    source[addr[facei]] += pnf[facei];
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
blockFvMatrix<Type>::blockFvMatrix
(
    GeometricField<Type, fvPatchField, volMesh>& psi,
    const dimensionSet& ds
)
:
    BlockLduMatrix<Type>(psi.mesh()),
    psi_(psi),
    dimensions_(ds),
    faceFluxCorrectionPtr_(NULL)
{
    if (debug)
    {
        Info<< "blockFvMatrix<Type>::blockFvMatrix(GeometricField<Type, fvPatchField, volMesh>&,"
               " const dimensionSet&) : "
               "constructing blockFvMatrix<Type> for field " << psi_.name()
            << endl;
    }

    // Initialise interfaces
    BlockLduMatrix<Type>::interfaces() = psi.boundaryField().blockInterfaces();
    
    psi_.boundaryField().updateCoeffs();
}


template<class Type>
blockFvMatrix<Type>::blockFvMatrix(const blockFvMatrix<Type>& fvm)
:
    BlockLduMatrix<Type>(fvm),
    psi_(fvm.psi_),
    dimensions_(fvm.dimensions_),
    faceFluxCorrectionPtr_(NULL)
{
    if (debug)
    {
        Info<< "blockFvMatrix<Type>::blockFvMatrix(const blockFvMatrix<Type>&) : "
            << "copying blockFvMatrix<Type> for field " << psi_.name()
            << endl;
    }

    if (fvm.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = new
        GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            *(fvm.faceFluxCorrectionPtr_)
        );
    }
}


#ifdef ConstructFromTmp
template<class Type>
blockFvMatrix<Type>::blockFvMatrix(const tmp<blockFvMatrix<Type> >& tfvm)
:
    BlockLduMatrix<Type>
    (
        const_cast<blockFvMatrix<Type>&>(tfvm()),
        tfvm.isTmp()
    ),
    psi_(tfvm().psi_),
    dimensions_(tfvm().dimensions_),
    faceFluxCorrectionPtr_(NULL)
{
    if (debug)
    {
        Info<< "blockFvMatrix<Type>::blockFvMatrix(const tmp<blockFvMatrix<Type> >&) : "
            << "copying blockFvMatrix<Type> for field " << psi_.name()
            << endl;
    }

    if (tfvm().faceFluxCorrectionPtr_)
    {
        if (tfvm.isTmp())
        {
            faceFluxCorrectionPtr_ = tfvm().faceFluxCorrectionPtr_;
            tfvm().faceFluxCorrectionPtr_ = NULL;
        }
        else
        {
            faceFluxCorrectionPtr_ = new
                GeometricField<Type, fvsPatchField, surfaceMesh>
                (
                    *(tfvm().faceFluxCorrectionPtr_)
                );
        }
    }

    tfvm.clear();
}
#endif


template<class Type>
blockFvMatrix<Type>::blockFvMatrix
(
    GeometricField<Type, fvPatchField, volMesh>& psi,
    Istream& is
)
:
    BlockLduMatrix<Type>(psi.mesh()),
    psi_(psi),
    dimensions_(is),
    faceFluxCorrectionPtr_(NULL)
{
    if (debug)
    {
        Info<< "blockFvMatrix<Type>::blockFvMatrix(GeometricField<Type, fvPatchField, volMesh>&,"
               " Istream&) : "
               "constructing blockFvMatrix<Type> for field " << psi_.name()
            << endl;
    }

    // Initialise source
    BlockLduMatrix<Type>::source() << is;
    
    // Initialise interfaces
    BlockLduMatrix<Type>::interfaces() = psi.boundaryField().blockInterfaces();
}


template<class Type>
blockFvMatrix<Type>::~blockFvMatrix()
{
    if (debug)
    {
        Info<< "blockFvMatrix<Type>::~blockFvMatrix<Type>() : "
            << "destroying blockFvMatrix<Type> for field " << psi_.name()
            << endl;
    }

    if (faceFluxCorrectionPtr_)
    {
        delete faceFluxCorrectionPtr_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
void blockFvMatrix<Type>::setReference
(
    const label celli,
    const Type& value,
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

                while (parCelli >= BlockLduMatrix<Type>::diag().size())
                {
                    // Out of range, pick a local cell
                    parCelli /= Pstream::nProcs();
                }

                BlockCoeff<Type>& diagCoeff = BlockLduMatrix<Type>::diag()[parCelli];
                
                BlockLduMatrix<Type>::source()[parCelli]
                    += BlockCoeff<Type>::multiply(diagCoeff, value);
                    
                diagCoeff += diagCoeff;
            }
        }
        else
        {
            // Serial run, standard practice
            BlockCoeff<Type>& diagCoeff = BlockLduMatrix<Type>::diag()[celli];
            
            BlockLduMatrix<Type>::source()[celli]
                += BlockCoeff<Type>::multiply(diagCoeff, value);
                
            diagCoeff += diagCoeff;
        }
    }
}


template<class Type>
tmp<CoeffField<Type> > blockFvMatrix<Type>::blockD() const
{
    tmp<CoeffField<Type> > tdiag
    (
        new CoeffField<Type>(BlockLduMatrix<Type>::diag())
    );
    
    addBoundaryDiag(tdiag());
    
    return tdiag;
}


template<class Type>
tmp<scalarField> blockFvMatrix<Type>::D() const
{
    const CoeffField<Type> diag(blockD());
    
    tmp<scalarField> tscalarDiag(new scalarField(diag.size()));
    scalarField& scalarDiag = tscalarDiag();
    
    if (diag.activeType() == blockCoeffBase::SCALAR)
    {
        scalarDiag = diag.asScalar();
    }
    else if (diag.activeType() == blockCoeffBase::LINEAR)
    {
        cmptAv(scalarDiag, diag.asLinear());
    }
    else if (diag.activeType() == blockCoeffBase::SQUARE)
    {
        contractScalar(scalarDiag, diag.asSquare());
    }
    
    return tscalarDiag;
}


template<class Type>
tmp<Field<Type> > blockFvMatrix<Type>::DD() const
{
    CoeffField<Type> diag(blockD());

    forAll(psi_.boundaryField(), patchI)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (!ptf.coupled() && ptf.size())
        {
            addToInternalField
            (
                BlockLduMatrix<Type>::lduAddr().patchAddr(patchI),
                BlockLduMatrix<Type>::interfaceDiag()[patchI],
                diag
            );
        }
    }

    tmp<Field<Type> > ttypeDiag(new Field<Type>());
    Field<Type>& typeDiag = ttypeDiag();

    if (diag.activeType() == blockCoeffBase::SCALAR)
    {
        expandLinear(typeDiag, diag.asScalar());
    }
    else if (diag.activeType() == blockCoeffBase::LINEAR)
    {
        typeDiag = diag.asLinear();
    }
    else if (diag.activeType() == blockCoeffBase::SQUARE)
    {
        contractLinear(typeDiag, diag.asSquare());
    }

    return ttypeDiag;
}


template<class Type>
tmp<volScalarField> blockFvMatrix<Type>::A() const
{
    tmp<volScalarField> tAphi
    (
        new volScalarField
        (
            IOobject
            (
                "A("+psi_.name()+')',
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions_/psi_.dimensions()/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tAphi().internalField() = D()/psi_.mesh().V();
    tAphi().correctBoundaryConditions();

    return tAphi;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > blockFvMatrix<Type>::H() const
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tHphi
    (
        new GeometricField<Type, fvPatchField, volMesh>
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
            zeroGradientFvPatchField<Type>::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& Hphi = tHphi();
    
    Hphi.internalField() = BlockLduMatrix<Type>::H(psi_.internalField());    
    Hphi.internalField() += BlockLduMatrix<Type>::source();
    addBoundarySource(Hphi.internalField());

    // Ivor Clifford 08/07/2011
    // NOTE - We only need to do this since blockFvMatrix::A() returns
    // SCALAR rather than LINEAR or full SQUARE values and we want to
    // remain compatible with the old fvMatrix approach.
    // I don't like this, we should be fixing blockFvMatrix::A() instead.
    
    // Correct for diagonal coefficient off-diagonal components.
    const Field<Type>& psii = psi_.internalField();
    
    Field<Type>& Hphii = Hphi.internalField();
    
    typedef typename BlockCoeff<Type>::squareTypeField squareTypeField;

    if (BlockLduMatrix<Type>::diag().activeType() == blockCoeffBase::LINEAR)
    {
        const Field<Type>& diag_ = BlockLduMatrix<Type>::diag().asLinear();
        
        forAll(psii, celli)
        {
            Hphii[celli] += cmptAv(diag_[celli])*psii[celli]
                - cmptMultiply(diag_[celli], psii[celli]);
        }
    }
    if (BlockLduMatrix<Type>::diag().activeType() == blockCoeffBase::SQUARE)
    {
        const squareTypeField& diag_ = BlockLduMatrix<Type>::diag().asSquare();
            
        forAll(psii, celli)
        {
            Hphii[celli] += contractScalar(diag_[celli])*psii[celli]
                - (diag_[celli] & psii[celli]);
        }
    }
    
    // Correct for boundary coefficient off-diagonal components.
    forAll(BlockLduMatrix<Type>::interfaceDiag(), patchI)
    {
        const CoeffField<Type>& interfaceDiag_
            = BlockLduMatrix<Type>::interfaceDiag()[patchI];
            
        const unallocLabelList& addr
            = BlockLduMatrix<Type>::lduAddr().patchAddr(patchI);
                
        if(interfaceDiag_.activeType() == blockCoeffBase::LINEAR)
        {
            const Field<Type>& coeffs = interfaceDiag_.asLinear();
            
            forAll(addr, facei)
            {
                label celli = addr[facei];
                
                Hphii[celli] += cmptAv(coeffs[facei])*psii[celli]
                    - cmptMultiply(coeffs[facei], psii[celli]);
            }
        }
        else if(interfaceDiag_.activeType() == blockCoeffBase::SQUARE)
        {
            const squareTypeField& coeffs = interfaceDiag_.asSquare();
            
            forAll(addr, facei)
            {
                label celli = addr[facei];
                
                Hphii[celli] += contractScalar(coeffs[facei])*psii[celli]
                    - (coeffs[facei] & psii[celli]);
            }
        }
    }
    
    Hphi.internalField() /= psi_.mesh().V();
    Hphi.correctBoundaryConditions();

    return tHphi;
}


template<class Type>
tmp<volScalarField> blockFvMatrix<Type>::H1() const
{
    tmp<volScalarField> tH1
    (
        new volScalarField
        (
            IOobject
            (
                "H(1)",
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions_/(dimVol*psi_.dimensions()),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& H1_ = tH1();

    H1_.internalField() = BlockLduMatrix<Type>::H1();

    H1_.internalField() /= psi_.mesh().V();
    H1_.correctBoundaryConditions();

    return tH1;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
blockFvMatrix<Type>::flux() const
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
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tfieldFlux
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
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
    GeometricField<Type, fvsPatchField, surfaceMesh>& fieldFlux = tfieldFlux();

    fieldFlux.internalField() = BlockLduMatrix<Type>::faceH(psi_.internalField());

    //- Boundary diagonal contribution
    forAll(psi_.boundaryField(), patchI)
    {
        Field<Type>& pFieldFlux = fieldFlux.boundaryField()[patchI];
        
        const CoeffField<Type>& pInterfaceDiag
            = BlockLduMatrix<scalar>::interfaceDiag()[patchI];
        
        if (pInterfaceDiag.activeType() == blockCoeffBase::SCALAR)
        {
            pFieldFlux = pInterfaceDiag.asScalar()
                * psi_.boundaryField()[patchI].patchInternalField();
        }
        else if (pInterfaceDiag.activeType() == blockCoeffBase::LINEAR)
        {
            pFieldFlux = cmptMultiply
            (
                pInterfaceDiag.asLinear(),
                psi_.boundaryField()[patchI].patchInternalField()
            );
        }
        else if (pInterfaceDiag.activeType() == blockCoeffBase::SQUARE)
        {
            pFieldFlux =
            (
                pInterfaceDiag.asSquare()
              & psi_.boundaryField()[patchI].patchInternalField()
            );
        }
        else
        {
            pFieldFlux = pTraits<Type>::zero;
        }
    }
        
    //- Boundary diagonal contribution
    forAll(psi_.boundaryField(), patchI)
    {        
        //- Coupled Boundaries off-diagonal contribution
        if (psi_.boundaryField()[patchI].coupled())
        {
            Field<Type>& pFieldFlux = fieldFlux.boundaryField()[patchI];
            
            CoeffField<Type>& pCoupleUpper
                = BlockLduMatrix<scalar>::coupleUpper()[patchI];
            
            if (pCoupleUpper.activeType() == blockCoeffBase::SCALAR)
            {
                pFieldFlux -= pCoupleUpper.asScalar()
                    * psi_.boundaryField()[patchI].patchNeighbourField();
            }
            else if (pCoupleUpper.activeType() == blockCoeffBase::LINEAR)
            {
                pFieldFlux -= cmptMultiply
                (
                    pCoupleUpper.asLinear(),
                    psi_.boundaryField()[patchI].patchNeighbourField()
                );
            }
            else if (pCoupleUpper.activeType() == blockCoeffBase::SQUARE)
            {
                pFieldFlux -=
                (
                    pCoupleUpper.asSquare()
                  & psi_.boundaryField()[patchI].patchNeighbourField()
                );
            }
        }
    }

    if (faceFluxCorrectionPtr_)
    {
        fieldFlux += *faceFluxCorrectionPtr_;
    }

    return tfieldFlux;
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void blockFvMatrix<Type>::operator=(const blockFvMatrix<Type>& fvmv)
{
    if (this == &fvmv)
    {
        FatalErrorIn("blockFvMatrix<Type>::operator=(const blockFvMatrix<Type>&)")
            << "attempted assignment to self"
            << abort(FatalError);
    }

    if (&psi_ != &(fvmv.psi_))
    {
        FatalErrorIn("blockFvMatrix<Type>::operator=(const blockFvMatrix<Type>&)")
            << "different fields"
            << abort(FatalError);
    }

    BlockLduMatrix<Type>::operator=(fvmv);

    if (faceFluxCorrectionPtr_ && fvmv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ = *fvmv.faceFluxCorrectionPtr_;
    }
    else if (fvmv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ =
            new GeometricField<Type, fvsPatchField, surfaceMesh>
        (*fvmv.faceFluxCorrectionPtr_);
    }
}


template<class Type>
void blockFvMatrix<Type>::operator=(const tmp<blockFvMatrix<Type> >& tfvmv)
{
    operator=(tfvmv());
    tfvmv.clear();
}


template<class Type>
void blockFvMatrix<Type>::negate()
{
    BlockLduMatrix<Type>::negate();

    if (faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_->negate();
    }
}


template<class Type>
void blockFvMatrix<Type>::operator+=(const blockFvMatrix<Type>& fvmv)
{
    checkMethod(*this, fvmv, "+=");

    dimensions_ += fvmv.dimensions_;
    BlockLduMatrix<Type>::operator+=(fvmv);
    
    if (faceFluxCorrectionPtr_ && fvmv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ += *fvmv.faceFluxCorrectionPtr_;
    }
    else if (fvmv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = new
        GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            *fvmv.faceFluxCorrectionPtr_
        );
    }
}


template<class Type>
void blockFvMatrix<Type>::operator+=(const tmp<blockFvMatrix<Type> >& tfvmv)
{
    operator+=(tfvmv());
    tfvmv.clear();
}


template<class Type>
void blockFvMatrix<Type>::operator-=(const blockFvMatrix<Type>& fvmv)
{
    checkMethod(*this, fvmv, "+=");

    dimensions_ -= fvmv.dimensions_;
    BlockLduMatrix<Type>::operator-=(fvmv);

    if (faceFluxCorrectionPtr_ && fvmv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ -= *fvmv.faceFluxCorrectionPtr_;
    }
    else if (fvmv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ =
            new GeometricField<Type, fvsPatchField, surfaceMesh>
        (-*fvmv.faceFluxCorrectionPtr_);
    }
}


template<class Type>
void blockFvMatrix<Type>::operator-=(const tmp<blockFvMatrix<Type> >& tfvmv)
{
    operator-=(tfvmv());
    tfvmv.clear();
}


template<class Type>
void blockFvMatrix<Type>::operator+=
(
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(*this, su, "+=");
    BlockLduMatrix<Type>::source() -= su.mesh().V()*su.internalField();
}


template<class Type>
void blockFvMatrix<Type>::operator+=
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    operator+=(tsu());
    tsu.clear();
}


template<class Type>
void blockFvMatrix<Type>::operator-=
(
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(*this, su, "-=");
    BlockLduMatrix<Type>::source() += su.mesh().V()*su.internalField();
}


template<class Type>
void blockFvMatrix<Type>::operator-=
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    operator-=(tsu());
    tsu.clear();
}


template<class Type>
void blockFvMatrix<Type>::operator+=
(
    const dimensioned<Type>& su
)
{
    BlockLduMatrix<Type>::source() -= psi().mesh().V()*su;
}


template<class Type>
void blockFvMatrix<Type>::operator-=
(
    const dimensioned<Type>& su
)
{
    BlockLduMatrix<Type>::source() += psi().mesh().V()*su;
}


template<class Type>
void blockFvMatrix<Type>::operator+=
(
    const zeroField&
)
{}


template<class Type>
void blockFvMatrix<Type>::operator-=
(
    const zeroField&
)
{}


template<class Type>
void blockFvMatrix<Type>::operator*=
(
    const DimensionedField<scalar, volMesh>& dsf
)
{
    dimensions_ *= dsf.dimensions();
    BlockLduMatrix<Type>::operator*=(dsf.field());
    
    if (faceFluxCorrectionPtr_)
    {
        FatalErrorIn("blockFvMatrix<Type>::operator*=(const DimensionedField<scalar, volMesh>&)")
            << "cannot scale a matrix containing a faceFluxCorrection"
            << abort(FatalError);
    }
}


template<class Type>
void blockFvMatrix<Type>::operator*=
(
    const tmp<DimensionedField<scalar, volMesh> >& tdsf
)
{
    operator*=(tdsf());
    tdsf.clear();
}


template<class Type>
void blockFvMatrix<Type>::operator*=
(
    const tmp<volScalarField>& tvsf
)
{
    operator*=(tvsf());
    tvsf.clear();
}


template<class Type>
void blockFvMatrix<Type>::operator*=
(
    const dimensioned<scalar>& ds
)
{
    dimensions_ *= ds.dimensions();
    BlockLduMatrix<Type>::operator*=(ds.value());

    if (faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ *= ds.value();
    }
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

template<class Type>
void checkMethod
(
    const blockFvMatrix<Type>& fvm1,
    const blockFvMatrix<Type>& fvm2,
    const char* op
)
{
    if (&fvm1.psi() != &fvm2.psi())
    {
        FatalErrorIn
        (
            "checkMethod(const blockFvMatrix<Type>&, const blockFvMatrix<Type>&)"
        )   << "incompatible fields for operation "
            << endl << "    "
            << "[" << fvm1.psi().name() << "] "
            << op
            << " [" << fvm2.psi().name() << "]"
            << abort(FatalError);
    }

    if (dimensionSet::debug && fvm1.dimensions() != fvm2.dimensions())
    {
        FatalErrorIn
        (
            "checkMethod(const blockFvMatrix<Type>&, const blockFvMatrix<Type>&)"
        )   << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fvm1.psi().name() << fvm1.dimensions()/dimVolume << " ] "
            << op
            << " [" << fvm2.psi().name() << fvm2.dimensions()/dimVolume << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void checkMethod
(
    const blockFvMatrix<Type>& fvm,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const char* op
)
{
    if (dimensionSet::debug && fvm.dimensions()/dimVolume != vf.dimensions())
    {
        FatalErrorIn
        (
            "checkMethod(const blockFvMatrix<Type>&, const GeometricField<Type, "
            "fvPatchField, volMesh>&)"
        )   <<  "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fvm.psi().name() << fvm.dimensions()/dimVolume << " ] "
            << op
            << " [" << vf.name() << vf.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void checkMethod
(
    const blockFvMatrix<Type>& fvm,
    const dimensioned<Type>& dt,
    const char* op
)
{
    if (dimensionSet::debug && fvm.dimensions()/dimVolume != dt.dimensions())
    {
        FatalErrorIn
        (
            "checkMethod(const blockFvMatrix<Type>&, const dimensioned<Type>&)"
        )   << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fvm.psi().name() << fvm.dimensions()/dimVolume << " ] "
            << op
            << " [" << dt.name() << dt.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
BlockSolverPerformance<Type> solve
(
    blockFvMatrix<Type>& fvm,
    const dictionary& solverControls
)
{
    return fvm.solve(solverControls);
}

template<class Type>
BlockSolverPerformance<Type> solve
(
    const tmp<blockFvMatrix<Type> >& tfvm,
    const dictionary& solverControls
)
{
    BlockSolverPerformance<Type> solverPerf = 
        const_cast<blockFvMatrix<Type>&>(tfvm()).solve(solverControls);

    tfvm.clear();

    return solverPerf;
}


template<class Type>
BlockSolverPerformance<Type> solve(blockFvMatrix<Type>& fvm)
{
    return fvm.solve();
}

template<class Type>
BlockSolverPerformance<Type> solve(const tmp<blockFvMatrix<Type> >& tfvm)
{
    BlockSolverPerformance<Type> solverPerf =
        const_cast<blockFvMatrix<Type>&>(tfvm()).solve();

    tfvm.clear();

    return solverPerf;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const blockFvMatrix<Type>& A,
    const blockFvMatrix<Type>& B
)
{
    checkMethod(A, B, "+");
    tmp<blockFvMatrix<Type> > tC(new blockFvMatrix<Type>(A));
    tC() += B;
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const tmp<blockFvMatrix<Type> >& tA,
    const blockFvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "+");
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC() += B;
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const blockFvMatrix<Type>& A,
    const tmp<blockFvMatrix<Type> >& tB
)
{
    checkMethod(A, tB(), "+");
    tmp<blockFvMatrix<Type> > tC(tB.ptr());
    tC() += A;
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const tmp<blockFvMatrix<Type> >& tA,
    const tmp<blockFvMatrix<Type> >& tB
)
{
    checkMethod(tA(), tB(), "+");
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC() += tB();
    tB.clear();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const blockFvMatrix<Type>& A
)
{
    tmp<blockFvMatrix<Type> > tC(new blockFvMatrix<Type>(A));
    tC().negate();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const tmp<blockFvMatrix<Type> >& tA
)
{
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC().negate();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const blockFvMatrix<Type>& A,
    const blockFvMatrix<Type>& B
)
{
    checkMethod(A, B, "-");
    tmp<blockFvMatrix<Type> > tC(new blockFvMatrix<Type>(A));
    tC() -= B;
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const tmp<blockFvMatrix<Type> >& tA,
    const blockFvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "-");
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC() -= B;
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const blockFvMatrix<Type>& A,
    const tmp<blockFvMatrix<Type> >& tB
)
{
    checkMethod(A, tB(), "-");
    tmp<blockFvMatrix<Type> > tC(tB.ptr());
    tC() -= A;
    tC().negate();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const tmp<blockFvMatrix<Type> >& tA,
    const tmp<blockFvMatrix<Type> >& tB
)
{
    checkMethod(tA(), tB(), "-");
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC() -= tB();
    tB.clear();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const blockFvMatrix<Type>& A,
    const blockFvMatrix<Type>& B
)
{
    checkMethod(A, B, "==");
    return (A - B);
}


template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const tmp<blockFvMatrix<Type> >& tA,
    const blockFvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "==");
    return (tA - B);
}


template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const blockFvMatrix<Type>& A,
    const tmp<blockFvMatrix<Type> >& tB
)
{
    checkMethod(A, tB(), "==");
    return (A - tB);
}


template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const tmp<blockFvMatrix<Type> >& tA,
    const tmp<blockFvMatrix<Type> >& tB
)
{
    checkMethod(tA(), tB(), "==");
    return (tA - tB);
}


template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const blockFvMatrix<Type>& A,
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(A, su, "+");
    tmp<blockFvMatrix<Type> > tC(new blockFvMatrix<Type>(A));
    tC().source() -= su.mesh().V()*su.internalField();
    return tC;
}

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const tmp<blockFvMatrix<Type> >& tA,
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC().source() -= su.mesh().V()*su.internalField();
    return tC;
}

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const blockFvMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(A, tsu(), "+");
    tmp<blockFvMatrix<Type> > tC(new blockFvMatrix<Type>(A));
    tC().source() -= tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const tmp<blockFvMatrix<Type> >& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC().source() -= tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const GeometricField<Type, fvPatchField, volMesh>& su,
    const blockFvMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<blockFvMatrix<Type> > tC(new blockFvMatrix<Type>(A));
    tC().source() -= su.mesh().V()*su.internalField();
    return tC;
}

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const GeometricField<Type, fvPatchField, volMesh>& su,
    const tmp<blockFvMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC().source() -= su.mesh().V()*su.internalField();
    return tC;
}

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu,
    const blockFvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "+");
    tmp<blockFvMatrix<Type> > tC(new blockFvMatrix<Type>(A));
    tC().source() -= tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu,
    const tmp<blockFvMatrix<Type> >& tA
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC().source() -= tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const blockFvMatrix<Type>& A,
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(A, su, "-");
    tmp<blockFvMatrix<Type> > tC(new blockFvMatrix<Type>(A));
    tC().source() += su.mesh().V()*su.internalField();
    return tC;
}

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const tmp<blockFvMatrix<Type> >& tA,
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC().source() += su.mesh().V()*su.internalField();
    return tC;
}

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const blockFvMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(A, tsu(), "-");
    tmp<blockFvMatrix<Type> > tC(new blockFvMatrix<Type>(A));
    tC().source() += tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const tmp<blockFvMatrix<Type> >& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC().source() += tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const GeometricField<Type, fvPatchField, volMesh>& su,
    const blockFvMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<blockFvMatrix<Type> > tC(new blockFvMatrix<Type>(A));
    tC().negate();
    tC().source() -= su.mesh().V()*su.internalField();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const GeometricField<Type, fvPatchField, volMesh>& su,
    const tmp<blockFvMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC().source() -= su.mesh().V()*su.internalField();
    return tC;
}

template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu,
    const blockFvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "-");
    tmp<blockFvMatrix<Type> > tC(new blockFvMatrix<Type>(A));
    tC().negate();
    tC().source() -= tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu,
    const tmp<blockFvMatrix<Type> >& tA
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC().source() -= tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const blockFvMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "+");
    tmp<blockFvMatrix<Type> > tC(new blockFvMatrix<Type>(A));
    tC().source() -= su.value()*A.psi().mesh().V();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const tmp<blockFvMatrix<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC().source() -= su.value()*tC().psi().mesh().V();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const dimensioned<Type>& su,
    const blockFvMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<blockFvMatrix<Type> > tC(new blockFvMatrix<Type>(A));
    tC().source() -= su.value()*A.psi().mesh().V();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator+
(
    const dimensioned<Type>& su,
    const tmp<blockFvMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC().source() -= su.value()*tC().psi().mesh().V();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const blockFvMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "-");
    tmp<blockFvMatrix<Type> > tC(new blockFvMatrix<Type>(A));
    tC().source() += su.value()*tC().psi().mesh().V();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const tmp<blockFvMatrix<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC().source() += su.value()*tC().psi().mesh().V();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const dimensioned<Type>& su,
    const blockFvMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<blockFvMatrix<Type> > tC(new blockFvMatrix<Type>(A));
    tC().negate();
    tC().source() -= su.value()*A.psi().mesh().V();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator-
(
    const dimensioned<Type>& su,
    const tmp<blockFvMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC().source() -= su.value()*tC().psi().mesh().V();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const blockFvMatrix<Type>& A,
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(A, su, "==");
    tmp<blockFvMatrix<Type> > tC(new blockFvMatrix<Type>(A));
    tC().source() += su.mesh().V()*su.internalField();
    return tC;
}

template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const tmp<blockFvMatrix<Type> >& tA,
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC().source() += su.mesh().V()*su.internalField();
    return tC;
}

template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const blockFvMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(A, tsu(), "==");
    tmp<blockFvMatrix<Type> > tC(new blockFvMatrix<Type>(A));
    tC().source() += tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const tmp<blockFvMatrix<Type> >& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "==");
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC().source() += tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const blockFvMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "==");
    tmp<blockFvMatrix<Type> > tC(new blockFvMatrix<Type>(A));
    tC().source() += A.psi().mesh().V()*su.value();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const tmp<blockFvMatrix<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC().source() += tC().psi().mesh().V()*su.value();
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const blockFvMatrix<Type>& A,
    const zeroField&
)
{
    return A;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator==
(
    const tmp<blockFvMatrix<Type> >& tA,
    const zeroField&
)
{
    return tA;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator*
(
    const volScalarField& vsf,
    const blockFvMatrix<Type>& A
)
{
    tmp<blockFvMatrix<Type> > tC(new blockFvMatrix<Type>(A));
    tC() *= vsf;
    return tC;
}

template<class Type>
tmp<blockFvMatrix<Type> > operator*
(
    const tmp<volScalarField>& tvsf,
    const blockFvMatrix<Type>& A
)
{
    tmp<blockFvMatrix<Type> > tC(new blockFvMatrix<Type>(A));
    tC() *= tvsf;
    return tC;
}

template<class Type>
tmp<blockFvMatrix<Type> > operator*
(
    const volScalarField& vsf,
    const tmp<blockFvMatrix<Type> >& tA
)
{
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC() *= vsf;
    return tC;
}

template<class Type>
tmp<blockFvMatrix<Type> > operator*
(
    const tmp<volScalarField>& tvsf,
    const tmp<blockFvMatrix<Type> >& tA
)
{
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC() *= tvsf;
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator*
(
    const dimensioned<scalar>& ds,
    const blockFvMatrix<Type>& A
)
{
    tmp<blockFvMatrix<Type> > tC(new blockFvMatrix<Type>(A));
    tC() *= ds;
    return tC;
}


template<class Type>
tmp<blockFvMatrix<Type> > operator*
(
    const dimensioned<scalar>& ds,
    const tmp<blockFvMatrix<Type> >& tA
)
{
    tmp<blockFvMatrix<Type> > tC(tA.ptr());
    tC() *= ds;
    return tC;
}


template<class Type>
tmp<GeometricField<Type,fvPatchField,volMesh> > operator&
(
    const blockFvMatrix<Type>& M,
    const DimensionedField<Type, volMesh>& psi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tMphi
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "M&"+psi.name(),
                psi.instance(),
                psi.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi.mesh(),
            M.dimensions()/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& Mphi = tMphi();

    CoeffField<Type> saveDiag(M.diag());
    M.addBoundaryDiag(M.diag());
    
    Mphi.internalField() += M.BlockLduMatrix<Type>::H(psi.internalField()) + M.source();
    M.addBoundarySource(Mphi.internalField());

    M.diag() = saveDiag;
    
    Mphi.internalField() /= -psi.mesh().V();
    Mphi.correctBoundaryConditions();

    return tMphi;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > operator&
(
    const blockFvMatrix<Type>& M,
    const tmp<DimensionedField<Type, volMesh> >& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tMpsi = M & tpsi();
    tpsi.clear();
    return tMpsi;
}

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > operator&
(
    const tmp<blockFvMatrix<Type> >& tM,
    const DimensionedField<Type, volMesh>& psi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tMpsi = tM() & psi;
    tM.clear();
    return tMpsi;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > operator&
(
    const tmp<blockFvMatrix<Type> >& tM,
    const tmp<DimensionedField<Type, volMesh> >& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tMpsi = tM() & tpsi();
    tM.clear();
    tpsi.clear();
    return tMpsi;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > operator&
(
    const blockFvMatrix<Type>& M,
    const tmp<GeometricField<Type,fvPatchField,volMesh> >& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tMpsi = M & tpsi();
    tpsi.clear();
    return tMpsi;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > operator&
(
    const tmp<blockFvMatrix<Type> >& tM,
    const tmp<GeometricField<Type,fvPatchField,volMesh> >& tpsi
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tMpsi = tM() & tpsi();
    tM.clear();
    tpsi.clear();
    return tMpsi;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Ostream& operator<<(Ostream& os, const blockFvMatrix<Type>& fvm)
{
    os  << static_cast<const BlockLduMatrix<Type>&>(fvm) << nl
        << "dimensions " << fvm.dimensions_ << endl;

    os.check("Ostream& operator<<(Ostream&, blockFvMatrix<Type>&");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * Solvers * * * * * * * * * * * * * * * * * //

#include "blockFvMatrixSolve.C"

// ************************************************************************* //
