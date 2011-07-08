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

#include "gaussBlockConvectionScheme.H"
#include "fvcSurfaceIntegrate.H"
#include "blockFvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class Type>
tmp<blockFvMatrix<Type> >
gaussBlockConvectionScheme<Type>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<surfaceScalarField> tweights = tinterpScheme_().weights(vf);
    const surfaceScalarField& weights = tweights();

    tmp<blockFvMatrix<Type> > tfvm
    (
        new blockFvMatrix<Type>
        (
            vf,
            faceFlux.dimensions()*vf.dimensions()
        )
    );
    blockFvMatrix<Type>& fvm = tfvm();

    fvm.lower() = -weights.internalField()*faceFlux.internalField();
    fvm.upper() = fvm.lower().asScalar() + faceFlux.internalField();
    fvm.negSumDiag();

    forAll(fvm.psi().boundaryField(), patchI)
    {
        const fvPatchField<Type>& psf = fvm.psi().boundaryField()[patchI];
        const fvsPatchScalarField& patchFlux = faceFlux.boundaryField()[patchI];
        const fvsPatchScalarField& pw = weights.boundaryField()[patchI];

        fvm.interfaceDiag()[patchI] = patchFlux*psf.valueInternalCoeffs(pw);
        fvm.interfaceSource()[patchI] = -patchFlux*psf.valueBoundaryCoeffs(pw);
        
        if(psf.coupled())
        {
            fvm.coupleUpper()[patchI] = -patchFlux*psf.valueUpperCoeffs(pw);
            fvm.coupleLower()[patchI] = -patchFlux*psf.valueLowerCoeffs(pw);
        }
    }

    if (tinterpScheme_().corrected())
    {
        fvm += fvc::surfaceIntegrate(faceFlux*tinterpScheme_().correction(vf));
    }

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
