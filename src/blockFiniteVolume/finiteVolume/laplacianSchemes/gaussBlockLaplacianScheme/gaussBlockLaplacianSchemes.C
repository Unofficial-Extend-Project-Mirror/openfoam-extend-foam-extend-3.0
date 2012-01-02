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

#include "gaussBlockLaplacianScheme.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    makeBlockFvLaplacianScheme(gaussBlockLaplacianScheme)
}
}

#define declareBlockFvmLaplacianScalarGamma(Type)                             \
                                                                              \
template<>                                                                    \
Foam::tmp<Foam::blockFvMatrix<Foam::Type> >                                   \
Foam::fv::gaussBlockLaplacianScheme<Foam::Type, Foam::scalar>::fvmLaplacian   \
(                                                                             \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& gamma,          \
    GeometricField<Type, fvPatchField, volMesh>& vf                           \
)                                                                             \
{                                                                             \
    const fvMesh& mesh = this->mesh();                                        \
                                                                              \
    GeometricField<scalar, fvsPatchField, surfaceMesh> gammaMagSf             \
    (                                                                         \
        gamma*mesh.magSf()                                                    \
    );                                                                        \
                                                                              \
    tmp<blockFvMatrix<Type> > tfvm                                            \
            = fvmLaplacianUncorrected(gammaMagSf, vf);                        \
    blockFvMatrix<Type>& fvm = tfvm();                                        \
                                                                              \
    if (this->tsnGradScheme_().corrected())                                   \
    {                                                                         \
        if (mesh.fluxRequired(vf.name()))                                     \
        {                                                                     \
            fvm.faceFluxCorrectionPtr() = new                                 \
            GeometricField<Type, fvsPatchField, surfaceMesh>                  \
            (                                                                 \
                gammaMagSf*this->tsnGradScheme_().correction(vf)              \
            );                                                                \
                                                                              \
            fvm.source() -=                                                   \
                mesh.V()*                                                     \
                fvc::div                                                      \
                (                                                             \
                    *fvm.faceFluxCorrectionPtr()                              \
                )().internalField();                                          \
        }                                                                     \
        else                                                                  \
        {                                                                     \
            fvm.source() -=                                                   \
                mesh.V()*                                                     \
                fvc::div                                                      \
                (                                                             \
                    gammaMagSf*this->tsnGradScheme_().correction(vf)          \
                )().internalField();                                          \
        }                                                                     \
    }                                                                         \
                                                                              \
    return tfvm;                                                              \
}                                                                             \
                                                                              \
template<>                                                                    \
Foam::tmp<Foam::GeometricField<Foam::Type, Foam::fvPatchField,                \
    Foam::volMesh> >                                                          \
Foam::fv::gaussBlockLaplacianScheme<Foam::Type, Foam::scalar>::fvcLaplacian   \
(                                                                             \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& gamma,          \
    const GeometricField<Type, fvPatchField, volMesh>& vf                     \
)                                                                             \
{                                                                             \
    const fvMesh& mesh = this->mesh();                                        \
                                                                              \
    surfaceVectorField Sn = mesh.Sf()/mesh.magSf();                           \
                                                                              \
    surfaceVectorField SfGamma = gamma * mesh.Sf();                           \
    GeometricField<scalar, fvsPatchField, surfaceMesh> SfGammaSn              \
        = SfGamma & Sn;                                                       \
    surfaceVectorField SfGammaCorr = SfGamma - SfGammaSn*Sn;                  \
                                                                              \
    tmp<GeometricField<Type, fvPatchField, volMesh> > tLaplacian              \
    (                                                                         \
        fvc::div                                                              \
        (                                                                     \
            SfGammaSn*this->tsnGradScheme_().snGrad(vf)                       \
          + gammaSnGradCorr(SfGammaCorr, vf)                                  \
        )                                                                     \
    );                                                                        \
                                                                              \
    tLaplacian().rename("laplacian(" + gamma.name() + ',' + vf.name() + ')'); \
                                                                              \
    return tLaplacian;                                                        \
}

declareBlockFvmLaplacianScalarGamma(scalar);
declareBlockFvmLaplacianScalarGamma(vector);
declareBlockFvmLaplacianScalarGamma(sphericalTensor);
declareBlockFvmLaplacianScalarGamma(symmTensor);
declareBlockFvmLaplacianScalarGamma(tensor);

#undef declareBlockFvmLaplacianScalarGamma

// ************************************************************************* //
