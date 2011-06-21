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

#include "volFields.H"
#include "surfaceFields.H"
#include "blockFvMatrix.H"
#include "blockLaplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace blockFvm
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<blockFvMatrix<Type> >
laplacian
(
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    surfaceScalarField Gamma
    (
        IOobject
        (
            "1",
            vf.time().constant(),
            vf.mesh(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        dimensionedScalar("1", dimless, 1.0)
    );

    return blockFvm::laplacian(Gamma, vf, name);
}


template<class Type>
tmp<blockFvMatrix<Type> >
laplacian
(
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    surfaceScalarField Gamma
    (
        IOobject
        (
            "1",
            vf.time().constant(),
            vf.mesh(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        dimensionedScalar("1", dimless, 1.0)
    );

    return blockFvm::laplacian
    (
        Gamma,
        vf,
        "laplacian(" + vf.name() + ')'
    );
}


template<class Type>
tmp<blockFvMatrix<Type> >
laplacian
(
    const zeroField&,
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return tmp<blockFvMatrix<Type> >
    (
        new blockFvMatrix<Type>(vf, dimensionSet(0, 0, -2, 0, 0))
    );
}


template<class Type>
tmp<blockFvMatrix<Type> >
laplacian
(
    const zeroField&,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return tmp<blockFvMatrix<Type> >
    (
        new blockFvMatrix<Type>(vf, dimensionSet(0, 0, -2, 0, 0))
    );
}


template<class Type>
tmp<blockFvMatrix<Type> >
laplacian
(
    const geometricOneField&,
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return blockFvm::laplacian(vf, name);
}


template<class Type>
tmp<blockFvMatrix<Type> >
laplacian
(
    const geometricOneField&,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return blockFvm::laplacian(vf);
}


template<class Type, class GType>
tmp<blockFvMatrix<Type> >
laplacian
(
    const dimensioned<GType>& gamma,
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    GeometricField<GType, fvsPatchField, surfaceMesh> Gamma
    (
        IOobject
        (
            gamma.name(),
            vf.instance(),
            vf.mesh(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        gamma
    );

    return blockFvm::laplacian(Gamma, vf, name);
}


template<class Type, class GType>
tmp<blockFvMatrix<Type> >
laplacian
(
    const dimensioned<GType>& gamma,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    GeometricField<GType, fvsPatchField, surfaceMesh> Gamma
    (
        IOobject
        (
            gamma.name(),
            vf.instance(),
            vf.mesh(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        gamma
    );

    return blockFvm::laplacian(Gamma, vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<blockFvMatrix<Type> >
laplacian
(
    const GeometricField<GType, fvPatchField, volMesh>& gamma,
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::blockLaplacianScheme<Type, GType>::New
    (
        vf.mesh(),
        vf.mesh().laplacianScheme(name)
    )().fvmLaplacian(gamma, vf);
}


template<class Type, class GType>
tmp<blockFvMatrix<Type> >
laplacian
(
    const tmp<GeometricField<GType, fvPatchField, volMesh> >& tgamma,
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    tmp<blockFvMatrix<Type> > Laplacian(blockFvm::laplacian(tgamma(), vf, name));
    tgamma.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<blockFvMatrix<Type> >
laplacian
(
    const GeometricField<GType, fvPatchField, volMesh>& gamma,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return blockFvm::laplacian
    (
        gamma,
        vf,
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp<blockFvMatrix<Type> >
laplacian
(
    const tmp<GeometricField<GType, fvPatchField, volMesh> >& tgamma,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<blockFvMatrix<Type> > Laplacian(blockFvm::laplacian(tgamma(), vf));
    tgamma.clear();
    return Laplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<blockFvMatrix<Type> >
laplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::blockLaplacianScheme<Type, GType>::New
    (
        vf.mesh(),
        vf.mesh().laplacianScheme(name)
    )().fvmLaplacian(gamma, vf);
}


template<class Type, class GType>
tmp<blockFvMatrix<Type> >
laplacian
(
    const tmp<GeometricField<GType, fvsPatchField, surfaceMesh> >& tgamma,
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    tmp<blockFvMatrix<Type> > tLaplacian = blockFvm::laplacian(tgamma(), vf, name);
    tgamma.clear();
    return tLaplacian;
}


template<class Type, class GType>
tmp<blockFvMatrix<Type> >
laplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return blockFvm::laplacian
    (
        gamma,
        vf,
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp<blockFvMatrix<Type> >
laplacian
(
    const tmp<GeometricField<GType, fvsPatchField, surfaceMesh> >& tGamma,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<blockFvMatrix<Type> > tfvm(blockFvm::laplacian(tGamma(), vf));
    tGamma.clear();
    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
