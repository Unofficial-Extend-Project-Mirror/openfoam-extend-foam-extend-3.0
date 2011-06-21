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
#include "blockFvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace blockFvm
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<blockFvMatrix<scalar> > Su
(
    const GeometricField<scalar, fvPatchField, volMesh>& su,
    GeometricField<scalar, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<blockFvMatrix<scalar> > tfvm
    (
        new blockFvMatrix<scalar>
        (
            vf,
            dimVol*su.dimensions()
        )
    );
    blockFvMatrix<scalar>& fvm = tfvm();

    fvm.source() -= mesh.V()*su.internalField()*vf;

    return tfvm;
}

tmp<blockFvMatrix<scalar> > Su
(
    const tmp<GeometricField<scalar, fvPatchField, volMesh> >& tsu,
    GeometricField<scalar, fvPatchField, volMesh>& vf
)
{
    tmp<blockFvMatrix<scalar> > tfvm = Su(tsu(), vf);
    tsu.clear();
    return tfvm;
}


zeroField Su
(
    const zeroField&,
    GeometricField<scalar, fvPatchField, volMesh>&
)
{
    return zeroField();
}


tmp<blockFvMatrix<scalar> > Sp
(
    const GeometricField<scalar, fvPatchField, volMesh>& sp,
    GeometricField<scalar, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<blockFvMatrix<scalar> > tfvm
    (
        new blockFvMatrix<scalar>
        (
            vf,
            dimVol*sp.dimensions()*vf.dimensions()
        )
    );
    blockFvMatrix<scalar>& fvm = tfvm();

    fvm.diag() += mesh.V()*sp.internalField();

    return tfvm;
}


tmp<blockFvMatrix<scalar> > Sp
(
    const tmp<GeometricField<scalar, fvPatchField, volMesh> >& tsp,
    GeometricField<scalar, fvPatchField, volMesh>& vf
)
{
    tmp<blockFvMatrix<scalar> > tfvm = Sp(tsp(), vf);
    tsp.clear();
    return tfvm;
}


tmp<blockFvMatrix<scalar> > Sp
(
    const dimensioned<scalar> sp,
    GeometricField<scalar, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<blockFvMatrix<scalar> > tfvm
    (
        new blockFvMatrix<scalar>
        (
            vf,
            dimVol*sp.dimensions()*vf.dimensions()
        )
    );
    blockFvMatrix<scalar>& fvm = tfvm();

    fvm.diag() += mesh.V()*sp.value();

    return tfvm;
}


zeroField Sp
(
    const zeroField&,
    GeometricField<scalar, fvPatchField, volMesh>&
)
{
    return zeroField();
}


tmp<blockFvMatrix<scalar> > SuSp
(
    const GeometricField<scalar, fvPatchField, volMesh>& sp,
    GeometricField<scalar, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<blockFvMatrix<scalar> > tfvm
    (
        new blockFvMatrix<scalar>
        (
            vf,
            dimVol*sp.dimensions()*vf.dimensions()
        )
    );
    blockFvMatrix<scalar>& fvm = tfvm();

    fvm.diag() += mesh.V()*max(sp.internalField(), pTraits<scalar>::zero);

    fvm.source() -= mesh.V()*(min(sp.internalField(), pTraits<scalar>::zero)
        * vf.internalField());

    return tfvm;
}


tmp<blockFvMatrix<scalar> > SuSp
(
    const tmp<GeometricField<scalar, fvPatchField, volMesh> >& tsp,
    GeometricField<scalar, fvPatchField, volMesh>& vf
)
{
    tmp<blockFvMatrix<scalar> > tfvm = SuSp(tsp(), vf);
    tsp.clear();
    return tfvm;
}


zeroField SuSp
(
    const zeroField&,
    GeometricField<scalar, fvPatchField, volMesh>&
)
{
    return zeroField();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockFvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
