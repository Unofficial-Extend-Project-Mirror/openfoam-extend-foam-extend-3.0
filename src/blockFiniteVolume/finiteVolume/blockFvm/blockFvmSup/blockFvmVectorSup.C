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
#include "expandTensor.H"
#include "expandDiagTensor.H"
#include "expandSphericalTensor.H"
#include "expandTensorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace blockFvm
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<blockFvMatrix<vector> > Su
(
    const GeometricField<tensor, fvPatchField, volMesh>& su,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<blockFvMatrix<vector> > tfvm
    (
        new blockFvMatrix<vector>
        (
            vf,
            dimVol*su.dimensions()
        )
    );
    blockFvMatrix<vector>& fvm = tfvm();

    fvm.source() -= mesh.V()*(su.internalField() & vf);

    return tfvm;
}


tmp<blockFvMatrix<vector> > Su
(
    const GeometricField<diagTensor, fvPatchField, volMesh>& su,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<blockFvMatrix<vector> > tfvm
    (
        new blockFvMatrix<vector>
        (
            vf,
            dimVol*su.dimensions()
        )
    );
    blockFvMatrix<vector>& fvm = tfvm();

    fvm.source() -= mesh.V()*(su.internalField() & vf);

    return tfvm;
}


tmp<blockFvMatrix<vector> > Su
(
    const GeometricField<sphericalTensor, fvPatchField, volMesh>& su,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<blockFvMatrix<vector> > tfvm
    (
        new blockFvMatrix<vector>
        (
            vf,
            dimVol*su.dimensions()
        )
    );
    blockFvMatrix<vector>& fvm = tfvm();

    fvm.source() -= mesh.V()*(su.internalField() & vf);

    return tfvm;
}


tmp<blockFvMatrix<vector> > Su
(
    const tmp<GeometricField<tensor, fvPatchField, volMesh> >& tsu,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    tmp<blockFvMatrix<vector> > tfvm = Su(tsu(), vf);
    tsu.clear();
    return tfvm;
}


tmp<blockFvMatrix<vector> > Su
(
    const tmp<GeometricField<diagTensor, fvPatchField, volMesh> >& tsu,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    tmp<blockFvMatrix<vector> > tfvm = Su(tsu(), vf);
    tsu.clear();
    return tfvm;
}


tmp<blockFvMatrix<vector> > Su
(
    const tmp<GeometricField<sphericalTensor, fvPatchField, volMesh> >& tsu,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    tmp<blockFvMatrix<vector> > tfvm = Su(tsu(), vf);
    tsu.clear();
    return tfvm;
}


zeroField Su
(
    const zeroField&,
    GeometricField<vector, fvPatchField, volMesh>&
)
{
    return zeroField();
}


tmp<blockFvMatrix<vector> > Sp
(
    const GeometricField<tensor, fvPatchField, volMesh>& sp,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<blockFvMatrix<vector> > tfvm
    (
        new blockFvMatrix<vector>
        (
            vf,
            dimVol*sp.dimensions()*vf.dimensions()
        )
    );
    blockFvMatrix<vector>& fvm = tfvm();

    fvm.diag() += mesh.V()*sp.internalField();

    return tfvm;
}


tmp<blockFvMatrix<vector> > Sp
(
    const GeometricField<diagTensor, fvPatchField, volMesh>& sp,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<blockFvMatrix<vector> > tfvm
    (
        new blockFvMatrix<vector>
        (
            vf,
            dimVol*sp.dimensions()*vf.dimensions()
        )
    );
    blockFvMatrix<vector>& fvm = tfvm();
    
    tmp<vectorField> tLinearSp(new vectorField(mesh.nCells()));
    
    contractLinear(tLinearSp(), sp.internalField());
    
    fvm.diag() += mesh.V()*tLinearSp;

    return tfvm;
}


tmp<blockFvMatrix<vector> > Sp
(
    const GeometricField<sphericalTensor, fvPatchField, volMesh>& sp,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<blockFvMatrix<vector> > tfvm
    (
        new blockFvMatrix<vector>
        (
            vf,
            dimVol*sp.dimensions()*vf.dimensions()
        )
    );
    blockFvMatrix<vector>& fvm = tfvm();

    tmp<scalarField> tScalarSp(new scalarField(mesh.nCells()));
    
    contractScalar(tScalarSp(), sp.internalField());

    fvm.diag() += mesh.V()*tScalarSp;

    return tfvm;
}


tmp<blockFvMatrix<vector> > Sp
(
    const tmp<GeometricField<tensor, fvPatchField, volMesh> >& tsp,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    tmp<blockFvMatrix<vector> > tfvm = Sp(tsp(), vf);
    tsp.clear();
    return tfvm;
}


tmp<blockFvMatrix<vector> > Sp
(
    const tmp<GeometricField<diagTensor, fvPatchField, volMesh> >& tsp,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    tmp<blockFvMatrix<vector> > tfvm = Sp(tsp(), vf);
    tsp.clear();
    return tfvm;
}


tmp<blockFvMatrix<vector> > Sp
(
    const tmp<GeometricField<sphericalTensor, fvPatchField, volMesh> >& tsp,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    tmp<blockFvMatrix<vector> > tfvm = Sp(tsp(), vf);
    tsp.clear();
    return tfvm;
}


tmp<blockFvMatrix<vector> > Sp
(
    const dimensioned<tensor> sp,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<blockFvMatrix<vector> > tfvm
    (
        new blockFvMatrix<vector>
        (
            vf,
            dimVol*sp.dimensions()*vf.dimensions()
        )
    );
    blockFvMatrix<vector>& fvm = tfvm();

    fvm.diag() += mesh.V()*sp.value();

    return tfvm;
}


tmp<blockFvMatrix<vector> > Sp
(
    const dimensioned<diagTensor> sp,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<blockFvMatrix<vector> > tfvm
    (
        new blockFvMatrix<vector>
        (
            vf,
            dimVol*sp.dimensions()*vf.dimensions()
        )
    );
    blockFvMatrix<vector>& fvm = tfvm();
        
    fvm.diag() += mesh.V()*contractLinear(sp.value());

    return tfvm;
}


tmp<blockFvMatrix<vector> > Sp
(
    const dimensioned<sphericalTensor> sp,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<blockFvMatrix<vector> > tfvm
    (
        new blockFvMatrix<vector>
        (
            vf,
            dimVol*sp.dimensions()*vf.dimensions()
        )
    );
    blockFvMatrix<vector>& fvm = tfvm();

    fvm.diag() += mesh.V()*contractScalar(sp.value());

    return tfvm;
}


zeroField Sp
(
    const zeroField&,
    GeometricField<vector, fvPatchField, volMesh>&
)
{
    return zeroField();
}


tmp<blockFvMatrix<vector> > SuSp
(
    const GeometricField<tensor, fvPatchField, volMesh>& sp,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<blockFvMatrix<vector> > tfvm
    (
        new blockFvMatrix<vector>
        (
            vf,
            dimVol*sp.dimensions()*vf.dimensions()
        )
    );
    blockFvMatrix<vector>& fvm = tfvm();

    fvm.diag() += mesh.V()*max(sp.internalField(), tensor::zero);

    fvm.source() -= mesh.V()*(min(sp.internalField(), tensor::zero)
        & vf.internalField());

    return tfvm;
}


tmp<blockFvMatrix<vector> > SuSp
(
    const GeometricField<diagTensor, fvPatchField, volMesh>& sp,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<blockFvMatrix<vector> > tfvm
    (
        new blockFvMatrix<vector>
        (
            vf,
            dimVol*sp.dimensions()*vf.dimensions()
        )
    );
    blockFvMatrix<vector>& fvm = tfvm();

    tmp<vectorField> tLinearSp(new vectorField(mesh.nCells()));    
    contractLinear(tLinearSp(), max(sp.internalField(), diagTensor::zero)());
    
    fvm.diag() += mesh.V()*tLinearSp;

    fvm.source() -= mesh.V()*(min(sp.internalField(), diagTensor::zero)
        & vf.internalField());

    return tfvm;
}


tmp<blockFvMatrix<vector> > SuSp
(
    const GeometricField<sphericalTensor, fvPatchField, volMesh>& sp,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<blockFvMatrix<vector> > tfvm
    (
        new blockFvMatrix<vector>
        (
            vf,
            dimVol*sp.dimensions()*vf.dimensions()
        )
    );
    blockFvMatrix<vector>& fvm = tfvm();

    tmp<scalarField> tScalarSp(new scalarField(mesh.nCells()));    
    contractScalar(tScalarSp(), max(sp.internalField(), sphericalTensor::zero)());

    fvm.diag() += mesh.V()*tScalarSp;

    fvm.source() -= mesh.V()*(min(sp.internalField(), sphericalTensor::zero)
        & vf.internalField());

    return tfvm;
}


tmp<blockFvMatrix<vector> > SuSp
(
    const tmp<GeometricField<tensor, fvPatchField, volMesh> >& tsp,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    tmp<blockFvMatrix<vector> > tfvm = SuSp(tsp(), vf);
    tsp.clear();
    return tfvm;
}


tmp<blockFvMatrix<vector> > SuSp
(
    const tmp<GeometricField<diagTensor, fvPatchField, volMesh> >& tsp,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    tmp<blockFvMatrix<vector> > tfvm = SuSp(tsp(), vf);
    tsp.clear();
    return tfvm;
}


tmp<blockFvMatrix<vector> > SuSp
(
    const tmp<GeometricField<sphericalTensor, fvPatchField, volMesh> >& tsp,
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    tmp<blockFvMatrix<vector> > tfvm = SuSp(tsp(), vf);
    tsp.clear();
    return tfvm;
}


zeroField SuSp
(
    const zeroField&,
    GeometricField<vector, fvPatchField, volMesh>&
)
{
    return zeroField();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockFvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
