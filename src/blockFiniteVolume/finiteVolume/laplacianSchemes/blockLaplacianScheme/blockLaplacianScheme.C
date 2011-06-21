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

#include "fv.H"
#include "HashTable.H"
#include "linear.H"
#include "blockFvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<blockLaplacianScheme<Type, GType> > blockLaplacianScheme<Type, GType>::New
(
    const fvMesh& mesh,
    Istream& schemeData
)
{
    if (fv::debug)
    {
        Info<< "blockLaplacianScheme<Type, GType>::New(const fvMesh&, Istream&) : "
               "constructing blockLaplacianScheme<Type, GType>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorIn
        (
            "blockLaplacianScheme<Type, GType>::New(const fvMesh&, Istream&)",
            schemeData
        )   << "Laplacian scheme not specified" << endl << endl
            << "Valid laplacian schemes are :" << endl
            << IstreamConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    word schemeName(schemeData);

    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(schemeName);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "blockLaplacianScheme<Type, GType>::New(const fvMesh&, Istream&)",
            schemeData
        )   << "unknown laplacian scheme " << schemeName << endl << endl
            << "Valid laplacian schemes are :" << endl
            << IstreamConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return cstrIter()(mesh, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, class GType>
blockLaplacianScheme<Type, GType>::~blockLaplacianScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<blockFvMatrix<Type> >
blockLaplacianScheme<Type, GType>::fvmLaplacian
(
    const GeometricField<GType, fvPatchField, volMesh>& gamma,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvmLaplacian(tinterpGammaScheme_().interpolate(gamma)(), vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
