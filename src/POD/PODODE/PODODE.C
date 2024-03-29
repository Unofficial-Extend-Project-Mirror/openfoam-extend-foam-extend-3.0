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

Class
    PODODE

\*---------------------------------------------------------------------------*/

#include "PODODE.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PODODE, 0);
    defineRunTimeSelectionTable(PODODE, dictionary);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::PODODE> Foam::PODODE::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    word PODTypeName = dict.lookup("type");

    Info<< "Selecting POD ODE model " << PODTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(PODTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "PODODE::New\n"
            "(\n"
            "    const fvMesh& mesh\n"
            "    const dictionary& dict\n"
            ")",
            dict
        )   << "Unknown PODODE type "
            << PODTypeName << endl << endl
            << "Valid  POD ODEs are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<PODODE>
        (cstrIter()(mesh, dict.subDict(PODTypeName + "Coeffs")));
}


// ************************************************************************* //
