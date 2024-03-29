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

\*---------------------------------------------------------------------------*/

#include "DevolatilisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::autoPtr<Foam::DevolatilisationModel<CloudType> >
Foam::DevolatilisationModel<CloudType>::New
(
    const dictionary& dict,
    CloudType& owner
)
{
    word DevolatilisationModelType(dict.lookup("DevolatilisationModel"));

    Info<< "Selecting DevolatilisationModel " << DevolatilisationModelType
        << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(DevolatilisationModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "DevolatilisationModel<CloudType>::New"
            "("
                "const dictionary&, "
                "CloudType&"
            ")"
        )   << "Unknown DevolatilisationModelType type "
            << DevolatilisationModelType
            << ", constructor not in hash table" << nl << nl
            << "    Valid DevolatilisationModel types are:" << nl
            << dictionaryConstructorTablePtr_->toc() << exit(FatalError);
    }

    return autoPtr<DevolatilisationModel<CloudType> >(cstrIter()(dict, owner));
}


// ************************************************************************* //
