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
    ggiLduInterface

Description
    An abstract base class for ggi coupled interfaces.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    ggiLduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef ggiLduInterface_H
#define ggiLduInterface_H

#include "lduInterface.H"
#include "primitiveFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ggiLduInterface Declaration
\*---------------------------------------------------------------------------*/

class ggiLduInterface
{
protected:

    // Protected data

        //- Transfer buffer
        mutable labelField labelTransferBuffer_;

        //- Field transfer buffer
        mutable scalarField fieldTransferBuffer_;


public:

    //- Runtime type information
    TypeName("ggiLduInterface");

    // Constructors

        //- Construct null
        ggiLduInterface()
        :
            labelTransferBuffer_(),
            fieldTransferBuffer_()
        {}


    // Destructor

        virtual ~ggiLduInterface();


    // Member Functions

        //- Is this the master side?
        virtual bool master() const = 0;

        //- Is this the fine level?
        virtual bool fineLevel() const = 0;

        //- Return shadow interface index
        virtual label shadowIndex() const = 0;

        //- Return shadow interface
        virtual const ggiLduInterface& shadowInterface() const = 0;

        //- Return zone size
        virtual label zoneSize() const = 0;

        //- Return zone addressing
        virtual const labelList& zoneAddressing() const = 0;

        //- Return fine level addressing
        virtual const labelListList& addressing() const = 0;

        //- Is the patch localised on a single processor
        virtual bool localParallel() const = 0;

        //- Return weights
        virtual const scalarListList& weights() const = 0;

        //- Return face transformation tensor
        virtual const tensorField& forwardT() const = 0;

        //- Return face reverse transformation tensor
        virtual const tensorField& reverseT() const = 0;


        // Transfer buffer access

            //- Return contents of the label transfer buffer
            const labelField& labelTransferBuffer() const
            {
                return labelTransferBuffer_;
            }

            //- Return contents of the field transfer buffer
            const scalarField& fieldTransferBuffer() const
            {
                return fieldTransferBuffer_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
