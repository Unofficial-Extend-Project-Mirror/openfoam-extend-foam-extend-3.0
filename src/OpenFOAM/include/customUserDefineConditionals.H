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

Description
    Compiler defines that depend on those set in customUserDefines.H are set
    here.  customUserDefines.H is the user-editable file; this one is not.

\*---------------------------------------------------------------------------*/

// VectorSpaceN define error checking
#ifdef VECTORSPACE_FULL
#   define VECTORSPACE_SET
#endif
#ifdef VECTORSPACE_PARTIAL
#   ifndef VECTORSPACE_UNROLLN
#       error VECTORSPACE_PARTIAL requires VECTORSPACE_UNROLLN to be set in \
            the customUserDefines.H file
#   elif VECTORSPACE_UNROLLN < 2
#       error VECTORSPACE_UNROLLN must be set to a value of 2 or greater
#   endif
#   ifdef VECTORSPACE_SET
#       error Only define one of: VECTORSPACE_FULL, VECTORSPACE_PARTIAL, \
            VECTORSPACE_FORLOOPS in the customUserDefines.H file
#   endif
#   define VECTORSPACE_SET
#endif
#ifdef VECTORSPACE_FORLOOPS
#   ifdef VECTORSPACE_SET
#       error Only define one of: VECTORSPACE_FULL, VECTORSPACE_PARTIAL, \
            VECTORSPACE_FORLOOPS in the customUserDefines.H file
#   endif
#   define VECTORSPACE_SET
#endif
#ifndef VECTORSPACE_SET
#   error One of : VECTORSPACE_FULL, VECTORSPACE_PARTIAL, VECTORSPACE_FORLOOPS\
        must be set in the customUserDefines.H file
#endif
#undef VECTORSPACE_SET

