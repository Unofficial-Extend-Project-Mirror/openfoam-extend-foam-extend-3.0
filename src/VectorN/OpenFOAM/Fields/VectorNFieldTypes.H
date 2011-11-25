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
    Include the header files for all the primitive types that Fields are
    instantiated for.

\*---------------------------------------------------------------------------*/

#ifndef VectorNFieldTypes_H
#define VectorNFieldTypes_H

#include "customUserDefines.H"

#define vectorTensorTypeDefs(size)                                            \
    typedef VectorN<scalar, size> vector##size;                               \
    typedef TensorN<scalar, size> tensor##size;                               \
    typedef SphericalTensorN<scalar, size> sphericalTensor##size;             \
    typedef DiagTensorN<scalar, size> diagTensor##size;                       \
    template<>                                                                \
    inline bool contiguous<vector##size>() {return true;}                     \
    template<>                                                                \
    inline bool contiguous<tensor##size>() {return true;}                     \
    template<>                                                                \
    inline bool contiguous<sphericalTensor##size>() {return true;}            \
    template<>                                                                \
    inline bool contiguous<diagTensor##size>() {return true;}

#include "scalar.H"
#include "contiguous.H"
#include "TensorN.H"

// Specialized implementations are loaded here
#ifdef VECTORN_2
#   include "vector2.H"
#   include "tensor2.H"
#   include "diagTensor2.H"
#endif
#ifdef VECTORN_4
#   include "vector4.H"
#   include "tensor4.H"
#   include "diagTensor4.H"
#endif

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Generic typedefs - created using a macroForLoop, with a conditional start

#define VectorNMacro vectorTensorTypeDefs
#define SKIP_VECTORN_2
#define SKIP_VECTORN_4
#include "VectorNRunMacro.H"
#undef SKIP_VECTORN_2
#undef SKIP_VECTORN_4

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#undef vectorTensorTypeDefs

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
