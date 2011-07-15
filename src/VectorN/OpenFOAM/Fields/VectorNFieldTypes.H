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

#define VECTORN_1_10

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

#include "contiguous.H"
#include "TensorN.H"
#include "scalar.H"

// Specialized implementations are loaded here
#ifdef VECTORN_1_10
#   include "vector2.H"
#   include "vector4.H"
#   include "tensor2.H"
#   include "tensor4.H"
#   include "diagTensor2.H"
#   include "diagTensor4.H"
#endif

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Generic typedefs
#ifdef VECTORN_1_10
vectorTensorTypeDefs(5)
vectorTensorTypeDefs(6)
vectorTensorTypeDefs(7)
vectorTensorTypeDefs(8)
vectorTensorTypeDefs(9)
vectorTensorTypeDefs(10)
#endif

#ifdef VECTORN_6_15
vectorTensorTypeDefs(6)
vectorTensorTypeDefs(7)
vectorTensorTypeDefs(8)
vectorTensorTypeDefs(9)
vectorTensorTypeDefs(10)
vectorTensorTypeDefs(11)
vectorTensorTypeDefs(12)
vectorTensorTypeDefs(13)
vectorTensorTypeDefs(14)
vectorTensorTypeDefs(15)
#endif

#ifdef VECTORN_11_20
vectorTensorTypeDefs(11)
vectorTensorTypeDefs(12)
vectorTensorTypeDefs(13)
vectorTensorTypeDefs(14)
vectorTensorTypeDefs(15)
vectorTensorTypeDefs(16)
vectorTensorTypeDefs(17)
vectorTensorTypeDefs(18)
vectorTensorTypeDefs(19)
vectorTensorTypeDefs(20)
#endif

#ifdef VECTORN_16_25
vectorTensorTypeDefs(16)
vectorTensorTypeDefs(17)
vectorTensorTypeDefs(18)
vectorTensorTypeDefs(19)
vectorTensorTypeDefs(20)
vectorTensorTypeDefs(21)
vectorTensorTypeDefs(22)
vectorTensorTypeDefs(23)
vectorTensorTypeDefs(24)
vectorTensorTypeDefs(25)
#endif

#ifdef VECTORN_21_30
vectorTensorTypeDefs(21)
vectorTensorTypeDefs(22)
vectorTensorTypeDefs(23)
vectorTensorTypeDefs(24)
vectorTensorTypeDefs(25)
vectorTensorTypeDefs(26)
vectorTensorTypeDefs(27)
vectorTensorTypeDefs(28)
vectorTensorTypeDefs(29)
vectorTensorTypeDefs(30)
#endif

#ifdef VECTORN_24_35
vectorTensorTypeDefs(26)
vectorTensorTypeDefs(27)
vectorTensorTypeDefs(28)
vectorTensorTypeDefs(29)
vectorTensorTypeDefs(30)
vectorTensorTypeDefs(31)
vectorTensorTypeDefs(32)
vectorTensorTypeDefs(33)
vectorTensorTypeDefs(34)
vectorTensorTypeDefs(35)
#endif

#ifdef VECTORN_31_40
vectorTensorTypeDefs(31)
vectorTensorTypeDefs(32)
vectorTensorTypeDefs(33)
vectorTensorTypeDefs(34)
vectorTensorTypeDefs(35)
vectorTensorTypeDefs(36)
vectorTensorTypeDefs(37)
vectorTensorTypeDefs(38)
vectorTensorTypeDefs(39)
vectorTensorTypeDefs(40)
#endif

#ifdef VECTORN_36_45
vectorTensorTypeDefs(36)
vectorTensorTypeDefs(37)
vectorTensorTypeDefs(38)
vectorTensorTypeDefs(39)
vectorTensorTypeDefs(40)
vectorTensorTypeDefs(41)
vectorTensorTypeDefs(42)
vectorTensorTypeDefs(43)
vectorTensorTypeDefs(44)
vectorTensorTypeDefs(45)
#endif

#ifdef VECTORN_41_50
vectorTensorTypeDefs(41)
vectorTensorTypeDefs(42)
vectorTensorTypeDefs(43)
vectorTensorTypeDefs(44)
vectorTensorTypeDefs(45)
vectorTensorTypeDefs(46)
vectorTensorTypeDefs(47)
vectorTensorTypeDefs(48)
vectorTensorTypeDefs(49)
vectorTensorTypeDefs(50)
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef VECTORN_1_10
#define forAllVectorNTypes(m, args...)                                        \
    m(vector2, Vector2, args)                                                 \
    m(vector4, Vector4, args)                                                 \
    m(vector5, Vector5, args)                                                 \
    m(vector6, Vector6, args)                                                 \
    m(vector7, Vector7, args)                                                 \
    m(vector8, Vector8, args)                                                 \
    m(vector9, Vector9, args)                                                 \
    m(vector10, Vector10, args)
#endif

#ifdef VECTORN_6_15
#define forAllVectorNTypes(m, args...)                                        \
    m(vector6, Vector6, args)                                                 \
    m(vector7, Vector7, args)                                                 \
    m(vector8, Vector8, args)                                                 \
    m(vector9, Vector9, args)                                                 \
    m(vector10, Vector10, args)                                               \
    m(vector11, Vector11, args)                                               \
    m(vector12, Vector12, args)                                               \
    m(vector13, Vector13, args)                                               \
    m(vector14, Vector14, args)                                               \
    m(vector15, Vector15, args)
#endif

#ifdef VECTORN_11_20
#define forAllVectorNTypes(m, args...)                                        \
    m(vector11, Vector11, args)                                               \
    m(vector12, Vector12, args)                                               \
    m(vector13, Vector13, args)                                               \
    m(vector14, Vector14, args)                                               \
    m(vector15, Vector15, args)                                               \
    m(vector16, Vector16, args)                                               \
    m(vector17, Vector17, args)                                               \
    m(vector18, Vector18, args)                                               \
    m(vector19, Vector19, args)                                               \
    m(vector20, Vector20, args)
#endif

#ifdef VECTORN_16_25
#define forAllVectorNTypes(m, args...)                                        \
    m(vector16, Vector16, args)                                               \
    m(vector17, Vector17, args)                                               \
    m(vector18, Vector18, args)                                               \
    m(vector19, Vector19, args)                                               \
    m(vector20, Vector20, args)                                               \
    m(vector21, Vector21, args)                                               \
    m(vector22, Vector22, args)                                               \
    m(vector23, Vector23, args)                                               \
    m(vector24, Vector24, args)                                               \
    m(vector25, Vector25, args)
#endif

#ifdef VECTORN_21_30
#define forAllVectorNTypes(m, args...)                                        \
    m(vector21, Vector21, args)                                               \
    m(vector22, Vector22, args)                                               \
    m(vector23, Vector23, args)                                               \
    m(vector24, Vector24, args)                                               \
    m(vector25, Vector25, args)                                               \
    m(vector26, Vector26, args)                                               \
    m(vector27, Vector27, args)                                               \
    m(vector28, Vector28, args)                                               \
    m(vector29, Vector29, args)                                               \
    m(vector30, Vector30, args)
#endif

#ifdef VECTORN_26_35
#define forAllVectorNTypes(m, args...)                                        \
    m(vector26, Vector26, args)                                               \
    m(vector27, Vector27, args)                                               \
    m(vector28, Vector28, args)                                               \
    m(vector29, Vector29, args)                                               \
    m(vector30, Vector30, args)                                               \
    m(vector31, Vector31, args)                                               \
    m(vector32, Vector32, args)                                               \
    m(vector33, Vector33, args)                                               \
    m(vector34, Vector34, args)                                               \
    m(vector35, Vector35, args)
#endif

#ifdef VECTORN_31_40
#define forAllVectorNTypes(m, args...)                                        \
    m(vector31, Vector31, args)                                               \
    m(vector32, Vector32, args)                                               \
    m(vector33, Vector33, args)                                               \
    m(vector34, Vector34, args)                                               \
    m(vector35, Vector35, args)                                               \
    m(vector36, Vector36, args)                                               \
    m(vector37, Vector37, args)                                               \
    m(vector38, Vector38, args)                                               \
    m(vector39, Vector39, args)                                               \
    m(vector40, Vector40, args)
#endif

#ifdef VECTORN_36_45
#define forAllVectorNTypes(m, args...)                                        \
    m(vector36, Vector36, args)                                               \
    m(vector37, Vector37, args)                                               \
    m(vector38, Vector38, args)                                               \
    m(vector39, Vector39, args)                                               \
    m(vector40, Vector40, args)                                               \
    m(vector41, Vector41, args)                                               \
    m(vector42, Vector42, args)                                               \
    m(vector43, Vector43, args)                                               \
    m(vector44, Vector44, args)                                               \
    m(vector45, Vector45, args)
#endif

#ifdef VECTORN_41_50
#define forAllVectorNTypes(m, args...)                                        \
    m(vector41, Vector41, args)                                               \
    m(vector42, Vector42, args)                                               \
    m(vector43, Vector43, args)                                               \
    m(vector44, Vector44, args)                                               \
    m(vector45, Vector45, args)                                               \
    m(vector46, Vector46, args)                                               \
    m(vector47, Vector47, args)                                               \
    m(vector48, Vector48, args)                                               \
    m(vector49, Vector49, args)                                               \
    m(vector50, Vector50, args)
#endif

#ifdef VECTORN_1_10
#define forAllTensorNTypes(m, args...)                                        \
    m(tensor2, Tensor2, args)                                                 \
    m(tensor4, Tensor4, args)                                                 \
    m(tensor5, Tensor5, args)                                                 \
    m(tensor6, Tensor6, args)                                                 \
    m(tensor7, Tensor7, args)                                                 \
    m(tensor8, Tensor8, args)                                                 \
    m(tensor9, Tensor9, args)                                                 \
    m(tensor10, Tensor10, args)
#endif

#ifdef VECTORN_6_15
#define forAllTensorNTypes(m, args...)                                        \
    m(tensor6, Tensor6, args)                                                 \
    m(tensor7, Tensor7, args)                                                 \
    m(tensor8, Tensor8, args)                                                 \
    m(tensor9, Tensor9, args)                                                 \
    m(tensor10, Tensor10, args)                                               \
    m(tensor11, Tensor11, args)                                               \
    m(tensor12, Tensor12, args)                                               \
    m(tensor13, Tensor13, args)                                               \
    m(tensor14, Tensor14, args)                                               \
    m(tensor15, Tensor15, args)
#endif

#ifdef VECTORN_11_20
#define forAllTensorNTypes(m, args...)                                        \
    m(tensor11, Tensor11, args)                                               \
    m(tensor12, Tensor12, args)                                               \
    m(tensor13, Tensor13, args)                                               \
    m(tensor14, Tensor14, args)                                               \
    m(tensor15, Tensor15, args)                                               \
    m(tensor16, Tensor16, args)                                               \
    m(tensor17, Tensor17, args)                                               \
    m(tensor18, Tensor18, args)                                               \
    m(tensor19, Tensor19, args)                                               \
    m(tensor20, Tensor20, args)
#endif

#ifdef VECTORN_16_25
#define forAllTensorNTypes(m, args...)                                        \
    m(tensor16, Tensor16, args)                                               \
    m(tensor17, Tensor17, args)                                               \
    m(tensor18, Tensor18, args)                                               \
    m(tensor19, Tensor19, args)                                               \
    m(tensor20, Tensor20, args)                                               \
    m(tensor21, Tensor21, args)                                               \
    m(tensor22, Tensor22, args)                                               \
    m(tensor23, Tensor23, args)                                               \
    m(tensor24, Tensor24, args)                                               \
    m(tensor25, Tensor25, args)
#endif

#ifdef VECTORN_21_30
#define forAllTensorNTypes(m, args...)                                        \
    m(tensor21, Tensor21, args)                                               \
    m(tensor22, Tensor22, args)                                               \
    m(tensor23, Tensor23, args)                                               \
    m(tensor24, Tensor24, args)                                               \
    m(tensor25, Tensor25, args)                                               \
    m(tensor26, Tensor26, args)                                               \
    m(tensor27, Tensor27, args)                                               \
    m(tensor28, Tensor28, args)                                               \
    m(tensor29, Tensor29, args)                                               \
    m(tensor30, Tensor30, args)
#endif

#ifdef VECTORN_26_35
#define forAllTensorNTypes(m, args...)                                        \
    m(tensor26, Tensor26, args)                                               \
    m(tensor27, Tensor27, args)                                               \
    m(tensor28, Tensor28, args)                                               \
    m(tensor29, Tensor29, args)                                               \
    m(tensor30, Tensor30, args)                                               \
    m(tensor31, Tensor31, args)                                               \
    m(tensor32, Tensor32, args)                                               \
    m(tensor33, Tensor33, args)                                               \
    m(tensor34, Tensor34, args)                                               \
    m(tensor35, Tensor35, args)
#endif

#ifdef VECTORN_31_40
#define forAllTensorNTypes(m, args...)                                        \
    m(tensor31, Tensor31, args)                                               \
    m(tensor32, Tensor32, args)                                               \
    m(tensor33, Tensor33, args)                                               \
    m(tensor34, Tensor34, args)                                               \
    m(tensor35, Tensor35, args)                                               \
    m(tensor36, Tensor36, args)                                               \
    m(tensor37, Tensor37, args)                                               \
    m(tensor38, Tensor38, args)                                               \
    m(tensor39, Tensor39, args)                                               \
    m(tensor40, Tensor40, args)
#endif

#ifdef VECTORN_36_45
#define forAllTensorNTypes(m, args...)                                        \
    m(tensor36, Tensor36, args)                                               \
    m(tensor37, Tensor37, args)                                               \
    m(tensor38, Tensor38, args)                                               \
    m(tensor39, Tensor39, args)                                               \
    m(tensor40, Tensor40, args)                                               \
    m(tensor41, Tensor41, args)                                               \
    m(tensor42, Tensor42, args)                                               \
    m(tensor43, Tensor43, args)                                               \
    m(tensor44, Tensor44, args)                                               \
    m(tensor45, Tensor45, args)
#endif

#ifdef VECTORN_41_50
#define forAllTensorNTypes(m, args...)                                        \
    m(tensor41, Tensor41, args)                                               \
    m(tensor42, Tensor42, args)                                               \
    m(tensor43, Tensor43, args)                                               \
    m(tensor44, Tensor44, args)                                               \
    m(tensor45, Tensor45, args)                                               \
    m(tensor46, Tensor46, args)                                               \
    m(tensor47, Tensor47, args)                                               \
    m(tensor48, Tensor48, args)                                               \
    m(tensor49, Tensor49, args)                                               \
    m(tensor50, Tensor50, args)
#endif

#ifdef VECTORN_1_10
#define forAllDiagTensorNTypes(m, args...)                                    \
    m(diagTensor2, DiagTensor2, args)                                         \
    m(diagTensor4, DiagTensor4, args)                                         \
    m(diagTensor5, DiagTensor5, args)                                         \
    m(diagTensor6, DiagTensor6, args)                                         \
    m(diagTensor7, DiagTensor7, args)                                         \
    m(diagTensor8, DiagTensor8, args)                                         \
    m(diagTensor9, DiagTensor9, args)                                         \
    m(diagTensor10, DiagTensor10, args)
#endif

#ifdef VECTORN_6_15
#define forAllDiagTensorNTypes(m, args...)                                    \
    m(diagTensor6, DiagTensor6, args)                                         \
    m(diagTensor7, DiagTensor7, args)                                         \
    m(diagTensor8, DiagTensor8, args)                                         \
    m(diagTensor9, DiagTensor9, args)                                         \
    m(diagTensor10, DiagTensor10, args)                                       \
    m(diagTensor11, DiagTensor11, args)                                       \
    m(diagTensor12, DiagTensor12, args)                                       \
    m(diagTensor13, DiagTensor13, args)                                       \
    m(diagTensor14, DiagTensor14, args)                                       \
    m(diagTensor15, DiagTensor15, args)
#endif

#ifdef VECTORN_11_20
#define forAllDiagTensorNTypes(m, args...)                                    \
    m(diagTensor11, DiagTensor11, args)                                       \
    m(diagTensor12, DiagTensor12, args)                                       \
    m(diagTensor13, DiagTensor13, args)                                       \
    m(diagTensor14, DiagTensor14, args)                                       \
    m(diagTensor15, DiagTensor15, args)                                       \
    m(diagTensor16, DiagTensor16, args)                                       \
    m(diagTensor17, DiagTensor17, args)                                       \
    m(diagTensor18, DiagTensor18, args)                                       \
    m(diagTensor19, DiagTensor19, args)                                       \
    m(diagTensor20, DiagTensor20, args)
#endif

#ifdef VECTORN_16_25
#define forAllDiagTensorNTypes(m, args...)                                    \
    m(diagTensor16, DiagTensor16, args)                                       \
    m(diagTensor17, DiagTensor17, args)                                       \
    m(diagTensor18, DiagTensor18, args)                                       \
    m(diagTensor19, DiagTensor19, args)                                       \
    m(diagTensor20, DiagTensor20, args)                                       \
    m(diagTensor21, DiagTensor21, args)                                       \
    m(diagTensor22, DiagTensor22, args)                                       \
    m(diagTensor23, DiagTensor23, args)                                       \
    m(diagTensor24, DiagTensor24, args)                                       \
    m(diagTensor25, DiagTensor25, args)
#endif

#ifdef VECTORN_21_30
#define forAllDiagTensorNTypes(m, args...)                                    \
    m(diagTensor21, DiagTensor21, args)                                       \
    m(diagTensor22, DiagTensor22, args)                                       \
    m(diagTensor23, DiagTensor23, args)                                       \
    m(diagTensor24, DiagTensor24, args)                                       \
    m(diagTensor25, DiagTensor25, args)                                       \
    m(diagTensor26, DiagTensor26, args)                                       \
    m(diagTensor27, DiagTensor27, args)                                       \
    m(diagTensor28, DiagTensor28, args)                                       \
    m(diagTensor29, DiagTensor29, args)                                       \
    m(diagTensor30, DiagTensor30, args)
#endif

#ifdef VECTORN_26_35
#define forAllDiagTensorNTypes(m, args...)                                    \
    m(diagTensor26, DiagTensor26, args)                                       \
    m(diagTensor27, DiagTensor27, args)                                       \
    m(diagTensor28, DiagTensor28, args)                                       \
    m(diagTensor29, DiagTensor29, args)                                       \
    m(diagTensor30, DiagTensor30, args)                                       \
    m(diagTensor31, DiagTensor31, args)                                       \
    m(diagTensor32, DiagTensor32, args)                                       \
    m(diagTensor33, DiagTensor33, args)                                       \
    m(diagTensor34, DiagTensor34, args)                                       \
    m(diagTensor35, DiagTensor35, args)
#endif

#ifdef VECTORN_31_40
#define forAllDiagTensorNTypes(m, args...)                                    \
    m(diagTensor31, DiagTensor31, args)                                       \
    m(diagTensor32, DiagTensor32, args)                                       \
    m(diagTensor33, DiagTensor33, args)                                       \
    m(diagTensor34, DiagTensor34, args)                                       \
    m(diagTensor35, DiagTensor35, args)                                       \
    m(diagTensor36, DiagTensor36, args)                                       \
    m(diagTensor37, DiagTensor37, args)                                       \
    m(diagTensor38, DiagTensor38, args)                                       \
    m(diagTensor39, DiagTensor39, args)                                       \
    m(diagTensor40, DiagTensor40, args)
#endif

#ifdef VECTORN_36_45
#define forAllDiagTensorNTypes(m, args...)                                    \
    m(diagTensor36, DiagTensor36, args)                                       \
    m(diagTensor37, DiagTensor37, args)                                       \
    m(diagTensor38, DiagTensor38, args)                                       \
    m(diagTensor39, DiagTensor39, args)                                       \
    m(diagTensor40, DiagTensor40, args)                                       \
    m(diagTensor41, DiagTensor41, args)                                       \
    m(diagTensor42, DiagTensor42, args)                                       \
    m(diagTensor43, DiagTensor43, args)                                       \
    m(diagTensor44, DiagTensor44, args)                                       \
    m(diagTensor45, DiagTensor45, args)
#endif

#ifdef VECTORN_41_50
#define forAllDiagTensorNTypes(m, args...)                                    \
    m(diagTensor41, DiagTensor41, args)                                       \
    m(diagTensor42, DiagTensor42, args)                                       \
    m(diagTensor43, DiagTensor43, args)                                       \
    m(diagTensor44, DiagTensor44, args)                                       \
    m(diagTensor45, DiagTensor45, args)                                       \
    m(diagTensor46, DiagTensor46, args)                                       \
    m(diagTensor47, DiagTensor47, args)                                       \
    m(diagTensor48, DiagTensor48, args)                                       \
    m(diagTensor49, DiagTensor49, args)                                       \
    m(diagTensor50, DiagTensor50, args)
#endif

#ifdef VECTORN_1_10
#define forAllSphericalTensorNTypes(m, args...)                               \
    m(sphericalTensor2, SphericalTensor2, args)                               \
    m(sphericalTensor4, SphericalTensor4, args)                               \
    m(sphericalTensor5, SphericalTensor5, args)                               \
    m(sphericalTensor6, SphericalTensor6, args)                               \
    m(sphericalTensor7, SphericalTensor7, args)                               \
    m(sphericalTensor8, SphericalTensor8, args)                               \
    m(sphericalTensor9, SphericalTensor9, args)                               \
    m(sphericalTensor10, SphericalTensor10, args)
#endif

#ifdef VECTORN_6_15
#define forAllSphericalTensorNTypes(m, args...)                               \
    m(sphericalTensor6, SphericalTensor6, args)                               \
    m(sphericalTensor7, SphericalTensor7, args)                               \
    m(sphericalTensor8, SphericalTensor8, args)                               \
    m(sphericalTensor9, SphericalTensor9, args)                               \
    m(sphericalTensor10, SphericalTensor10, args)                             \
    m(sphericalTensor11, SphericalTensor11, args)                             \
    m(sphericalTensor12, SphericalTensor12, args)                             \
    m(sphericalTensor13, SphericalTensor13, args)                             \
    m(sphericalTensor14, SphericalTensor14, args)                             \
    m(sphericalTensor15, SphericalTensor15, args)
#endif

#ifdef VECTORN_11_20
#define forAllSphericalTensorNTypes(m, args...)                               \
    m(sphericalTensor11, SphericalTensor11, args)                             \
    m(sphericalTensor12, SphericalTensor12, args)                             \
    m(sphericalTensor13, SphericalTensor13, args)                             \
    m(sphericalTensor14, SphericalTensor14, args)                             \
    m(sphericalTensor15, SphericalTensor15, args)                             \
    m(sphericalTensor16, SphericalTensor16, args)                             \
    m(sphericalTensor17, SphericalTensor17, args)                             \
    m(sphericalTensor18, SphericalTensor18, args)                             \
    m(sphericalTensor19, SphericalTensor19, args)                             \
    m(sphericalTensor20, SphericalTensor20, args)
#endif

#ifdef VECTORN_16_25
#define forAllSphericalTensorNTypes(m, args...)                               \
    m(sphericalTensor16, SphericalTensor16, args)                             \
    m(sphericalTensor17, SphericalTensor17, args)                             \
    m(sphericalTensor18, SphericalTensor18, args)                             \
    m(sphericalTensor19, SphericalTensor19, args)                             \
    m(sphericalTensor20, SphericalTensor20, args)                             \
    m(sphericalTensor21, SphericalTensor21, args)                             \
    m(sphericalTensor22, SphericalTensor22, args)                             \
    m(sphericalTensor23, SphericalTensor23, args)                             \
    m(sphericalTensor24, SphericalTensor24, args)                             \
    m(sphericalTensor25, SphericalTensor25, args)
#endif

#ifdef VECTORN_21_30
#define forAllSphericalTensorNTypes(m, args...)                               \
    m(sphericalTensor21, SphericalTensor21, args)                             \
    m(sphericalTensor22, SphericalTensor22, args)                             \
    m(sphericalTensor23, SphericalTensor23, args)                             \
    m(sphericalTensor24, SphericalTensor24, args)                             \
    m(sphericalTensor25, SphericalTensor25, args)                             \
    m(sphericalTensor26, SphericalTensor26, args)                             \
    m(sphericalTensor27, SphericalTensor27, args)                             \
    m(sphericalTensor28, SphericalTensor28, args)                             \
    m(sphericalTensor29, SphericalTensor29, args)                             \
    m(sphericalTensor30, SphericalTensor30, args)
#endif

#ifdef VECTORN_26_35
#define forAllSphericalTensorNTypes(m, args...)                               \
    m(sphericalTensor26, SphericalTensor26, args)                             \
    m(sphericalTensor27, SphericalTensor27, args)                             \
    m(sphericalTensor28, SphericalTensor28, args)                             \
    m(sphericalTensor29, SphericalTensor29, args)                             \
    m(sphericalTensor30, SphericalTensor30, args)                             \
    m(sphericalTensor31, SphericalTensor31, args)                             \
    m(sphericalTensor32, SphericalTensor32, args)                             \
    m(sphericalTensor33, SphericalTensor33, args)                             \
    m(sphericalTensor34, SphericalTensor34, args)                             \
    m(sphericalTensor35, SphericalTensor35, args)
#endif

#ifdef VECTORN_31_40
#define forAllSphericalTensorNTypes(m, args...)                               \
    m(sphericalTensor31, SphericalTensor31, args)                             \
    m(sphericalTensor32, SphericalTensor32, args)                             \
    m(sphericalTensor33, SphericalTensor33, args)                             \
    m(sphericalTensor34, SphericalTensor34, args)                             \
    m(sphericalTensor35, SphericalTensor35, args)                             \
    m(sphericalTensor36, SphericalTensor36, args)                             \
    m(sphericalTensor37, SphericalTensor37, args)                             \
    m(sphericalTensor38, SphericalTensor38, args)                             \
    m(sphericalTensor39, SphericalTensor39, args)                             \
    m(sphericalTensor40, SphericalTensor40, args)
#endif

#ifdef VECTORN_36_45
#define forAllSphericalTensorNTypes(m, args...)                               \
    m(sphericalTensor36, SphericalTensor36, args)                             \
    m(sphericalTensor37, SphericalTensor37, args)                             \
    m(sphericalTensor38, SphericalTensor38, args)                             \
    m(sphericalTensor39, SphericalTensor39, args)                             \
    m(sphericalTensor40, SphericalTensor40, args)                             \
    m(sphericalTensor41, SphericalTensor41, args)                             \
    m(sphericalTensor42, SphericalTensor42, args)                             \
    m(sphericalTensor43, SphericalTensor43, args)                             \
    m(sphericalTensor44, SphericalTensor44, args)                             \
    m(sphericalTensor45, SphericalTensor45, args)
#endif

#ifdef VECTORN_41_50
#define forAllSphericalTensorNTypes(m, args...)                               \
    m(sphericalTensor41, SphericalTensor41, args)                             \
    m(sphericalTensor42, SphericalTensor42, args)                             \
    m(sphericalTensor43, SphericalTensor43, args)                             \
    m(sphericalTensor44, SphericalTensor44, args)                             \
    m(sphericalTensor45, SphericalTensor45, args)                             \
    m(sphericalTensor46, SphericalTensor46, args)                             \
    m(sphericalTensor47, SphericalTensor47, args)                             \
    m(sphericalTensor48, SphericalTensor48, args)                             \
    m(sphericalTensor49, SphericalTensor49, args)                             \
    m(sphericalTensor50, SphericalTensor50, args)
#endif

#ifdef VECTORN_1_10
#define forAllVectorTensorNTypes(m, args...)                                  \
    m(tensor2, diagTensor2, sphericalTensor2, vector2, scalar, args)          \
    m(tensor4, diagTensor4, sphericalTensor4, vector4, scalar, args)          \
    m(tensor5, diagTensor5, sphericalTensor5, vector5, scalar, args)          \
    m(tensor6, diagTensor6, sphericalTensor6, vector6, scalar, args)          \
    m(tensor7, diagTensor7, sphericalTensor7, vector7, scalar, args)          \
    m(tensor8, diagTensor8, sphericalTensor8, vector8, scalar, args)          \
    m(tensor9, diagTensor9, sphericalTensor9, vector9, scalar, args)          \
    m(tensor10, diagTensor10, sphericalTensor10, vector10, scalar, args)
#endif

#ifdef VECTORN_6_15
#define forAllVectorTensorNTypes(m, args...)                                  \
    m(tensor6, diagTensor6, sphericalTensor6, vector6, scalar, args)          \
    m(tensor7, diagTensor7, sphericalTensor7, vector7, scalar, args)          \
    m(tensor8, diagTensor8, sphericalTensor8, vector8, scalar, args)          \
    m(tensor9, diagTensor9, sphericalTensor9, vector9, scalar, args)          \
    m(tensor10, diagTensor10, sphericalTensor10, vector10, scalar, args)      \
    m(tensor11, diagTensor11, sphericalTensor11, vector11, scalar, args)      \
    m(tensor12, diagTensor12, sphericalTensor12, vector12, scalar, args)      \
    m(tensor13, diagTensor13, sphericalTensor13, vector13, scalar, args)      \
    m(tensor14, diagTensor14, sphericalTensor14, vector14, scalar, args)      \
    m(tensor15, diagTensor15, sphericalTensor15, vector15, scalar, args)
#endif

#ifdef VECTORN_11_20
#define forAllVectorTensorNTypes(m, args...)                                  \
    m(tensor11, diagTensor11, sphericalTensor11, vector11, scalar, args)      \
    m(tensor12, diagTensor12, sphericalTensor12, vector12, scalar, args)      \
    m(tensor13, diagTensor13, sphericalTensor13, vector13, scalar, args)      \
    m(tensor14, diagTensor14, sphericalTensor14, vector14, scalar, args)      \
    m(tensor15, diagTensor15, sphericalTensor15, vector15, scalar, args)      \
    m(tensor16, diagTensor16, sphericalTensor16, vector16, scalar, args)      \
    m(tensor17, diagTensor17, sphericalTensor17, vector17, scalar, args)      \
    m(tensor18, diagTensor18, sphericalTensor18, vector18, scalar, args)      \
    m(tensor19, diagTensor19, sphericalTensor19, vector19, scalar, args)      \
    m(tensor20, diagTensor20, sphericalTensor20, vector20, scalar, args)
#endif

#ifdef VECTORN_16_25
#define forAllVectorTensorNTypes(m, args...)                                  \
    m(tensor16, diagTensor16, sphericalTensor16, vector16, scalar, args)      \
    m(tensor17, diagTensor17, sphericalTensor17, vector17, scalar, args)      \
    m(tensor18, diagTensor18, sphericalTensor18, vector18, scalar, args)      \
    m(tensor19, diagTensor19, sphericalTensor19, vector19, scalar, args)      \
    m(tensor20, diagTensor20, sphericalTensor20, vector20, scalar, args)      \
    m(tensor21, diagTensor21, sphericalTensor21, vector21, scalar, args)      \
    m(tensor22, diagTensor22, sphericalTensor22, vector22, scalar, args)      \
    m(tensor23, diagTensor23, sphericalTensor23, vector23, scalar, args)      \
    m(tensor24, diagTensor24, sphericalTensor24, vector24, scalar, args)      \
    m(tensor25, diagTensor25, sphericalTensor25, vector25, scalar, args)
#endif

#ifdef VECTORN_21_30
#define forAllVectorTensorNTypes(m, args...)                                  \
    m(tensor21, diagTensor21, sphericalTensor21, vector21, scalar, args)      \
    m(tensor22, diagTensor22, sphericalTensor22, vector22, scalar, args)      \
    m(tensor23, diagTensor23, sphericalTensor23, vector23, scalar, args)      \
    m(tensor24, diagTensor24, sphericalTensor24, vector24, scalar, args)      \
    m(tensor25, diagTensor25, sphericalTensor25, vector25, scalar, args)      \
    m(tensor26, diagTensor26, sphericalTensor26, vector26, scalar, args)      \
    m(tensor27, diagTensor27, sphericalTensor27, vector27, scalar, args)      \
    m(tensor28, diagTensor28, sphericalTensor28, vector28, scalar, args)      \
    m(tensor29, diagTensor29, sphericalTensor29, vector29, scalar, args)      \
    m(tensor30, diagTensor30, sphericalTensor30, vector30, scalar, args)
#endif

#ifdef VECTORN_26_35
#define forAllVectorTensorNTypes(m, args...)                                  \
    m(tensor26, diagTensor26, sphericalTensor26, vector26, scalar, args)      \
    m(tensor27, diagTensor27, sphericalTensor27, vector27, scalar, args)      \
    m(tensor28, diagTensor28, sphericalTensor28, vector28, scalar, args)      \
    m(tensor29, diagTensor29, sphericalTensor29, vector29, scalar, args)      \
    m(tensor30, diagTensor30, sphericalTensor30, vector30, scalar, args)      \
    m(tensor31, diagTensor31, sphericalTensor31, vector31, scalar, args)      \
    m(tensor32, diagTensor32, sphericalTensor32, vector32, scalar, args)      \
    m(tensor33, diagTensor33, sphericalTensor33, vector33, scalar, args)      \
    m(tensor34, diagTensor34, sphericalTensor34, vector34, scalar, args)      \
    m(tensor35, diagTensor35, sphericalTensor35, vector35, scalar, args)
#endif

#ifdef VECTORN_31_40
#define forAllVectorTensorNTypes(m, args...)                                  \
    m(tensor31, diagTensor31, sphericalTensor31, vector31, scalar, args)      \
    m(tensor32, diagTensor32, sphericalTensor32, vector32, scalar, args)      \
    m(tensor33, diagTensor33, sphericalTensor33, vector33, scalar, args)      \
    m(tensor34, diagTensor34, sphericalTensor34, vector34, scalar, args)      \
    m(tensor35, diagTensor35, sphericalTensor35, vector35, scalar, args)      \
    m(tensor36, diagTensor36, sphericalTensor36, vector36, scalar, args)      \
    m(tensor37, diagTensor37, sphericalTensor37, vector37, scalar, args)      \
    m(tensor38, diagTensor38, sphericalTensor38, vector38, scalar, args)      \
    m(tensor39, diagTensor39, sphericalTensor39, vector39, scalar, args)      \
    m(tensor40, diagTensor40, sphericalTensor40, vector40, scalar, args)
#endif

#ifdef VECTORN_36_45
#define forAllVectorTensorNTypes(m, args...)                                  \
    m(tensor36, diagTensor36, sphericalTensor36, vector36, scalar, args)      \
    m(tensor37, diagTensor37, sphericalTensor37, vector37, scalar, args)      \
    m(tensor38, diagTensor38, sphericalTensor38, vector38, scalar, args)      \
    m(tensor39, diagTensor39, sphericalTensor39, vector39, scalar, args)      \
    m(tensor40, diagTensor40, sphericalTensor40, vector40, scalar, args)      \
    m(tensor41, diagTensor41, sphericalTensor41, vector41, scalar, args)      \
    m(tensor42, diagTensor42, sphericalTensor42, vector42, scalar, args)      \
    m(tensor43, diagTensor43, sphericalTensor43, vector43, scalar, args)      \
    m(tensor44, diagTensor44, sphericalTensor44, vector44, scalar, args)      \
    m(tensor45, diagTensor45, sphericalTensor45, vector45, scalar, args)
#endif

#ifdef VECTORN_41_50
#define forAllVectorTensorNTypes(m, args...)                                  \
    m(tensor41, diagTensor41, sphericalTensor41, vector41, scalar, args)      \
    m(tensor42, diagTensor42, sphericalTensor42, vector42, scalar, args)      \
    m(tensor43, diagTensor43, sphericalTensor43, vector43, scalar, args)      \
    m(tensor44, diagTensor44, sphericalTensor44, vector44, scalar, args)      \
    m(tensor45, diagTensor45, sphericalTensor45, vector45, scalar, args)      \
    m(tensor46, diagTensor46, sphericalTensor46, vector46, scalar, args)      \
    m(tensor47, diagTensor47, sphericalTensor47, vector47, scalar, args)      \
    m(tensor48, diagTensor48, sphericalTensor48, vector48, scalar, args)      \
    m(tensor49, diagTensor49, sphericalTensor49, vector49, scalar, args)      \
    m(tensor50, diagTensor50, sphericalTensor50, vector50, scalar, args)
#endif

#undef vectorTensorTypeDefs

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //