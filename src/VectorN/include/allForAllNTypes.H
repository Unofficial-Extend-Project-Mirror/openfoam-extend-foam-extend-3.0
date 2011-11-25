// This inline include is a short handcombination of:
//  - forAllVectorNTypes.H
//  - forAllTensorNTypes.H
//  - forAllDiagTensorNTypes.H
//  - forAllSphericalTensorNTypes.H
// It does not include: forAllVectorTensorFunctions.H
//
// Must be preceded with the line:
// #define MacroArgs macroName
// where 'macroName' is whatever m (below) is replaced with

#define VectorNMacro(number, m)                                               \
    m(vector##number, Vector##number, args)                                   \
    m(tensor##number, Tensor##number, args)                                   \
    m(diagTensor##number, DiagTensor##number, args)                           \
    m(sphericalTensor##number, SphericalTensor##number, args)

#include "VectorNRunMacro.H"
