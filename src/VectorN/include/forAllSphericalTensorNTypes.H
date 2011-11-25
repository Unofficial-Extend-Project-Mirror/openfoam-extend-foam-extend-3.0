// Must be preceded with the line:
// #define MacroArgs macroName
// where 'macroName' is whatever m (below) is replaced with

#define VectorNMacro(number, m)                                               \
    m(sphericalTensor##number, SphericalTensor##number, args)

#include "VectorNRunMacro.H"
