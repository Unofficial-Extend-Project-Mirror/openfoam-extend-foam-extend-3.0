// Must be preceded with the line:
// #define MacroArgs macroName
// where 'macroName' is whatever m (below) is replaced with

#define VectorNMacro(number, m)                                               \
    m                                                                         \
    (                                                                         \
        tensor##number, diagTensor##number, sphericalTensor##number,          \
        vector##number, scalar, args                                          \
    )

#include "VectorNRunMacro.H"
