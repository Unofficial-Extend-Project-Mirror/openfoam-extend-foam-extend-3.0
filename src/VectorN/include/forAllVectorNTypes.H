// Must be preceded with the line:
// #define MacroArgs macroName
// where 'macroName' is whatever m (below) is replaced with

#define VectorNMacro(number, m)                                               \
    m(vector##number, Vector##number, args)

#include "VectorNRunMacro.H"
