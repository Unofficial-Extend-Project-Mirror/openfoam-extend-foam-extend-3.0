#include "customUserDefines.H"

#ifdef VECTORSPACE_FULL
#   include "VectorSpaceMfull.H"
#elif defined VECTORSPACE_PARTIAL
#   include "VectorSpaceMpartial.H"
#else /*VECTORSPACE_FORLOOPS is defined */
#   include "VectorSpaceMforLoops.H"
#endif
