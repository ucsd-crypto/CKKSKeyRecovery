#ifndef PALISADE_UTILS
#define PALISADE_UTILS
// Helper functions for palisade

#include "utils/inttypes.h"
#include "lattice/elemparams.h"
#include "lattice/ilparams.h"
#include "lattice/ildcrtparams.h"
#include "lattice/ilelement.h"

using namespace lbcrypto;

void printDCRTPoly(DCRTPoly const& p, size_t num);

void printNativePoly(NativePoly const& p, size_t num);

#endif
