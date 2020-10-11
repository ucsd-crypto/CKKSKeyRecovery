#ifndef PALISADE_UTILS
#define PALISADE_UTILS
// Helper functions for palisade

#include "palisadecore.h"
#include "lattice/ilparams.h"
#include "lattice/ildcrtparams.h"
#include "lattice/poly.h"
#include "lattice/dcrtpoly.h"

using namespace lbcrypto;

void printDCRTPoly(DCRTPoly const& p, size_t num);

void printNativePoly(NativePoly const& p, size_t num);

void printPoly(Poly const& p, size_t num);

#endif
