#ifndef NTL_UTILS_H
#define NTL_UTILS_H
// NTL related helper functions, e.g. print ZZ_p

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#include <vector>
#include <complex>

void showVec(NTL::ZZ const* vals, long size);
void showVec(NTL::ZZ_p const* vals, long size);
void showVec(std::vector<std::complex<double>> const* vals, long size);

NTL::ZZ getZZBal(NTL::ZZ const& zz, NTL::ZZ const& modulus);
void showVecBal(NTL::ZZ_p const* vals, long size);

void showPoly(NTL::ZZ_pE const* poly, long size);
void showPoly(NTL::ZZ_pX const* poly, long size);

NTL::ZZ maxElm(NTL::ZZ const* aX, int n, NTL::ZZ const& modQ);

#endif  // NTL_UTILS_H
