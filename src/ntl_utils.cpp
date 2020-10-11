#include "ntl_utils.h"

#include <iostream>

void showVec(NTL::ZZ const* vals, long size) {
   std::cout << "[";
   std::cout << vals[0];
   for (long i = 1; i < size; ++i) {
      std::cout << ", " << vals[i];
   }
   std::cout << "]" << std::endl;
}


void showVec(NTL::ZZ_p const* vals, long size) {
   std::cout << "{" << NTL::ZZ_p::modulus() << "} [";
   std::cout << vals[0];
   for (long i = 1; i < size; ++i) {
      std::cout << ", " << vals[i];
   }
   std::cout << "]" << std::endl;
}

void showVec(std::vector<std::complex<double>>  const* vals, long size) {
   std::cout << "[";
   std::cout << (*vals)[0];
   for (long i = 1; i < size; ++i) {
      std::cout << ", " << (*vals)[i];
   }
   std::cout << "]" << std::endl;
}

NTL::ZZ getZZBal(NTL::ZZ const& zz, NTL::ZZ const& modulus) {
   NTL::ZZ res = zz;
   if (res >= modulus / 2) {
      res = res - modulus;
   }
   return res;
}

void showVecBal(NTL::ZZ_p const* vals, long size) {
   std::cout << "{" << NTL::ZZ_p::modulus() << "} [";
   std::cout << getZZBal(NTL::rep(vals[0]), NTL::ZZ_p::modulus());
   for (long i = 1; i < size; ++i) {
      std::cout << ", " << getZZBal(NTL::rep(vals[i]), NTL::ZZ_p::modulus());
   }
   std::cout << "]" << std::endl;
}

NTL::ZZ maxElm(NTL::ZZ const* aX, int n, NTL::ZZ const& modQ) {
   NTL::ZZ m = NTL::ZZ::zero();
   NTL::ZZ hQ = modQ / 2;
   for (int i=0; i<n; i++) {
      NTL::ZZ d = aX[i] % modQ;
      if (d >= hQ) {
         d = d - modQ;
      }
      d = abs(d);
      if (m<d) {
         m = d;
      }
   }
   return m;
}

