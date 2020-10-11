#ifndef HELIB_UTILS_H
#define HELIB_UTILS_H
/* Some helper functions for using helib */

#include <vector>
#include <NTL/ZZ_pX.h>
#include <helib/helib.h>

#include "ntl_utils.h"

double largestCxNorm(std::vector<std::complex<double>> const& vec);

void showDCRT(helib::DoubleCRT * dcrt, long size);

void showVec(helib::zzX * vals, long size);

void showCtxtScale(helib::Ctxt const& c, char const* str);

typedef std::complex<double> cx_double;

void mul(std::vector<cx_double> & out, std::vector<cx_double> const& in0, std::vector<cx_double> const& in1);

void add(std::vector<cx_double> & out, std::vector<cx_double> const& in0, std::vector<cx_double> const& in1);

std::vector<cx_double> diff(std::vector<cx_double> const& in0, std::vector<cx_double> const& in1);

extern NTL::ZZX decrypted;
void copyTo(NTL::ZZX * tgt, NTL::ZZX const* src);

#include <iostream>
#include <iterator>

template<typename T>
std::ostream &operator <<(std::ostream &os, const std::vector<T> &v) {
   using namespace std;
   os << "[";
   copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
   os << "]";
   return os;
}


#endif  // HELIB_UTILS_H





