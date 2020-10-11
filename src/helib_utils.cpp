#include "helib_utils.h"

#include <iostream>

double largestCxNorm(std::vector<std::complex<double>> const& vec) {
   double m = 0;
   for (auto& x : vec) {
      if (m < std::abs(x))
         m = std::abs(x);
   }
   return m;
}


void showDCRT(helib::DoubleCRT * dcrt, long size) {
   helib::Context const& context = dcrt->getContext();
   helib::IndexMap<NTL::vec_long> const& map = dcrt->getMap();
   helib::IndexSet const& s = map.getIndexSet();
   NTL::ZZ Q;
   context.productOfPrimes(Q, s);
   std::cout << "[[" << Q << "]]" << std::endl;

   for (long i : s) {
      NTL::vec_long const& row = map[i];
      NTL::zz_pX tmp;
      context.ithModulus(i).iFFT(tmp, row); // inverse FFT
      long phim = row.length();
      long pi = context.ithPrime(i); // the i'th modulus
      std::cout << "[" << pi << "] ";
      for (size_t j = 0; j < phim && j < size; j++) {
         std::cout << tmp.rep[j] << ", ";
      }
      std::cout << std::endl;
   }
}

void showVec(helib::zzX * vals, long size) {
   std::cout << "[";
   std::cout << (*vals)[0];
   for (long i = 1; i < size; ++i) {
      std::cout << ", " << (*vals)[i];
   }
   std::cout << "]" << std::endl;
}

void showCtxtScale(helib::Ctxt const& c, char const* str) {
   std::cout << str << " rf = 2^" << NTL::log(NTL::fabs(c.getRatFactor()))/log(2.0)
             << ", pm = " << c.getPtxtMag()
             << ", log |q_l| = " << c.logOfPrimeSet()/log(2.0) << std::endl;
}

void add(std::vector<cx_double> & out, std::vector<cx_double> const& in0, std::vector<cx_double> const& in1) {
   for (size_t i = 0; i < in0.size(); i++) {
      out[i].real(in0[i].real() + in1[i].real());
      out[i].imag(in0[i].imag() + in1[i].imag());
   }
}

void mul(std::vector<cx_double> & out, std::vector<cx_double> const& in0, std::vector<cx_double> const& in1) {
   for (size_t i = 0; i < in0.size(); i++) {
      out[i].real(in0[i].real() * in1[i].real());
      out[i].imag(in0[i].imag() * in1[i].imag());
   }
}

std::vector<cx_double> diff(std::vector<cx_double> const& in0, std::vector<cx_double> const& in1) {
   std::vector<cx_double> out(in0.size());
   for (size_t i = 0; i < std::min(in0.size(), in1.size()); i++) {
      out[i].real(in0[i].real() - in1[i].real());
      out[i].imag(in0[i].imag() - in1[i].imag());
   }
   return out;
}

NTL::ZZX decrypted;

void copyTo(NTL::ZZX * tgt, NTL::ZZX const* src) {
   *tgt = *src;
}


