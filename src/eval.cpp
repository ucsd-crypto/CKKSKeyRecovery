#include "eval.h"
#include <cmath>
#include <cstring>

#include <NTL/ZZ.h>
#include <iostream>

void evalPlainAdd(std::vector<cx_double> & res,
                  std::vector<cx_double> const& in0, std::vector<cx_double> const& in1) {
   size_t len = std::min(in0.size(), in1.size());
   res.resize(len);
   for (size_t i = 0; i < len; i++) {
      res[i] = in0[i] + in1[i];
   }
}


void evalPlainMul(std::vector<cx_double> & res,
                  std::vector<cx_double> const& in0, std::vector<cx_double> const& in1) {
   size_t len = std::min(in0.size(), in1.size());
   res.resize(len);
   for (size_t i = 0; i < len; i++) {
      res[i] = in0[i] * in1[i];
   }
}

void evalPlainNegate(std::vector<cx_double> & res, std::vector<cx_double> const& in) {
   size_t len = in.size();
   res.resize(len);
   for (size_t i = 0; i < len; i++) {
      res[i] = -in[i];
   }
}

void evalPlainInverse(std::vector<cx_double> & res, std::vector<cx_double> const& in) {
   size_t len = in.size();
   res.resize(len);
   for (size_t i = 0; i < len; i++) {
      res[i] = 1.0 / in[i];
   }
}

void evalPlainPowerOf2(std::vector<cx_double> & res, std::vector<cx_double> const& in, size_t logDeg) {
   res = in;                    // copy all the numbers
   for (size_t j = 0; j < logDeg; j++) {
      for (size_t i = 0; i < in.size(); i++) {
         res[i] = res[i] * res[i];
      }
   }
}

void evalPlainPower(std::vector<cx_double> & res, std::vector<cx_double> const& in, size_t deg) {
   size_t logDeg = (size_t)floor(std::log2((double)deg));
   size_t remDeg = deg - (1 << logDeg);
   evalPlainPowerOf2(res, in, logDeg);
   if (remDeg > 0) {
      std::vector<cx_double> tmp(in.size());
      evalPlainPower(tmp, in, remDeg);
      evalPlainMul(res, res, tmp);
   }
}

void evalPlainAddi(std::vector<cx_double> & res, std::vector<cx_double> const& in, double c) {
   res.resize(in.size());
   for (size_t i = 0; i < in.size(); i++) {
      res[i] = in[i] + c;
   }
}

void evalPlainMuli(std::vector<cx_double> & res, std::vector<cx_double> const& in, double c) {
   res.resize(in.size());
   for (size_t i = 0; i < in.size(); i++) {
      res[i] = in[i] * c;
   }
}

std::map<SpecialFunction::FuncName, std::vector<double>>
SpecialFunction::coeffsOf = {
   { FuncName::LOG, {0,1,-0.5,1./3,-1./4,1./5,-1./6,1./7,-1./8,1./9,-1./10} },
   { FuncName::EXP, {1,1,0.5,1./6,1./24,1./120,1./720,1./5040,1./40320,1./362880,1./3628800 } },
   { FuncName::SIGMOID, {1./2,1./4,0,-1./48,0,1./480,0,-17./80640,0,31./1451520,0} }
};


void evalPlainFunc(std::vector<cx_double> & res, std::vector<cx_double> const& in, std::vector<double> const& coeff, int evDeg) {
   res.resize(in.size(), 0);
   res = in;                          // x

   evalPlainMuli(res, res, coeff[1]); // c_1 x
   evalPlainAddi(res, res, coeff[0]); // c_1 x + c_0

   const int deg = evDeg == -1 ? coeff.size()-1 : evDeg;
   const int logDeg = (int)floor(std::log2(deg));
   std::vector<std::vector<cx_double>> basis(logDeg+1); // x^(2^i) for i=0..logDeg
   basis[0] = in;                                       // x^(2^0)
   for (int j = 0, i = 1; j < logDeg; j++, i++) {
      evalPlainPowerOf2(basis[i], in, i); // x^(2^i)
   }

   for (int i = 2; i <= deg; i++) {
      int k = floor(std::log2(i));
      int r = i - (1 << k);
      std::vector<cx_double> tmp = basis[k]; // x^[2^k]
      while (r > 0) {
         k = floor(std::log2(r));
         r = r - (1 << k);
         evalPlainMul(tmp, tmp, basis[k]);
      }
      evalPlainMuli(tmp, tmp, coeff[i]); // c_i * x^i
      evalPlainAdd(res, res, tmp);
   }
}

void evalPlainFunc(std::vector<cx_double> & res, std::vector<cx_double> const& in, SpecialFunction::FuncName name, int evDeg) {
   std::vector<double> const& coeff = SpecialFunction::coeffsOf[name];
   evalPlainFunc(res, in, coeff, evDeg);
}

void evalPlainFunc(std::vector<cx_double> & res, cx_double * in, size_t len, SpecialFunction::FuncName name, int evDeg) {
   std::vector<cx_double> vin(in, in+len);
   evalPlainFunc(res, vin, name, evDeg);
}

double largestElm(std::vector<std::complex<double>> const& vec) {
   double m = 0;
   for (auto& x : vec) {
      if (m < std::abs(x.real()))
         m = std::abs(x.real());
      if (m < std::abs(x.imag()))
         m = std::abs(x.imag());
   }
   return m;
}

void evalPlainVariance(std::vector<cx_double> & res, std::vector<cx_double> const& in) {
   evalPlainMul(res, in, in);
   cx_double sum = 0;
   for (auto const& x : res) {
      sum += x;
   }
   for (size_t i = 0; i < res.size(); i++) {
      res[i] = sum/((double)res.size());
   }
}

double maxDiff(std::vector<cx_double> const& in0, std::vector<cx_double> const& in1) {
   size_t len = std::min(in0.size(), in1.size());
   std::vector<cx_double> tmp(len);
   evalPlainNegate(tmp, in1);
   evalPlainAdd(tmp, in0, tmp);
   return largestElm(tmp);
}

double relError(std::vector<cx_double> const& in0, std::vector<cx_double> const& in1) {
   size_t len = std::min(in0.size(), in1.size());
   std::vector<cx_double> diff(len);
   evalPlainNegate(diff, in1);
   evalPlainAdd(diff, diff, in0);
   double res = 0;
   for (size_t i = 0; i < len; i++) {
      double tmp = std::fabs(diff[i].real() / in1[i].real());
      if (res < tmp) {
         res = tmp;
      }
      tmp = std::fabs(diff[i].imag() / in1[i].imag());
      if (res < tmp) {
         res = tmp;
      }
   }
   return res;
}

void randomComplexVector(std::vector<cx_double>& array, size_t n, double rad) {
   if (rad <= 0) {
      rad = 1.0;                // default radius = 1
   }
   array.resize(n);             // allocate space
   for (auto& x : array) {
      long bits = NTL::RandomLen_long(32);         // 32 random bits
      double r = std::sqrt(bits & 0xffff) / 256.0; // sqrt(uniform[0,1])
      double theta =
         2.0L * M_PI * ((bits >> 16) & 0xffff) / 65536.0; // uniform(0,2pi)
      x = std::polar(rad * r, theta);
   }
}

cx_double * randomComplexVector(size_t n, double rad) {
   std::vector<cx_double> vec(n);
   cx_double * pvec = new cx_double[n];
   randomComplexVector(vec, n, rad);
   for (size_t i = 0; i < n; i++) {
      pvec[i] = vec[i];
   }
   return pvec;
}

void randomRealVector(std::vector<cx_double>& array, size_t n, double B) {
   B = fabs(B);
   array.resize(n);             // allocate space
   for (auto& x : array) {
      long bits = NTL::RandomLen_long(32);         // 32 random bits
      double r = std::sqrt(bits & 0xffff) / 256.0; // sqrt(uniform[0,1])
      double sign = ((bits >> 16) & 0xffff) > 32767 ? 1.0 : -1.0;
      x.real(B * r * sign);
      x.imag(0);
   }
}
cx_double * randomRealVector(size_t n, double rad) {
   std::vector<cx_double> vec(n);
   cx_double * pvec = new cx_double[n];
   randomRealVector(vec, n, rad);
   for (size_t i = 0; i < n; i++) {
      pvec[i] = vec[i];
   }
   return pvec;
}

HomomorphicComputation parseHC(char const* v) {
   HomomorphicComputation hc = HC_NOOP;              // Default noop
   if (!v) {
      return hc;
   }
   if (!strcmp(v, "variance")) {
      hc = HC_VARIANCE;
   } else if(!strcmp(v, "sigmoid")) {
      hc = HC_SIGMOID;
   } else if(!strcmp(v, "exp")) {
      hc = HC_EXP;
   }
   return hc;
}




char const* hcString(HomomorphicComputation hc) {
   switch (hc) {
   case HC_VARIANCE : return "variance";
   case HC_SIGMOID  : return "sigmoid";
   case HC_EXP      : return "exp";
   default          : return "noop";
   }
   return "noop";
}

void copyTo(std::complex<double> * dst, std::complex<double> const* src, size_t len) {
   for (size_t i = 0; i < len; i++) {
      dst[i] = src[i];
   }
}

bool isEqual(std::complex<double> const* m0, std::complex<double> const* m1, size_t len) {
   for (size_t i = 0; i < len; i++) {
      if (m0[i] != m1[i]) {
         std::cout.precision(10);
         std::cout << "different @ " << i << " : "
                   << std::scientific << m0[i] << ", " << m1[i] << std::endl;
         return false;
      }
   }
   return true;
}
