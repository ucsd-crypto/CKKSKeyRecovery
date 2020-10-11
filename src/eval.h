#ifndef EVAL_H
#define EVAL_H
// Evaluation of functions on plaintext numbers

#include <complex>
#include <vector>
#include <map>

enum HomomorphicComputation {
   HC_NOOP,
   HC_VARIANCE,
   HC_SIGMOID,
   HC_EXP
};
char const* hcString(HomomorphicComputation hc);
HomomorphicComputation parseHC(char const* v);

typedef std::complex<double> cx_double;

void evalPlainNegate(std::vector<cx_double> & res, std::vector<cx_double> const& in);
void evalPlainInverse(std::vector<cx_double> & res, std::vector<cx_double> const& in);

void evalPlainAdd(std::vector<cx_double> & res,
                  std::vector<cx_double> const& in0, std::vector<cx_double> const& in1);
void evalPlainMul(std::vector<cx_double> & res,
                  std::vector<cx_double> const& in0, std::vector<cx_double> const& in1);

void evalPlainPowerOf2(std::vector<cx_double> & res, std::vector<cx_double> const& in, size_t logDeg);
void evalPlainPower(std::vector<cx_double> & res, std::vector<cx_double> const& in, size_t deg);

void evalPlainAddi(std::vector<cx_double> & res, std::vector<cx_double> const& in, double c);
void evalPlainMuli(std::vector<cx_double> & res, std::vector<cx_double> const& in, double c);

struct SpecialFunction {
   enum FuncName {
      LOG,
      EXP,
      SIGMOID
   };
   static std::map<FuncName, std::vector<double>> coeffsOf;
   // each coefficient vector is indexed from 0, where the entry at index i is the
   // coefficient of X^i
};

void evalPlainFunc(std::vector<cx_double> & res, std::vector<cx_double> const& in, std::vector<double> const& coeff, int deg = -1);
void evalPlainFunc(std::vector<cx_double> & res, std::vector<cx_double> const& in, SpecialFunction::FuncName name, int deg = -1);
void evalPlainFunc(std::vector<cx_double> & res, cx_double * in, size_t len, SpecialFunction::FuncName name, int deg = -1);

double largestElm(std::vector<cx_double> const& vec);
double maxDiff(std::vector<cx_double> const& in0, std::vector<cx_double> const& in1);
double relError(std::vector<cx_double> const& in0, std::vector<cx_double> const& in1);

void evalPlainVariance(std::vector<cx_double> & res, std::vector<cx_double> const& in);

void randomComplexVector(std::vector<cx_double>& array, size_t n, double rad = 1.0);
cx_double * randomComplexVector(size_t n, double rad = 1.0);
void randomRealVector(std::vector<cx_double>& array, size_t n, double B = 1.0);
cx_double * randomRealVector(size_t n, double B = 1.0);

#endif  // EVAL_H




                   
