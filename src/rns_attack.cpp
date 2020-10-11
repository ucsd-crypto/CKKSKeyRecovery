#include <Numb.h>
#include <Context.h>
#include <Ciphertext.h>
#include <EvaluatorUtils.h>
#include <Scheme.h>
#include <SchemeAlgo.h>
#include <SecretKey.h>
#include <StringUtils.h>
#include <string.h>
#include "eval.h"
#include "ntl_utils.h"

using namespace std;

// Compute the L-infty norm of a - b, where balanced representation is used
uint64_t maxDiff(uint64_t const* a, uint64_t const* b, int n, uint64_t modQ) {
   uint64_t m = 0;
   uint64_t hQ = modQ / 2;
   int64_t ip = static_cast<int64_t>(modQ); // assume modQ < 2^63
   for (int i=0; i<n; i++) {
      int64_t ia = a[i] < hQ ? static_cast<int64_t>(a[i]) : static_cast<int64_t>(a[i]) - ip;
      int64_t ib = b[i] < hQ ? static_cast<int64_t>(b[i]) : static_cast<int64_t>(b[i]) - ip;
      int64_t d = std::abs(ia - ib);
      if (m < static_cast<uint64_t>(d)) {
         m = static_cast<uint64_t>(d);
      }
   }
   return m;
}

// Homomorphic computation
void evalVariance(Scheme & scheme, Ciphertext & ct, long len, Ciphertext & ctRes) {
   std::cout << "Compute variance" << std::endl;
   ctRes = scheme.mult(ct, ct);             // Enc((m*p)^2)
   scheme.reScaleByAndEqual(ctRes, 1);
   for (int i=2; i<=len; i*=2) {
      Ciphertext tmp = scheme.leftRotateFast(ctRes, len/i);  // tmp = ctRes << len/i
      scheme.addAndEqual(ctRes, tmp);                        // ctRes += tmp
   }
   scheme.multByConstAndEqual(ctRes, 1/(double)len);
   scheme.reScaleByAndEqual(ctRes, 1);
}

void evalSigmoid(Scheme & scheme, Ciphertext & ct, int evalDeg, Ciphertext & ctRes) {
   std::cout << "Compute sigmoid to degree " << evalDeg << std::endl;
   SchemeAlgo algo(scheme);
   ctRes = algo.sigmoid(ct, evalDeg);
}

void evalExp(Scheme & scheme, Ciphertext & ct, int evalDeg, Ciphertext & ctRes) {
   std::cout << "Compute exponent to degree " << evalDeg << std::endl;
   SchemeAlgo algo(scheme);
   ctRes = algo.exponent(ct, evalDeg);
}

// Key Recovery Attack
int attack(HomomorphicComputation hc, long logN, long L, long logp, double ptBound, long evalDeg) {
  // Generate plaintext input
  long slots = 1L << (logN-1);                             // use all slots
  complex<double> *val_input = randomRealVector(slots, ptBound);

  // Key generation
  Context context(logN, logp, L, L+1);
  SecretKey secretKey(context);
  Scheme scheme(secretKey, context);
  if (hc == HC_VARIANCE) {
     // Add rotation keys when computing variance
     for (int i=1; i<=slots/2; i*=2) {
        scheme.addLeftRotKey(secretKey, i); // for left shift 1<<i
     }
  }

  // Encrypt chosen messages to level L
  Ciphertext ctxt_input = scheme.encrypt(val_input, slots, L);

  // Homomorphic computation
  Ciphertext ctxt_res;
  switch (hc) {
  case HC_VARIANCE :
     evalVariance(scheme, ctxt_input, slots, ctxt_res); // compute sum(m^2 for m in val_input)
     break;
  case HC_SIGMOID :
     evalSigmoid(scheme, ctxt_input, evalDeg, ctxt_res); // compute the logistic function
     break;
  case HC_EXP :
     evalExp(scheme, ctxt_input, evalDeg, ctxt_res);
     break;
  default :
     ctxt_res = ctxt_input;
  }
  
  // Approximate decryption
  Plaintext ptxt_res = scheme.decryptMsg(secretKey, ctxt_res); // NTT form
  complex<double>* val_res = scheme.decode(ptxt_res);

  // Check homomorphic computation error
  std::vector<complex<double>> ptIn(val_input, val_input+slots), heRes(val_res, val_res+slots);
  std::vector<complex<double>> ptRes(slots);
  switch (hc) {
  case HC_VARIANCE :
     evalPlainVariance(ptRes, ptIn);
     break;
  case HC_SIGMOID :
     evalPlainFunc(ptRes, ptIn, SpecialFunction::SIGMOID);
     break;
  case HC_EXP :
     evalPlainFunc(ptRes, ptIn, SpecialFunction::EXP);
     break;
  default :
     ptRes.assign(val_input, val_input+slots);
  }
  std::cout << "exact value = " << ptRes[0] << ", computation error = " << abs(ptRes[0] - heRes[0])
            << ", relative error = " << relError(ptRes, heRes) << std::endl;


  // **************************************************
  // Key recovery attack
  // **************************************************

  // Recover the encryption error by encoding the approximate plaintext
  Plaintext ptxt_enc = scheme.encode(val_res, slots, ctxt_res.l);

  // Check for encoding error
  uint64_t* ptxt_res_coeff = new uint64_t[context.N]();
  uint64_t* ptxt_enc_coeff = new uint64_t[context.N]();
  copy(ptxt_res.mx, ptxt_res.mx + context.N, ptxt_res_coeff);
  copy(ptxt_enc.mx, ptxt_enc.mx + context.N, ptxt_enc_coeff);
  context.qiINTTAndEqual(ptxt_res_coeff, 0); // to coeff representation
  context.qiINTTAndEqual(ptxt_enc_coeff, 0);    // to coeff representation, only 1st tower
  uint64_t encoding_error = maxDiff(ptxt_res_coeff, ptxt_enc_coeff, context.N, context.qVec[0]);
  std::cout << "encoding error = " << encoding_error << std::endl;

  // rhs = m' - b
  uint64_t* rRhs = new uint64_t[ctxt_res.l << context.logN]();
  context.sub(rRhs, ptxt_enc.mx, ctxt_res.bx, ctxt_res.l); // What about using ptxt_res.mx ?

  // a^{-1}
  uint64_t* raInv = new uint64_t[ctxt_res.l << context.logN]();
#pragma omp parallel for
  for (int i = 0; i < ctxt_res.l; i++) {
     uint64_t qi = context.qVec[i];
     for (int j = 0; j < context.N; j++) {
        raInv[j + i*context.N] = invMod(ctxt_res.ax[j + i*context.N], qi);
     }
  }

  // guess sk
  uint64_t* rss = new uint64_t[ctxt_res.l << context.logN]();
  context.mul(rss, raInv, rRhs, ctxt_res.l, 0);

  // for debugging, convert to coeff representation
  uint64_t* ss_coeff = new uint64_t[ctxt_res.l << context.logN]();
  uint64_t* sk_coeff = new uint64_t[ctxt_res.l << context.logN]();
  copy(rss, rss + context.N * ctxt_res.l, ss_coeff);
  copy(secretKey.sx, secretKey.sx + context.N * ctxt_res.l, sk_coeff);
  context.INTTAndEqual(ss_coeff, ctxt_res.l, 0);
  context.INTTAndEqual(sk_coeff, ctxt_res.l, 0);

  return !memcmp(rss, secretKey.sx, sizeof(uint64_t) * ctxt_res.l * context.N);
}

int main(int argc, char* argv[]) {
   // Parameters //
   HomomorphicComputation hc = argc>1 ? parseHC(argv[1]) : HC_NOOP;  // Default noop
   long logN         = argc>2 ? atoi(argv[2]) : 15;  // Ring dimension
   long L            = argc>3 ? atoi(argv[3]) : 10;  // Total level of computation
   long logp         = argc>4 ? atoi(argv[4]) : 20;  // 20 bit precision
   double plainBound = argc>5 ? atof(argv[5]) : 1.0; // plaintext size
   long evalDeg      = argc>6 ? atof(argv[6]) : 5;   // max degree to evaluate

   std::cout << "Running attack on " << hcString(hc) << std::endl
             << "  N = " << (1 << logN) << ", L = " << L << ", logp = " << logp << ", |plaintext| = " << plainBound << std::endl;

   srand(time(NULL));
   if (attack(hc, logN, L, logp, plainBound, evalDeg)) {
      cout << "Found key!" << endl;
   }
   else cout << "Attack failed!" << endl;
  
   return 0;
}
