#include "HEAAN.h"
#include "StringUtils.h"
#include <NTL/ZZX.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ.h>
#include <NTL/GF2X.h>
#include <NTL/GF2E.h>
#include <NTL/GF2.h>
#include "eval.h"
#include "ntl_utils.h"

using namespace std;
using namespace NTL;

// Convert vector of integers to NTL polynomial
void NTLpolyQ(ZZ_pE & pX, const ZZ* vecZ, const int n) {
  ZZX zX;
  for (int i=0; i<n; i++) {
    SetCoeff(zX,i,vecZ[i]);
  }
  pX = conv<ZZ_pE>(conv<ZZ_pX>(zX));
}

#define Conv2toZ(p) (conv<ZZX>(conv<GF2X>(p)))
#define ConvZtoQ(p) (conv<ZZ_pE>(conv<ZZ_pX>(p)))
#define Conv2toQ(p) (ConvZtoQ(Conv2toZ(p)))
#define ConvQtoZ(q) (conv<ZZX>(conv<ZZ_pX>(q)))
#define ConvZto2(q) (conv<GF2E>(conv<GF2X>(q)))
#define ConvQto2(q) (ConvZto2(ConvQtoZ(q)))
#define ConvQto2X(q) (conv<GF2X>(ConvQtoZ(q)))

// Compute the L-infty norm of aX - bX
ZZ maxDiff(ZZ const* aX, ZZ const* bX, int n, ZZ const& modQ) {
   ZZ m = ZZ::zero();
   ZZ hQ = modQ / 2;
   for (int i=0; i<n; i++) {
      ZZ d = (aX[i] - bX[i]) % modQ;
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

// Homomorphic computation
void evalVariance(Scheme & scheme, Ciphertext & ct, long logp, long len, Ciphertext & ctRes) {
   std::cout << "Compute variance" << std::endl;
   Ciphertext ctMult;
   scheme.mult(ctMult, ct, ct);             // Enc((m*p)^2)
   scheme.reScaleBy(ctRes, ctMult, logp);
   for (int i=2; i<=len; i*=2) {
      Ciphertext tmp;
      scheme.leftRotateFast(tmp, ctRes, len/i);  // tmp = ciphertext << len/i
      scheme.addAndEqual(ctRes, tmp);            // ciphertext += tmp
   }
   scheme.multByConstAndEqual(ctRes, 1/(double)len, logp);
   scheme.reScaleBy(ctRes, ctRes, logp);
}

void evalSigmoid(Scheme & scheme, Ciphertext & ct, long logp, int evalDeg, Ciphertext & ctRes) {
   std::cout << "Compute sigmoid" << std::endl;
   SchemeAlgo algo(scheme);
   algo.function(ctRes, ct, SIGMOID, logp, evalDeg);
}

void evalExp(Scheme & scheme, Ciphertext & ct, long logp, int evalDeg, Ciphertext & ctRes) {
   std::cout << "Compute exponent" << std::endl;
   SchemeAlgo algo(scheme);
   algo.function(ctRes, ct, EXPONENT, logp, evalDeg);
}

// Key Recovery Attack
int attack(Ring ring, long n, long logp, long logq, double ptBound, int evalDeg, HomomorphicComputation hc) {
  // Generate plaintext input
  complex<double> *val_input = randomComplexVector(n, ptBound); // use all slots

  // Key generation
  SecretKey secretKey(ring);
  Scheme scheme(secretKey, ring);

  if (hc == HC_VARIANCE) {
     // Add rotation keys when computing variance
     for (int i=1; i<=n/2; i*=2) {
        scheme.addLeftRotKey(secretKey, i); // for left shift 1<<i
     }
  }

  // Encrypt chosen messages, until ciphertext if found with ax invertible mod 2
  Ciphertext ctxt_res, ctxt_input;
  int aParity = 0;
  while (!aParity) {
    scheme.encrypt(ctxt_input, val_input, n, logp, logq); // Enc(m*p) for scaling factor p

    switch (hc) {
    case HC_VARIANCE :
       evalVariance(scheme, ctxt_input, logp, n, ctxt_res); // compute sum(m^2 for m in val_input)
       break;
    case HC_SIGMOID :
       evalSigmoid(scheme, ctxt_input, logp, evalDeg, ctxt_res);
       break;
    case HC_EXP :
       evalExp(scheme, ctxt_input, logp, evalDeg, ctxt_res);
       break;
    default :
       ctxt_res.copy(ctxt_input); // noop, just copy the input ciphertext
    }

    for (int i=0; i<n; i++)
      aParity = (aParity + ctxt_res.ax[i]) % 2;

    std::cout << "parity = " << aParity << std::endl;
  }
  
  // Approximate decryption
  Plaintext ptxt_res;
  scheme.decryptMsg(ptxt_res, secretKey, ctxt_res);
  complex<double>* val_res = scheme.decode(ptxt_res);

  // Check homomorphic computation error
  std::vector<complex<double>> heRes(val_res, val_res+n), ptIn(val_input, val_input+n), ptRes(n);
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
     ptRes.assign(val_input, val_input+n); // noop, just copy the plaintext input
  }
  std::cout << "computation error = " << maxDiff(ptRes, heRes)
            << ", relative error = " << relError(ptRes, heRes) << std::endl;

  // Recover the encryption error by encoding the approximate plaintext
  Plaintext error;
  ring.encode(error.mx, val_res, n, ctxt_res.logp);

  // Check for encoding error
  ZZ ql = power2_ZZ(ctxt_res.logq); // The modulus of current level (in ctxt_res)
  std::cout << "encoding error = " << maxDiff(ptxt_res.mx, error.mx, n, ql) << std::endl
            << "m' norm = " << maxElm(ptxt_res.mx, n, ql) << std::endl;
  std::cout << "q_l = " << ql << std::endl;

  // Initialize rings ZZ_pE = ZZ[X] / (ql,1+X^N) and GF2X = ZZ[X] / (2,1+X^N)
  ZZ_p::init(ql);
  ZZX mX = ZZX(INIT_MONO,N,1) + 1; // X^N + 1
  ZZ_pE::init(conv<ZZ_pX>(mX));

  GF2X mX2 = conv<GF2X>(mX);
  GF2E::init(mX2);
  
  ZZ_pE aX, eX, sX, bX;

  NTLpolyQ(aX,ctxt_res.ax,N);
  NTLpolyQ(bX,ctxt_res.bx,N);
  NTLpolyQ(eX,error.mx,N);
  NTLpolyQ(sX,secretKey.sx,N);

  GF2E aX2 = ConvQto2(aX);
  ZZ_pE s0X = Conv2toQ(ConvQto2(eX-bX)/aX2);
  ZZX hX = ConvQtoZ(eX - bX + aX*s0X) / 2;
  ZZ_pE s1X = Conv2toQ(ConvZto2(hX) / aX2);
  ZZ_pE ssX = 2*s1X - s0X;
  return (ssX == sX);
}


int main(int argc, char* argv[]) {
  // Parameters //
  long logn = logNh; // use all available slots
  long n = 1 << logn;
  long logq = 300;
  HomomorphicComputation hc = argc>2 ? parseHC(argv[2]) : HC_NOOP;
  long logp = argc>3 ? atoi(argv[3]) : 20;          // 20 bit precision
  double plainBound = argc>4 ? atof(argv[4]) : 1.0; // plaintext size
  long evalDeg = argc>5 ? atoi(argv[5]) : 3;        // max degree to evaluate

  long numThread = 8; 
  SetNumThreads(numThread);
  int iter = 1;
  if (argc == 2) {
    iter = atoi(argv[1]);
  }
  
  // Key Generation
  Ring ring;

  int success=0;
  for (int i=0; i<iter; i++) {
     std::cout << "Running attack on " << hcString(hc) << std::endl
               << " N = " << N << ", logq = " << logq << ", logp = " << logp << ", |plaintext| = " << plainBound << std::endl;

    srand(time(NULL));
    if (attack(ring, n, logp, logq, plainBound, evalDeg, hc)) {
      cout << "Found key!" << endl;
      success++;
    }
    else cout << "Attack failed!" << endl;
  }
  
  cout << "Attack worked " << success << " times out of " << iter << endl;
  return 0;
}
