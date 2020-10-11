#include <helib/helib.h>
#include <helib/norms.h>
#include <helib/debugging.h>

#include <NTL/ZZ_pX.h>
#include <NTL/BasicThreadPool.h>
#include <iostream>
#include "helib_utils.h"
#include "eval.h"

// Homomorphic computations
helib::Ctxt evalVariance(helib::EncryptedArrayCx const& ea, helib::Ctxt const& ct, size_t n) {
   std::cout << "Compute variance" << std::endl;

   helib::Ctxt ctRes(ct);              // copy x
   ctRes.multiplyBy(ct);               // x^2
   ctRes.dropSmallAndSpecialPrimes();  // drop moduli p_i added in modUp
   for (int i=2; i<=n; i*=2) {
      helib::Ctxt tmp(ctRes);
      ea.rotate(tmp, n/i);             // tmp = ctRes >> n/i
      tmp.dropSmallAndSpecialPrimes(); // drop moduli p_i added in modUp

      showCtxtScale(tmp, "rotate ");
      ctRes += tmp;
   }
   ctRes.multByConstantCKKS(1/(double)n);
   return ctRes;
}

// ctRes[i] = encryption of x^(2^i), where ct = encryption of x, for 0 <= i <= logDeg
void evalPowerOf2(std::vector<helib::Ctxt*> & ctRes, helib::Ctxt const& ct, int logDeg) {
   ctRes.resize(logDeg+1);
   ctRes[0] = new helib::Ctxt(ct);             // x^(2^0)
   for (int i = 1; i <= logDeg; i++) {
      ctRes[i] = new helib::Ctxt(*ctRes[i-1]); // x^(2^{i-1})
      ctRes[i]->multiplyBy(*ctRes[i-1]);       // x^(2^i)
      ctRes[i]->dropSmallAndSpecialPrimes();
      showCtxtScale(*ctRes[i], "powerOf2 ");
   }
}

// A workaround for multiplying by a constant c, where helib would hit a division by 0 error if c<0
void multByConstantCKKSFix(helib::Ctxt & ct, double c) {
   if (c<0) {
      ct.multByConstantCKKS(-c);
      ct.negate();
   } else {
      ct.multByConstantCKKS(c);
   }
}
   

// Evaluate a polynomial function up to degree evalDeg
helib::Ctxt evalFunction(helib::Ctxt const& ct, size_t n,
                         std::vector<double> const& coeff, int evalDeg = -1) {
   int deg = evalDeg == -1 ? coeff.size()-1 : std::min((size_t)evalDeg, coeff.size() - 1); // assume coeff is not empty
   int logDeg = std::floor(std::log2((double)deg));
   std::cout << "evalFunction " << coeff << " to degree " << deg << std::endl;
   std::vector<helib::Ctxt*> ctPow2s(logDeg+1);
   evalPowerOf2(ctPow2s, ct, logDeg);
   helib::Ctxt ctRes(ct);                    // copy x
   multByConstantCKKSFix(ctRes, coeff[1]);   // c_1 * x
   ctRes.addConstantCKKS(coeff[0]);          // c_1 * x + c_0
   showCtxtScale(ctRes, "c_1 * x + c_0 ");

   for (int i = 2; i <= deg; i++) {
      if (fabs(coeff[i]) < 1e-27) {
         continue;                           // Too small, skip this term
      }
      int k = std::floor(std::log2((double)i));
      int r = i - (1 << k);                  // i = 2^k + r
      helib::Ctxt tmp(*ctPow2s[k]);          // x^(2^k)
      while (r > 0) {
         k = std::floor(std::log2((double)r));
         r = r - (1 << k);
         tmp.multiplyBy(*ctPow2s[k]);
         tmp.dropSmallAndSpecialPrimes();
      }
      multByConstantCKKSFix(tmp, coeff[i]);  // c_i * x^i
      showCtxtScale(tmp,   "c_i * x^i ");
      ctRes += tmp;
      showCtxtScale(ctRes, "add c_i * x^i ");
   }

   for (auto &x : ctPow2s) {
      delete x;
   }
   return ctRes;
}

void test(int logm, int logp, int logQ, double B, int evalDeg, HomomorphicComputation hc) {
   // B is the radius of plaintext numbers
   long m = pow(2,logm);  // Zm*
   long r = logp;         // bit precision
   long L = logQ;         // Number of bits of Q

   // Setup the context
   helib::Context context(m, -1, r); // p = -1 => complex field, ie m = p-1

   // context.scale = 10; // used for sampling error bound
   helib::buildModChain(context, L, /*c=*/2); // 2 columns in key switching key
   helib::SecKey secretKey(context);
   secretKey.GenSecKey();

   if (hc == HC_VARIANCE) {
      helib::addSome1DMatrices(secretKey); // add rotation keys for variance computation
   }
   helib::PubKey publicKey(secretKey);
   helib::EncryptedArrayCx const& ea(context.ea->getCx());
   long n = ea.size();        // # slots

   ea.getPAlgebra().printout();
   std::cout << "r = " << context.alMod.getR() << std::endl;
   std::cout << "ctxtPrimes=" << context.ctxtPrimes
             << ", ciphertext modulus bits=" << context.bitSizeOfQ() << std::endl
             << std::endl;

#ifdef HELIB_DEBUG
   helib::setupDebugGlobals(&secretKey, context.ea);
#endif

   // Initialize the plaintext vector
   std::vector<std::complex<double>> v1, v2; // v1 holds the plaintext input, v2 holds the decryption result
   ea.random(v1,B);                          // generate a random array of size m/2
   std::cout << "v : size = " << v1.size() << ", infty norm = " << largestCxNorm(v1) << std::endl;

   // Encryption
   helib::Ctxt c_v(publicKey);       // Ctxt::parts contains the ciphertext polynomials
   ea.encrypt(c_v, publicKey, v1);

   // Homomorphic computation
   std::vector<double> coeff(11);
   helib::Ctxt c_res(publicKey);
   switch (hc) {
   case HC_VARIANCE :
      c_res = evalVariance(ea, c_v, v1.size());
      break;
   case HC_SIGMOID :
      coeff = SpecialFunction::coeffsOf[SpecialFunction::FuncName::SIGMOID];
      c_res = evalFunction(c_v, n, coeff, evalDeg); // compute the logistic function
      break;
   case HC_EXP :
      coeff = SpecialFunction::coeffsOf[SpecialFunction::FuncName::EXP];
      c_res = evalFunction(c_v, n, coeff, evalDeg); // compute the exponential function
      break;
   default :
      c_res = c_v;              // just copy the input ciphertext
   }
   showCtxtScale(c_res, "result ");
   long logExtraScaling = std::ceil(log2(ea.encodeRoundingError() / 3.5));
   helib::IndexSet s1 = c_res.getPrimeSet();
   while(NTL::log(c_res.getRatFactor())/log(2.0) > r + logExtraScaling + 10 && s1.card() > 1) {
      s1.remove(s1.last());
      c_res.modDownToSet(s1);
      showCtxtScale(c_res, "modDown");
      s1 = c_res.getPrimeSet();
   }

   // Decryption
   ea.decrypt(c_res, secretKey, v2);

   // Check homomorphic computation error
   std::vector<std::complex<double>> ptRes(v2.size());
   switch (hc) {
   case HC_VARIANCE :
      evalPlainVariance(ptRes, v1);
      break;
   case HC_SIGMOID :
      evalPlainFunc(ptRes, v1, SpecialFunction::SIGMOID, evalDeg);
      break;
   case HC_EXP :
      evalPlainFunc(ptRes, v1, SpecialFunction::EXP, evalDeg);
      break;
   default :
      ptRes = v1;
   }
   std::cout << "computation error = " << maxDiff(ptRes, v2) // abs(ptRes[0] - v2[v2.size()-1])
             << ", relative error = " << relError(v2, ptRes) << std::endl; // maxDiff(ptRes, v2)/largestElm(ptRes)

   // Key recovery attack ************************************************** //

   // Now let's try to recover sk
   NTL::xdouble scalingFactor = c_res.getRatFactor();

   // Here we use a modified encoding function to round directly into ZZX,
   // instead of rounding to a helib::zzX, which is a vector of long so it could
   // cause integer overflow
   NTL::ZZX mPrimeX;
   helib::CKKS_embedInSlots(mPrimeX, v2, context.zMStar, NTL::to_double(scalingFactor));

   // Check if encoding recovers the decrypted ptxt (before modulo reduction)
   NTL::ZZX encodingDiffX = mPrimeX - helib::decrypted_ptxt_;
   NTL::xdouble mPrimeNorm = helib::coeffsL2Norm(mPrimeX);
   std::cout << "encoding error = " << helib::largestCoeff(encodingDiffX) << std::endl;
   std::cout << "m' norm = " << mPrimeNorm << ", bits = " << NTL::log(mPrimeNorm)/std::log(2) << std::endl;
   
   helib::DoubleCRT ctxtb = c_res.getPart(0); // Ctxt::getPart() is added to helib to access
   helib::DoubleCRT ctxta = c_res.getPart(1); // the individual parts

   NTL::ZZX ctxtbX, ctxtaX;
   ctxtb.toPoly(ctxtbX,true);
   ctxta.toPoly(ctxtaX,true);

   NTL::ZZ Q;
   context.productOfPrimes(Q, c_res.getPrimeSet());
   std::cout << "sk.Q = " << Q << std::endl;

   NTL::ZZ_p::init(Q);

   NTL::ZZ_pX ctxtb_pX, ctxta_pX, mPrime_pX, phim_pX;
   NTL::conv(ctxtb_pX, ctxtbX);
   NTL::conv(ctxta_pX, ctxtaX);
   NTL::conv(mPrime_pX, mPrimeX);
   NTL::conv(phim_pX, context.zMStar.getPhimX());

   NTL::ZZ_pX ss_pX, ctxtaInv_pX;
   NTL::ZZ_pX c_pX = mPrime_pX - ctxtb_pX; // c = m' - cb = ca * s
   NTL::InvMod(ctxtaInv_pX, ctxta_pX, phim_pX);    // ca^{-1} mod (X^{m/2} + 1)
   NTL::MulMod(ss_pX, c_pX, ctxtaInv_pX, phim_pX); // c * ca^{-1} = s

   helib::DoubleCRT sk = secretKey.sKeys[0];
   NTL::ZZX skX;
   sk.toPoly(skX,true);
   NTL::ZZ_pX sk_pX;
   NTL::conv(sk_pX, skX);

   bool foundKey = (ss_pX == sk_pX);
   std::cout << (foundKey ? "Found key!" : "Attack failed") << std::endl;
}


int main(int argc, char * argv[]) {
   HomomorphicComputation hc = argc>1 ? parseHC(argv[1]) : HC_NOOP;  // Default noop
   long logQ = 300;
   long logm = argc>2 ? atoi(argv[2])+1 : 17;    // m = 2N
   long logp = argc>3 ? atoi(argv[3])   : 20;    // 20 bit precision
   double plainBound = argc>4 ? atof(argv[4]) : 1.0; // plaintext size
   long evalDeg = argc>5 ? atoi(argv[5]) : -1;       // default to all degrees

   std::cout << "Running helib attack for " << hcString(hc)
             << ", N = 2^" << logm-1
             << ", logp = " << logp
             << ", |plaintext| = " << plainBound
             << ", evalDeg = " << evalDeg << std::endl;
   
   NTL::SetNumThreads(8);
   test(/*logm=*/logm,
        /*logp=*/logp,
        /*logQ=*/logQ,
        /*B=*/plainBound,
        /*evalDeg=*/evalDeg,
        /*hc=*/hc);
   return 0;
}
