#ifndef SEAL_UTILS_H
#define SEAL_UTILS_H
// Some helper functions to use with SEAL

#include <vector>
#include <iostream>
#include <seal/seal.h>
#include <seal/modulus.h>
#include <seal/util/iterator.h>

using namespace seal;

void print_parameters(std::shared_ptr<seal::SEALContext> context);

// compute a^{-1}, where a is a double-CRT polynomial whose evaluation representation
// is in aEvalRep. The double-CRT representation in SEAL is stored as a flat array of
// length coeff_count * modulus_count:
//    [ 0 .. coeff_count-1 , coeff_count .. 2*coeff_count-1, ... ]
//      ^--- a (mod p0)    , ^--- a (mod p1),              , ...
// return if the inverse exists, and result is also in evaluation representation
bool inv_dcrtpoly(util::ConstCoeffIter aEvalRep, std::size_t coeff_count, std::vector<Modulus> const& coeff_modulus,
                  util::CoeffIter result);

// compute a*b, where both a and b are in evaluation representation
void mul_dcrtpoly(util::ConstCoeffIter a, util::ConstCoeffIter b,
                  std::size_t coeff_count, std::vector<Modulus> const& coeff_modulus,
                  util::CoeffIter result);

// compute a+b, where both a and b are in the same representation
void add_dcrtpoly(util::ConstCoeffIter a, util::ConstCoeffIter b,
                  std::size_t coeff_count, std::vector<Modulus> const& coeff_modulus,
                  util::CoeffIter result);

// compute a-b, where both a and b are in the same representation
void sub_dcrtpoly(util::ConstCoeffIter a, util::ConstCoeffIter b,
                  std::size_t coeff_count, std::vector<Modulus> const& coeff_modulus,
                  util::CoeffIter result);

// assign result = a
void assign_dcrtpoly(util::ConstCoeffIter a, std::size_t coeff_count, std::size_t coeff_modulus_count,
                     util::CoeffIter result);

void to_eval_rep(util::CoeffIter a, size_t coeff_count, size_t coeff_modulus_count, util::NTTTables const* small_ntt_tables);

void to_coeff_rep(util::CoeffIter a, size_t coeff_count, size_t coeff_modulus_count, util::NTTTables const* small_ntt_tables);

long double infty_norm(util::ConstCoeffIter a, SEALContext::ContextData const* context_data);

long double l2_norm(util::ConstCoeffIter a, SEALContext::ContextData const* context_data);

std::string poly_to_string(std::uint64_t const* value, EncryptionParameters const& parms);

void print_poly(std::uint64_t const* value, EncryptionParameters const& parms, size_t max_count=0);

template<typename T>
std::ostream &operator <<(std::ostream &os, const std::vector<T> &v) {
   using namespace std;
   os << "[";
   copy(v.begin(), v.end(), ostream_iterator<T>(os, ", "));
   os << "]";
   return os;
}

#endif  // SEAL_UTILS_H
