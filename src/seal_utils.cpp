#include "seal_utils.h"
#include <seal/util/uintarithsmallmod.h>
#include <seal/util/polyarithsmallmod.h>

#include <iostream>

void print_parameters(std::shared_ptr<seal::SEALContext> context) {
    auto &context_data = *context->key_context_data();
    std::cout << "Encryption parameters :" << std::endl;
    std::cout << "   poly_modulus_degree: " <<
        context_data.parms().poly_modulus_degree() << std::endl;
    // Print the size of the true (product) coefficient modulus.
    std::cout << "   coeff_modulus size: ";

    std::cout << context_data.total_coeff_modulus_bit_count() << " (";
    auto coeff_modulus = context_data.parms().coeff_modulus();
    std::size_t coeff_mod_count = coeff_modulus.size();
    for (std::size_t i = 0; i < coeff_mod_count - 1; i++)
    {
        std::cout << coeff_modulus[i].bit_count() << " + ";
    }

    std::cout << coeff_modulus.back().bit_count();
    std::cout << ") bits" << std::endl;
}

bool inv_dcrtpoly(util::ConstCoeffIter operand, std::size_t coeff_count, std::vector<Modulus> const& coeff_modulus,
                  util::CoeffIter result) {
   bool * has_inv = new bool[coeff_modulus.size()];
   std::fill_n(has_inv, coeff_modulus.size(), true);
#pragma omp parallel for
   for (size_t j = 0; j < coeff_modulus.size(); j++) {
      for (size_t i = 0; i < coeff_count && has_inv[j]; i++) {
         uint64_t inv = 0;
         if (util::try_invert_uint_mod(operand[i + (j * coeff_count)], coeff_modulus[j], inv)) {
            result[i + (j * coeff_count)] = inv;
         } else {
            has_inv[j] = false;
         }
      }
   }
   for (size_t j = 0; j < coeff_modulus.size(); j++) {
      if (!has_inv[j]) return false;
   }
   delete [] has_inv;
   return true;
}

void mul_dcrtpoly(util::ConstCoeffIter a, util::ConstCoeffIter b, std::size_t coeff_count,
                  std::vector<Modulus> const& coeff_modulus, util::CoeffIter result) {
#pragma omp parallel for
   for (size_t j = 0; j < coeff_modulus.size(); j++) {
      util::dyadic_product_coeffmod(a + (j * coeff_count),
                                    b + (j * coeff_count),
                                    coeff_count,
                                    coeff_modulus[j],
                                    result + (j * coeff_count));
   }
}

void add_dcrtpoly(util::ConstCoeffIter a, util::ConstCoeffIter b, std::size_t coeff_count,
                  std::vector<Modulus> const& coeff_modulus, util::CoeffIter result) {
#pragma omp parallel for
   for (size_t j = 0; j < coeff_modulus.size(); j++) {
      util::add_poly_coeffmod(a + (j * coeff_count),
                              b + (j * coeff_count),
                              coeff_count,
                              coeff_modulus[j],
                              result + (j * coeff_count));
   }
}

void sub_dcrtpoly(util::ConstCoeffIter a, util::ConstCoeffIter b, std::size_t coeff_count,
                  std::vector<Modulus> const& coeff_modulus, util::CoeffIter result) {
#pragma omp parallel for
   for (size_t j = 0; j < coeff_modulus.size(); j++) {
      util::sub_poly_coeffmod(a + (j * coeff_count),
                              b + (j * coeff_count),
                              coeff_count,
                              coeff_modulus[j],
                              result + (j * coeff_count));
   }
}

void assign_dcrtpoly(util::ConstCoeffIter a, std::size_t coeff_count, std::size_t coeff_modulus_count,
                     util::CoeffIter result) {
#pragma omp parallel for
   for (size_t i = 0; i < coeff_modulus_count; i++) {
      util::set_poly(a + (i * coeff_count), coeff_count, 1, result + (i * coeff_count));
   }
}

void to_eval_rep(util::CoeffIter a, size_t coeff_count, size_t coeff_modulus_count, util::NTTTables const* small_ntt_tables) {
#pragma omp parallel for
   for (size_t j = 0; j < coeff_modulus_count; j++) {
      util::ntt_negacyclic_harvey(a + (j * coeff_count), small_ntt_tables[j]); // ntt form
   }
}

void to_coeff_rep(util::CoeffIter a, size_t coeff_count, size_t coeff_modulus_count, util::NTTTables const* small_ntt_tables) {
#pragma omp parallel for
   for (size_t j = 0; j < coeff_modulus_count; j++) {
      util::inverse_ntt_negacyclic_harvey(a + (j * coeff_count), small_ntt_tables[j]); // non-ntt form
   }
}

long double infty_norm(util::ConstCoeffIter a, SEALContext::ContextData const* context_data) {
   auto &ciphertext_parms = context_data->parms();
   auto &coeff_modulus = ciphertext_parms.coeff_modulus();
   size_t coeff_mod_count = coeff_modulus.size();
   size_t coeff_count = ciphertext_parms.poly_modulus_degree();
   auto decryption_modulus = context_data->total_coeff_modulus();
   auto upper_half_threshold = context_data->upper_half_threshold();

   long double max = 0;

   auto aCopy(util::allocate_zero_poly(coeff_count, coeff_mod_count, MemoryManager::GetPool()));
   assign_dcrtpoly(a, coeff_count, coeff_mod_count, aCopy.get());

   // CRT-compose the polynomial
   context_data->rns_tool()->base_q()->compose_array(aCopy.get(), coeff_count, MemoryManager::GetPool());

   long double two_pow_64 = powl(2.0, 64);

   for (std::size_t i = 0; i < coeff_count; i++) {
      long double coeff = 0.0, cur_pow = 1.0;
      if (util::is_greater_than_or_equal_uint(aCopy.get() + (i * coeff_mod_count),
                                              upper_half_threshold, coeff_mod_count)) {
         for (std::size_t j = 0; j < coeff_mod_count; j++, cur_pow *= two_pow_64) {
            if (aCopy[i * coeff_mod_count + j] > decryption_modulus[j]) {
               auto diff = aCopy[i * coeff_mod_count + j] - decryption_modulus[j];
               coeff += diff ? static_cast<long double>(diff) * cur_pow : 0.0;
            } else {
               auto diff = decryption_modulus[j] - aCopy[i * coeff_mod_count + j];
               coeff -= diff ? static_cast<long double>(diff) * cur_pow : 0.0;
            }
         }
      } else {
         for (std::size_t j = 0; j < coeff_mod_count; j++, cur_pow *= two_pow_64) {
            auto curr_coeff = aCopy[i * coeff_mod_count + j];
            coeff += curr_coeff ? static_cast<long double>(curr_coeff) * cur_pow : 0.0;
         }
      }

      if (fabsl(coeff) > max) {
         max = fabsl(coeff);
      }
   }

   return max;
}

long double l2_norm(util::ConstCoeffIter a, SEALContext::ContextData const* context_data) {
   auto &ciphertext_parms = context_data->parms();
   auto &coeff_modulus = ciphertext_parms.coeff_modulus();
   size_t coeff_mod_count = coeff_modulus.size();
   size_t coeff_count = ciphertext_parms.poly_modulus_degree();
   auto decryption_modulus = context_data->total_coeff_modulus();
   auto upper_half_threshold = context_data->upper_half_threshold();

   long double sum = 0;

   auto aCopy(util::allocate_zero_poly(coeff_count, coeff_mod_count, MemoryManager::GetPool()));
   assign_dcrtpoly(a, coeff_count, coeff_mod_count, aCopy.get());

   // CRT-compose the polynomial
   context_data->rns_tool()->base_q()->compose_array(aCopy.get(), coeff_count, MemoryManager::GetPool());

   long double two_pow_64 = powl(2.0, 64);

   for (std::size_t i = 0; i < coeff_count; i++) {
      long double coeff = 0.0, cur_pow = 1.0;
      if (util::is_greater_than_or_equal_uint(aCopy.get() + (i * coeff_mod_count),
                                              upper_half_threshold, coeff_mod_count)) {
         for (std::size_t j = 0; j < coeff_mod_count; j++, cur_pow *= two_pow_64) {
            if (aCopy[i * coeff_mod_count + j] > decryption_modulus[j]) {
               auto diff = aCopy[i * coeff_mod_count + j] - decryption_modulus[j];
               coeff += diff ? static_cast<long double>(diff) * cur_pow : 0.0;
            } else {
               auto diff = decryption_modulus[j] - aCopy[i * coeff_mod_count + j];
               coeff -= diff ? static_cast<long double>(diff) * cur_pow : 0.0;
            }
         }
      } else {
         for (std::size_t j = 0; j < coeff_mod_count; j++, cur_pow *= two_pow_64) {
            auto curr_coeff = aCopy[i * coeff_mod_count + j];
            coeff += curr_coeff ? static_cast<long double>(curr_coeff) * cur_pow : 0.0;
         }
      }

      sum += coeff * coeff;
   }

   return sqrtl(sum);
}

std::string poly_to_string(std::uint64_t const* value, EncryptionParameters const& parms) {
   auto coeff_modulus = parms.coeff_modulus();
   size_t coeff_mod_count = coeff_modulus.size();
   size_t coeff_count = parms.poly_modulus_degree();
   std::ostringstream result;
   for (size_t i = 0; i < coeff_mod_count; i++) {
      auto mod = coeff_modulus[i].value();
      if (i>0) {
         result << std::endl;
      }
      result << "[" << mod << "]: ";
      for (size_t j = 0; j < coeff_count; j++) {
         std::uint64_t v = *value;
         if (v >= mod/2) {
            result << "-" << mod-v;
         } else {
            result << v;
         }
         result << (j==coeff_count?"":", ");
         value++;
      }
   }
   return result.str();
}


void print_poly(std::uint64_t const* value, EncryptionParameters const& parms, size_t max_count) {
   auto coeff_modulus = parms.coeff_modulus();
   size_t coeff_mod_count = coeff_modulus.size();
   size_t coeff_count = parms.poly_modulus_degree();
   for (size_t i = 0; i < coeff_mod_count; i++) {
      auto mod = coeff_modulus[i].value();
      std::uint64_t const* v = value + i*coeff_count;
      if (i>0) {
         std::cout << std::endl;
      }
      std::cout << "[" << mod << "]: ";
      for (size_t j = 0; j < coeff_count && (max_count == 0 || j < max_count); j++) {
         if (*v >= mod/2) {
            std::cout << "-" << mod-(*v);
         } else {
            std::cout << *v;
         }
         std::cout << (j==coeff_count?"":", ");
         v++;
      }
   }
   std::cout.flush();
}
