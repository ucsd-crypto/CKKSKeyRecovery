#include "palisade_utils.h"

#include <iostream>

void printNativePoly(NativePoly const& p, size_t num) {
   NativePoly pc = p;
   pc.SetFormat(COEFFICIENT);
   std::cout << "Poly [" << pc.GetParams()->GetModulus() << "] " << std::endl;
   NativePoly::Vector const& vi = pc.GetValues();
   for (size_t j = 0; j < vi.GetLength() && j < num; j++) {
      std::cout << vi[j] << ", ";
   }
   std::cout << std::endl;
}

void printDCRTPoly(DCRTPoly const& p, size_t num) {
   DCRTPoly pc = p;
   pc.SetFormat(COEFFICIENT);
   auto params = pc.GetParams();
   std::cout << "DCRTPoly [" << params->GetModulus() << "] " << std::endl;
   for (size_t i = 0; i < params->GetParams().size(); i++) {
      auto modParam = params->GetParams()[i];
      std::cout << " * [" << modParam->GetModulus() << "] ";
      DCRTPoly::PolyType const& pi = pc.GetElementAtIndex(i);
      DCRTPoly::PolyType::Vector const& vi = pi.GetValues();
      for (size_t j = 0; j < vi.GetLength() && j < num; j++) {
         std::cout << vi[j] << ", ";
      }
      std::cout << std::endl;
   }
}

void printPoly(Poly const& p, size_t num) {
   Poly pc = p;
   pc.SetFormat(COEFFICIENT);
   std::cout << "Poly [" << pc.GetParams()->GetModulus() << "] " << std::endl;
   Poly::Vector const& vi = pc.GetValues();
   for (size_t j = 0; j < vi.GetLength() && j < num; j++) {
      std::cout << vi[j] << ", ";
   }
   std::cout << std::endl;
}


