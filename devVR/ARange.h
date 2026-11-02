#pragma once

#include <cmath>

#include "NTL/ZZ.h"

// How far the congruence tests (x - a)^n = x^n - a (mod n, x^r - 1) run.
//
//   Theorem        : 1 <= a <= floor(sqrt(phi(r)) * log2 n)   (Lenstra / AKS bound)
//   PhiMinusOne    : 1 <= a <= phi(r) - 1                     (dissertation, Algorithm 2.4.2)
//   LambdaMinusOne : 1 <= a <= lambda(r) - 1                  (Egan, 2005; dissertation 5.2.4)
//
// Only Theorem is backed by the AKS/Lenstra proof. The other two modes are kept
// to reproduce the dissertation's experiments; lambda(r) - 1 can fall below the
// theorem's bound, so a PRIME verdict in that mode is NOT certified.
enum class ARange { Theorem, PhiMinusOne, LambdaMinusOne };

inline const char* ARangeName(ARange mode) {
    switch (mode) {
        case ARange::Theorem:        return "theorem: floor(sqrt(phi(r)) log2 n)";
        case ARange::PhiMinusOne:    return "phi(r) - 1 [not certified]";
        case ARange::LambdaMinusOne: return "lambda(r) - 1 [not certified]";
    }
    return "?";
}

// floor(sqrt(phi) * log2 n), capped at n - 1 (a must stay below n).
inline long TheoremABound(const NTL::ZZ& n, long phi) {
    double log2n = NTL::log(n) / std::log(2.0);   // NTL::log(ZZ) is the natural log
    long bound = static_cast<long>(std::floor(std::sqrt(static_cast<double>(phi)) * log2n));
    if (n <= bound) {
        bound = NTL::to_long(n - 1);
    }
    return bound;
}

inline long ABound(ARange mode, const NTL::ZZ& n, long phi, long lambda) {
    switch (mode) {
        case ARange::Theorem:        return TheoremABound(n, phi);
        case ARange::PhiMinusOne:    return phi - 1;
        case ARange::LambdaMinusOne: return lambda - 1;
    }
    return TheoremABound(n, phi);
}
