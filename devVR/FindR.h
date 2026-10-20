#pragma once

#include "NTL/ZZ.h"

// Returns true if the multiplicative order of a modulo r is greater than k,
// i.e. a^i != 1 (mod r) for every 1 <= i <= k. Requires 0 <= a < r and
// gcd(a, r) = 1. Uses word-sized modular arithmetic: r is small in AKS.
inline bool OrderExceeds(long a, long r, long k) {
    long x = 1;
    for (long i = 1; i <= k; ++i) {
        x = NTL::MulMod(x, a, r);
        if (x == 1) {
            return false;
        }
    }
    return true;
}
