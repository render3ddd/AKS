#pragma once

#include "NTL/ZZ.h"

// Euler's phi function

inline NTL::ZZ Euler(long r){
    long eu = 1;

    for (long p = 2; p * p <= r; p += 2) {
        if (r % p == 0) {
            eu *= p - 1;
            r /= p;

            while (r % p == 0) {
                eu *= p;
                r /= p;
            }
        }

        if(p == 2) {
            --p;
        }
    }

    NTL::ZZ eu1 = NTL::to_ZZ(eu);

    // now r is prime or 1
    if (r == 1) {
        return eu1;
    }
    else {
        return eu1 * (r - 1) ;
    }
}