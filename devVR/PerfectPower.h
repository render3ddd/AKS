#include <iostream>

#include "NTL/ZZ.h" // NTL Libraries

// floor(n^(1/k)) for n >= 1 and k >= 2, using integer Newton iteration only.
// Starting from x0 = 2^ceil(bits(n)/k) >= n^(1/k), the iterates decrease
// monotonically to the floor of the root.
inline NTL::ZZ IntegerRoot(const NTL::ZZ& n, long k) {
    NTL::ZZ x = NTL::power2_ZZ((NTL::NumBits(n) + k - 1) / k);
    for (;;) {
        NTL::ZZ y = ((k - 1) * x + n / NTL::power(x, k - 1)) / k;
        if (y >= x) {
            return x;
        }
        x = y;
    }
}

// function to calculate if n = a^b
// takes input ZZ n and returns 1 if n is a perfect power, 0 otherwise
//
// Previously Newton's method ran in double precision and assigned the root
// to a long, so perfect powers with a base >= 2^63 were missed and some
// inputs never terminated. Every step is now exact. Only prime exponents
// b <= floor(log2 n) need testing: a^(pq) = (a^p)^q.

int PerfectPower(const NTL::ZZ& n){
    if (n < 4) {
        return 0;                       // 1, 2 and 3 are handled by the caller
    }

    const long maxB = NTL::NumBits(n) - 1;  // floor(log2 n): a >= 2 implies b <= log2 n

    for (long b = 2; b <= maxB; b = NTL::NextPrime(b + 1)) {
        NTL::ZZ a = IntegerRoot(n, b);
        if (NTL::power(a, b) == n) {
            // if n is a perfect power
            std::cout << "n = a^b\n" << "b = " << b << "\n" << "a = " << a << "\n\n";
            return 1;
        }
    }

    return 0; // n is not a perfect power
}
