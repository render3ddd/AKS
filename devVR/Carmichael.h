#include <vector>
#include <algorithm>
#include <cstdio>

#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"
NTL_CLIENT

// Carmichael's reduced totient lambda(n) and helpers.
//
// REFACTOR 1.8: trial division now stops at sqrt(n) (it ran to n/2), only
// the actual prime factors are recorded (every prime <= n/2 used to be
// pushed with exponent 0), the LCM is taken with gcds instead of matching
// factor lists by position, isPrime(0/1) is false, and subLambda always
// initialises its result.

bool isPrime(long n) {
    if (n < 2) {
        return false;
    }
    for (long d = 2; d * d <= n; ++d) {
        if (n % d == 0) {
            return false;
        }
    }
    return true;
}

// Prime factorisation n = prod p[i]^e[i] (p ascending), by trial division up to sqrt(n).
void primeFactors(long n, vector<long>& p, vector<long>& e) {
    for (long d = 2; d * d <= n; ++d) {
        if (n % d == 0) {
            long k = 0;
            while (n % d == 0) {
                n /= d;
                ++k;
            }
            p.push_back(d);
            e.push_back(k);
        }
    }
    if (n > 1) {                        // what is left is a prime factor > sqrt(original n)
        p.push_back(n);
        e.push_back(1);
    }
}

static unsigned long long gcdULL(unsigned long long a, unsigned long long b) {
    while (b != 0) {
        unsigned long long t = a % b;
        a = b;
        b = t;
    }
    return a;
}

unsigned long long LCM(vector<long>& n) {
    unsigned long long lcm = 1;
    for (long v : n) {
        unsigned long long u = static_cast<unsigned long long>(v);
        lcm = lcm / gcdULL(lcm, u) * u;
    }
    return lcm;
}

unsigned long long Phi(unsigned long long r){
    unsigned long long eu = 1;

    for (unsigned long long p = 2; p * p <= r; p += 2) {
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

    unsigned long long eu1 = eu;

    // now r is prime or 1
    if (r == 1) {
        return eu1;
    }
    else {
        return eu1 * (r - 1) ;
    }
}

// lambda(a^b) for a prime a and b >= 0.
long subLambda(const long& a, const long& b) {

    unsigned long long lambda = 1;          // lambda(a^0) = lambda(1) = 1

    if (b == 0) {
        lambda = 1;
    }
    else if (b <= 2 || a >= 3) {
        lambda = Phi(power_long(a,b));
    }
    else {                                  // a == 2 and b >= 3
        lambda = Phi(power_long(a,b))/2;
    }

    return lambda;
}

ZZ Carmichael(const long& n) {

    if (n <= 1) {
        return ZZ(1);
    }

    std::vector<long> p;
    std::vector<long> e;

    primeFactors(n,p,e);

    std::vector<long> comp;

    for (long i = 0; i < static_cast<long>(p.size()); ++i) {
        comp.push_back(subLambda(p[i],e[i]));
    }

    return ZZ(LCM(comp));
}
