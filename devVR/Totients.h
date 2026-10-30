#pragma once

// Euler's phi(r) and Carmichael's lambda(r) from a single factorisation of r
// (REFACTOR 2.5). Trial division only up to sqrt(r); r is small in AKS.

struct Totients {
    long phi;
    long lambda;
};

inline long GcdLong(long a, long b) {
    while (b != 0) {
        long t = a % b;
        a = b;
        b = t;
    }
    return a;
}

inline Totients PhiAndLambda(long r) {
    Totients t = {1, 1};
    if (r <= 1) {
        return t;                                   // phi(1) = lambda(1) = 1
    }

    // Each prime power p^k exactly dividing r contributes
    //   phi(p^k)    = p^(k-1) (p - 1)
    //   lambda(p^k) = phi(p^k), except lambda(2^k) = 2^(k-2) for k >= 3
    // and lambda(r) is the lcm of the lambda(p^k).
    auto addPrimePower = [&t](long p, long k) {
        long pk1 = 1;                               // p^(k-1)
        for (long i = 1; i < k; ++i) {
            pk1 *= p;
        }
        long phiPk = pk1 * (p - 1);
        long lambdaPk = (p == 2 && k >= 3) ? phiPk / 2 : phiPk;
        t.phi *= phiPk;
        t.lambda = t.lambda / GcdLong(t.lambda, lambdaPk) * lambdaPk;
    };

    for (long p = 2; p * p <= r; ++p) {
        if (r % p == 0) {
            long k = 0;
            while (r % p == 0) {
                r /= p;
                ++k;
            }
            addPrimePower(p, k);
        }
    }
    if (r > 1) {
        addPrimePower(r, 1);                        // remaining prime factor > sqrt(r)
    }
    return t;
}
