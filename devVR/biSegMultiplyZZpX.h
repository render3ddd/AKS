
#include <cmath>
#include <vector>
#include <algorithm>

#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
NTL_CLIENT

ZZ_p getMaxCoeff(const ZZ_pX& f) {
    /*
        Returns the largest integer coefficient, f_i, of the input
        polynomial, f(x) (mod p), in O(D) time where D is the number
        of terms in f(x) (mod p), including those with coefficients
        equal to zero.
    */

    ZZ_p f_i = ConstTerm(f);
    ZZ comp = rep(ConstTerm(f));

    for (long i = 1; i <= deg(f); ++i) {
        if (comp < rep(coeff(f,i))) {
            f_i = coeff(f,i);
            comp = rep(coeff(f,i));
        }
        else {
            continue;
        }
    }

    return f_i;
}

ZZ evaluate(const ZZ_pX& f, const ZZ& x) {
    /*
        Returns the integer result when the polynomial f(y) (mod p)
        is evaluated for y = x in O(D) time where D is the number of
        terms in f(x) (mod p), including those with coefficients
        equal to zero.
    */

    long fDeg = deg(f);
    ZZ ans = rep(ConstTerm(f));

    auto xpow = ZZ(1);
    for (long i = 1; i <= fDeg; ++i) {
        xpow = xpow*x;
        ans += rep(coeff(f,i))*xpow;
    }

    return ans;
}

ZZ_pX ZZpXmultiply(const ZZ_pX& f, const ZZ_pX& g) {
    /*
        Returns the polynomial product s(x) =  f(x) * g(x) (mod p),
        using binary segmentation, in O(D) time where D is the number
        of terms in s(x), including those with coefficients equal to zero.
    */

    if (IsZero(f) || IsZero(g)) {                           // deg = -1 would give a product with no terms
        return ZZ_pX();
    }

    // ZZ termsF = to_ZZ(deg(f)+1);
    // ZZ termsG = to_ZZ(deg(g)+1);
    long termsF = deg(f)+1;
    long termsG = deg(g)+1;

    // REFACTOR 2.2: coefficients of a ZZ_pX lie in [0, p-1], so every
    // coefficient of s(x) is a sum of at most min(U,V) products, each at most
    // (p-1)^2. That bound is exact a priori and needs no scan for max(f_j) and
    // max(g_k) (the old bound also used max(U,V), which over-counts terms).
    ZZ minTermFG = to_ZZ(std::min(termsF,termsG));
    ZZ rhs = minTermFG * sqr(ZZ_p::modulus() - 1);          // rhs = min(U,V) * (p-1)^2 >= |s_i|

    // Every coefficient c_i of s(x) satisfies |c_i| <= rhs. For the signed
    // digit recovery below to be exact we need |c_i| < X/2, so (per JHD's
    // feedback on Algorithm 5.1.2) choose b with 2^(b-1) > rhs, leaving room
    // for the sign bit, and substitute X = 2^b, not 2^b - 1.
    long b = NumBits(rhs) + 1;                              // smallest b with 2^(b - 1) > rhs
    ZZ X = power2_ZZ(b);                                    // X = 2^b
    ZZ halfX = power2_ZZ(b - 1);                            // X/2 = 2^(b - 1)
    ZZ F = evaluate(f,X);                                   // evaluate f(x) for x = 2^b
    ZZ G = evaluate(g,X);                                   // evaluate g(x) for x = 2^b

    ZZ m = F * G;                                           // Integer multiply

    long fgTrm = termsF + termsG - 1;                       // Terms in polynomial product
    std::vector<ZZ_p> s(fgTrm);                             // Store coefficients and constant (std::vector: VLAs are not standard C++)

    // Extract the base-X digits of m: s_i = floor( m/X^i ) mod X, read as a
    // signed digit in (-X/2, X/2]. The comparison must be made on the integer
    // digit itself, before it is reduced mod p.
    ZZ digit;
    for (long i = 0; i < fgTrm; ++i) {
        if (i > 0) {
            m = m / X;                                      // N.B. -- NTL "/" operator floors result by default
        }
        digit = m % X;                                      // 0 <= digit < X

        if (digit > halfX) {                                // negative coefficient: borrow from the next digit
            digit -= X;
            m += X;
        }
        s[i] = to_ZZ_p(digit);                              // reduce mod p only after the sign has been recovered
    }

    ZZ_pX polyModProd;                                      // Base-b digits of m are desired coefficients
    polyModProd.SetLength(fgTrm);                           // set the length of the underlying coefficient vector to number of Terms in polynomial product

    for (long j = 0; j < fgTrm; ++j) {
        SetCoeff(polyModProd,j,s[j]);
    }

    polyModProd.normalize();                                // remove leading zeros on coefficients

    return polyModProd;
}

ZZ_pX ZZpPowMod(ZZ_pX a, ZZ n, const ZZ_pX& b) {
    /*
        Calculates a(x)^n (mod b(x), p) in O(log n) time.
    */

    ZZ_pX ans;                          // Initialise answer
    SetCoeff(ans,0,1);

    while (n > 0) {
        if (n % 2 == 1) {               // if (n is odd) then
            ans = ZZpXmultiply(ans,a);
            ans %= b;
        }
        a = ZZpXmultiply(a,a);          // a = a^2 (mod b)
        // a = sqr(a);                     // a = a^2 (mod b)
        a %= b;

        n /= 2;                         // n = n/2
    }

    return ans;
}
