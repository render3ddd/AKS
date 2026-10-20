
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

    // ZZ termsF = to_ZZ(deg(f)+1);
    // ZZ termsG = to_ZZ(deg(g)+1);
    long termsF = deg(f)+1;
    long termsG = deg(g)+1;

    ZZ maxTermFG = to_ZZ(std::max(termsF,termsG));
    ZZ_p maxCoeffF = getMaxCoeff(f);
    ZZ_p maxCoeffG = getMaxCoeff(g);
    // ZZ maxCoeff = ZZ_p::modulus();

    ZZ rhs = maxTermFG * rep(maxCoeffF) * rep(maxCoeffG);   // rhs = max(U,V) * max(f_j) * max(g_k)
    // ZZ rhs = maxTermFG * sqr(maxCoeff);                     // rhs = max(U,V) * max(f_j) * max(g_k)

    long b = 1;
    ZZ lhs = (power2_ZZ(b - 1));                            // lhs = 2^(b - 1)
    while (lhs <= rhs) {                                    // Choose b such that 2^(b - 1) > max(U,V) * max(f_j) * max(g_k)
        b = b + 1;
        lhs = (power2_ZZ(b - 1));
    }

    ZZ X = lhs;
    ZZ F = evaluate(f,X);                                   // evaluate f(x) for x = 2^(b - 1)
    ZZ G = evaluate(g,X);                                   // evaluate g(x) for x = 2^(b - 1)

    ZZ m = F * G;                                           // Integer multiply

    // long fgDeg = deg(f) + deg(g);                           // Degree of polynomial product
    // long fgTrm = deg(f) + deg(g) + 1;                       // Terms in polynomial product
    long fgTrm = termsF + termsG - 1;                       // Terms in polynomial product
    ZZ_p s[fgTrm];                                          // Store coefficients and constant in an array

    s[0] = to_ZZ_p(m % lhs);                                // Reassemble coefficients into signal
    if (rep(s[0]) > lhs/2) {                                // Extract next b bits: s_i = floor( m/(2^(b - 1))^i ) mod 2^(b - 1)
        s[0] = s[0] - to_ZZ_p(lhs);
        m = m + lhs;
    }
    for (long i = 1; i < fgTrm; ++i) {
        m = m / lhs;                                        // N.B. -- NTL "/" operator floors result by default
        s[i] = to_ZZ_p(m % lhs);

        if (rep(s[i]) > lhs/2) {
            s[i] = s[i] - to_ZZ_p(lhs);
            m = m + lhs;
        }
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
