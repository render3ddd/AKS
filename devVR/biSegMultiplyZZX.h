#pragma once

#include "NTL/ZZX.h"

#include <math.h>
#include <vector>
#include <algorithm>

#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
NTL_CLIENT

inline ZZ powMod(ZZ a, ZZ n, const ZZ& b) {
    // calculates a^n (mod b) in O(log n)

    ZZ ans = ZZ(1);         // Initialise answer

    while (n > 0) {
        if (n % 2 == 1) {   // if (n is odd) then
            ans = (ans * a) % b;
        }
        a = (a * a) % b;    // a = a^2 (mod b)
        n /= 2;             // n = n/2
    }

    return ans;
}

inline ZZ getMaxCoeff(const ZZX& f) {
    /*
        Returns the largest absolute value of the integer coefficients,
        |f_i|, of the input polynomial, f(x), in O(D) time where D is the
        number of terms in f(x), including those with coefficients equal
        to zero. The binary segmentation size bound must use magnitudes:
        with signed coefficients, max{f_i} alone can be 1 (or negative)
        while |f_i| is large, e.g. f(x) = x - 5.
    */

    ZZ f_i = abs(ConstTerm(f));

    for (long i = 1; i <= deg(f); ++i) {
        if (f_i < abs(coeff(f,i))) {
            f_i = abs(coeff(f,i));
        }
    }

    return f_i;
}


// ZZ evaluate(const ZZX& f, const ZZ& x) {
//     /*
//         Returns the integer result when the polynomial f(y) is
//         evaluated for y = x in O(D) time where D is the number of
//         terms in f(x), including those with coefficients equal to
//         zero.
//     */

//     long fDeg = deg(f);
//     ZZ ans = ConstTerm(f);

//     for (long i = 1; i <= fDeg; ++i) {
//         ans += coeff(f,i)*power(x,i);
//     }

//     return ans;
// }

inline ZZ evaluate(const ZZX& f, const ZZ& x) {
    /*
        Returns the integer result when the polynomial f(y) (mod p)
        is evaluated for y = x in O(D) time where D is the number of
        terms in f(x) (mod p), including those with coefficients
        equal to zero.
    */

    long fDeg = deg(f);
    ZZ ans = ConstTerm(f);

    ZZ xpow = ZZ(1);
    for (long i = 1; i <= fDeg; ++i) {
        xpow = xpow*x;
        ans += coeff(f,i)*xpow;
    }

    return ans;
}

inline ZZX polyMultiply(const ZZX& f, const ZZX& g) {
    /*
        Returns the polynomial product s(x) = f(x) * g(x), using
        binary segmentation, in O(D) time where D is the number of
        terms in s(x), including those with coefficients equal to zero.
    */

    if (IsZero(f) || IsZero(g)) {                       // deg = -1 would give a product with no terms
        return ZZX();
    }

    long termsF = deg(f)+1;
    long termsG = deg(g)+1;

    ZZ maxTermFG = to_ZZ(std::max(termsF,termsG));
    ZZ maxCoeffF = getMaxCoeff(f);                      // max |f_j|
    ZZ maxCoeffG = getMaxCoeff(g);                      // max |g_k|

    ZZ rhs = maxTermFG * maxCoeffF * maxCoeffG;         // rhs = max(U,V) * max|f_j| * max|g_k| >= |s_i|

    // For the signed digit recovery below to be exact we need |s_i| < X/2,
    // so (per JHD's feedback on Algorithm 5.1.2) choose b with
    // 2^(b-1) > rhs, leaving room for the sign bit, and substitute
    // X = 2^b, not 2^b - 1.
    long b = 1;
    while (power2_ZZ(b - 1) <= rhs) {                   // Choose b such that 2^(b - 1) > max(U,V) * max|f_j| * max|g_k|
        b = b + 1;
    }

    ZZ X = power2_ZZ(b);                                // X = 2^b
    ZZ halfX = power2_ZZ(b - 1);                        // X/2 = 2^(b - 1)
    ZZ F = evaluate(f,X);                               // evaluate f(x) for x = 2^b
    ZZ G = evaluate(g,X);                               // evaluate g(x) for x = 2^b

    ZZ m = F * G;                                       // Integer multiply

    long fgTrm = termsF + termsG - 1;                   // Terms in polynomial product
    std::vector<ZZ> s(fgTrm);                           // Store coefficients and constant (std::vector: VLAs are not standard C++)

    // Extract the base-X digits of m: s_i = floor( m/X^i ) mod X, read as a
    // signed digit in (-X/2, X/2].
    for (long i = 0; i < fgTrm; ++i) {
        if (i > 0) {
            m = m / X;                                  // N.B. -- NTL "/" operator floors result by default
        }
        s[i] = m % X;                                   // 0 <= s_i < X

        if (s[i] > halfX) {                             // negative coefficient: borrow from the next digit
            s[i] -= X;
            m += X;
        }
    }

    ZZX polyProduct;                                    // Base-X digits of m are desired coefficients
    polyProduct.SetLength(fgTrm);                       // set the length of the underlying coefficient vector to number of Terms in polynomial product

    for (long j = 0; j < fgTrm; ++j) {
        SetCoeff(polyProduct,j,s[j]);
    }

    polyProduct.normalize();                            // remove leading zeros on coefficients

    return polyProduct;
}

inline ZZX polyMulMod(const ZZX& f, const ZZX& g, const ZZX& h) {
    /*
        Returns the polynomial product s(x) = f(x) * g(x) (mod h(x)),
        using binary segmentation, in O(D) time where D is the number
        of terms in s(x), including those with coefficients equal to zero.
    */

    return polyMultiply(f,g) % h;
}

inline ZZX polyMulModN(const ZZX& f, const ZZX& g, const ZZ& n) {
    /*
        Returns the polynomial product s(x) = f(x) * g(x) (mod n), using
        binary segmentation, in O(D) time where D is the number of terms
        in s(x), including those with coefficients equal to zero.
    */

    ZZX polyProdModN = polyMultiply(f,g);

    for (long j = 0; j <= deg(polyProdModN); ++j) {
        SetCoeff(polyProdModN,j,(coeff(polyProdModN,j) % n));     // s[j] (mod n)
    }

    polyProdModN.normalize();                           // remove leading zeros on coefficients

    return polyProdModN;
}

inline ZZX polyPowMod(ZZX a, ZZ n, const ZZX& b) {
    /*
        Calculates a(x)^n (mod b(x)) in O((log n)*O(D)) time.
    */

    ZZX ans;                    // Initialise answer
    SetCoeff(ans,0,1);
    a %= b;

    while (n > 0) {
        if (n % 2 == 1) {       // if (n is odd) then
            ans = polyMultiply(ans,a) % b;
        }
        a = polyMultiply(a,a) % b;  // a = a^2 (mod b) -- line 6 of Algorithm 5.1.3 was missing "mod b"
        n /= 2;                 // n = n/2
    }

    return ans;
}