
#include <math.h>
#include <vector>
#include <algorithm>

#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
NTL_CLIENT

ZZ powMod(ZZ a, ZZ n, const ZZ& b) {
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

ZZ getMaxCoeff(const ZZX& f) {
    /*
        Returns the largest integer coefficient, f_i, of the input
        polynomial, f(x), in O(D) time where D is the number of terms
        in f(x), including those with coefficients equal to zero.
    */

    ZZ f_i = ConstTerm(f);

    for (long i = 1; i <= deg(f); ++i) {
        if (f_i < coeff(f,i)) {
            f_i = coeff(f,i);
        }
        else {
            continue;
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

ZZ evaluate(const ZZX& f, const ZZ& x) {
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

ZZX polyMultiply(const ZZX& f, const ZZX& g) {
    /*
        Returns the polynomial product s(x) = f(x) * g(x), using
        binary segmentation, in O(D) time where D is the number of
        terms in s(x), including those with coefficients equal to zero.
    */

    ZZ termsF = to_ZZ(deg(f)+1);
    ZZ termsG = to_ZZ(deg(g)+1);

    ZZ maxTermFG = to_ZZ(std::max(termsF,termsG));
    ZZ maxCoeffF = getMaxCoeff(f);
    ZZ maxCoeffG = getMaxCoeff(g);

    ZZ rhs = maxTermFG * maxCoeffF * maxCoeffG;         // rhs = max(U,V) * max(f_j) * max(g_k)

    long b = 1;
    ZZ lhs = (power2_ZZ(b - 1));                        // lhs = 2^(b - 1)
    while (lhs <= rhs) {                                // Choose b such that 2^(b - 1) > max(U,V) * max(f_j) * max(g_k)
        b = b + 1;
        lhs = (power2_ZZ(b - 1));
    }

    ZZ X = lhs;
    ZZ F = evaluate(f,X);                               // evaluate f(x) for x = 2^b - 1
    ZZ G = evaluate(g,X);                               // evaluate g(x) for x = 2^b - 1

    ZZ m = F * G;                                       // Integer multiply

    // long fgDeg = deg(f) + deg(g);                       // Degree of polynomial product
    long fgTrm = deg(f) + deg(g) + 1;                   // Terms in polynomial product
    ZZ s[fgTrm];                                        // Store coefficients and constant in an array

    s[0] = m % lhs;                                     // Reassemble coefficients into signal
    if (s[0] > lhs/2) {                                 // Extract next b bits: s_i = floor( m/(2b−1)^i ) mod 2^(b - 1)
        s[0] = s[0] - lhs;
        m = m + lhs;
    }
    for (long i = 1; i < fgTrm; ++i) {
        m = m / lhs;                                    // N.B. -- NTL "/" operator floors result by default
        s[i] = m % lhs;

        if (s[i] > lhs/2) {
            s[i] = s[i] - lhs;
            m = m + lhs;
        }
    }

    ZZX polyProduct;                                    // Base-b digits of m are desired coefficients
    polyProduct.SetLength(fgTrm);                       // set the length of the underlying coefficient vector to number of Terms in polynomial product

    for (long j = 0; j < fgTrm; ++j) {
        SetCoeff(polyProduct,j,s[j]);
    }

    polyProduct.normalize();                            // remove leading zeros on coefficients

    return polyProduct;
}

ZZX polyMulMod(const ZZX& f, const ZZX& g, const ZZX& h) {
    /*
        Returns the polynomial product s(x) = f(x) * g(x) (mod h(x)),
        using binary segmentation, in O(D) time where D is the number
        of terms in s(x), including those with coefficients equal to zero.
    */

    ZZ termsF = to_ZZ(deg(f)+1);
    ZZ termsG = to_ZZ(deg(g)+1);

    ZZ maxTermFG = to_ZZ(std::max(termsF,termsG));
    ZZ maxCoeffF = getMaxCoeff(f);
    ZZ maxCoeffG = getMaxCoeff(g);

    ZZ rhs = maxTermFG * maxCoeffF * maxCoeffG;         // rhs = max(U,V) * max(f_j) * max(g_k)

    long b = 1;
    ZZ lhs = (power2_ZZ(b)) - 1;                        // lhs = 2^(b - 1)
    while (lhs <= rhs) {                                 // Choose b such that 2^(b - 1) > max(U,V) * max(x_i) * max(y_k)
        b = b + 1;
        lhs = (power2_ZZ(b)) - 1;
    }

    ZZ X = lhs;
    ZZ F = evaluate(f,X);                               // evaluate f(x) for x = 2^(b - 1)
    ZZ G = evaluate(g,X);                               // evaluate g(x) for x = 2^(b - 1)

    ZZ m = F * G;                                       // Integer multiply

    // long fgDeg = deg(f) + deg(g);                       // Degree of polynomial product
    long fgTrm = deg(f) + deg(g) + 1;                   // Terms in polynomial product
    ZZ s[fgTrm];                                        // Store coefficients and constant in an array

    s[0] = m % lhs;                                     // Reassemble coefficients into signal
    if (s[0] > lhs/2) {                                 // Extract next b bits: s_i = floor( m/(2b−1)^i ) mod 2^(b - 1)
        s[0] = s[0] - lhs;
        m = m + lhs;
    }
    for (long i = 1; i < fgTrm; ++i) {
        m = m / lhs;                                    // N.B. -- NTL "/" operator floors result by default
        s[i] = m % lhs;

        if (s[i] > lhs/2) {
            s[i] = s[i] - lhs;
            m = m + lhs;
        }
    }

    ZZX polyProduct;                                    // Base-b digits of m are desired coefficients
    polyProduct.SetLength(fgTrm);                       // set the length of the underlying coefficient vector to number of Terms in polynomial product

    for (long j = 0; j < fgTrm; ++j) {
        SetCoeff(polyProduct,j,s[j]);
    }

    polyProduct.normalize();                            // remove leading zeros on coefficients

    ZZX polyProdMod = polyProduct % h;

    return polyProdMod;
}

ZZX polyMulModN(const ZZX& f, const ZZX& g, const ZZ& n) {
    /*
        Returns the polynomial product s(x) = f(x) * g(x) (mod n), using
        binary segmentation, in O(D) time where D is the number of terms
        in s(x), including those with coefficients equal to zero.
    */

    ZZ termsF = to_ZZ(deg(f)+1);
    ZZ termsG = to_ZZ(deg(g)+1);

    ZZ maxTermFG = to_ZZ(std::max(termsF,termsG));
    ZZ maxCoeffF = getMaxCoeff(f);
    ZZ maxCoeffG = getMaxCoeff(g);

    ZZ rhs = maxTermFG * maxCoeffF * maxCoeffG;     // rhs = max(U,V) * max(f_j) * max(g_k)
    long b = 1;
    ZZ lhs = (power2_ZZ(b)) - 1;                    // lhs = 2^(b - 1)
    while (lhs <= rhs) {                            // Choose b such that 2^(b - 1) > max(U,V) * max(x_i) * max(y_k)
        b = b + 1;
        lhs = (power2_ZZ(b)) - 1;
    }

    ZZ X = lhs;
    ZZ F = evaluate(f,X);                           // evaluate f(x) for x = 2^(b - 1)
    ZZ G = evaluate(g,X);                           // evaluate g(x) for x = 2^(b - 1)

    ZZ m = F * G;                                   // Integer multiply

    // long fgDeg = deg(f) + deg(g);                   // Degree of polynomial product
    long fgTrm = deg(f) + deg(g) + 1;               // Terms in polynomial product
    ZZ s[fgTrm];                                    // Store coefficients and constant in an array

    s[0] = m % lhs;                                 // Reassemble coefficients into signal
    if (s[0] > lhs/2) {                             // Extract next b bits: s_i = floor( m/(2b−1)^i ) mod 2^(b - 1)
        s[0] = s[0] - lhs;
        m = m + lhs;
    }
    for (long i = 1; i < fgTrm; ++i) {
        m = m / lhs;                                // N.B. -- NTL "/" operator floors result by default
        s[i] = m % lhs;

        if (s[i] > lhs/2) {
            s[i] = s[i] - lhs;
            m = m + lhs;
        }
    }

    ZZX polyProdModN;                               // Base-b digits of m are desired coefficients
    polyProdModN.SetLength(fgTrm);                  // set the length of the underlying coefficient vector to number of Terms in polynomial product

    for (long j = 0; j < fgTrm; ++j) {
        SetCoeff(polyProdModN,j,(s[j] % n));        // s[j] (mod n)
    }

    polyProdModN.normalize();                       // remove leading zeros on coefficients

    return polyProdModN;
}

ZZX polyPowMod(ZZX a, ZZ n, ZZX b) {
    /*
        Calculates a(x)^n (mod b(x)) in O((log n)*O(D)) time.
    */

    ZZX ans;                    // Initialise answer
    SetCoeff(ans,0,1);

    while (n > 0) {
        if (n % 2 == 1) {       // if (n is odd) then
            ans = polyMultiply(ans,a) % b;
        }
        a = polyMultiply(a,a);  // a = a^2 (mod b)
        n /= 2;                 // n = n/2
    }

    return ans;
}
