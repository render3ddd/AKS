
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
        Returns the polynomial product s(x) = f(x) * g(x) (mod p) using
        binary segmentation (Kronecker substitution) in O(M(D*b)) time, where
        D is the number of terms in s(x) and b the bits per coefficient slot.

        REFACTOR 2.1: the previous version evaluated f(2^b) and g(2^b) with
        D multiplications of ever-larger integers, and extracted the digits
        with D full-width "% X" and "/ X" operations, so it was O(D^2 * b).
        Because X = 2^b, substitution and extraction are just bit packing:
        each coefficient occupies its own fixed-width slot. Coefficients of
        a ZZ_pX lie in [0, p), so every product coefficient is
        non-negative and fits in an unsigned slot (no sign bit needed).
    */

    if (IsZero(f) || IsZero(g)) {                           // deg = -1 would give a product with no terms
        return ZZ_pX();
    }

    long termsF = deg(f)+1;
    long termsG = deg(g)+1;

    ZZ maxTermFG = to_ZZ(std::max(termsF,termsG));
    ZZ_p maxCoeffF = getMaxCoeff(f);
    ZZ_p maxCoeffG = getMaxCoeff(g);

    ZZ rhs = maxTermFG * rep(maxCoeffF) * rep(maxCoeffG);   // rhs = max(U,V) * max(f_j) * max(g_k) >= every s_i >= 0

    long slotBytes = (NumBits(rhs) + 7) / 8;                // byte-aligned slot of b = 8*slotBytes bits, 2^b > rhs
    if (slotBytes < 1) {
        slotBytes = 1;
    }

    // F = f(2^b), G = g(2^b): place coefficient i in bytes [i*slotBytes, (i+1)*slotBytes), little-endian
    auto pack = [slotBytes](const ZZ_pX& h, long terms) {
        std::vector<unsigned char> buf(terms * slotBytes, 0);
        for (long i = 0; i <= deg(h); ++i) {
            BytesFromZZ(&buf[i * slotBytes], rep(h.rep[i]), slotBytes);
        }
        ZZ z;
        ZZFromBytes(z, buf.data(), buf.size());
        return z;
    };

    ZZ m = pack(f, termsF) * pack(g, termsG);               // one big-integer multiply

    long fgTrm = termsF + termsG - 1;                       // Terms in polynomial product
    std::vector<unsigned char> out(fgTrm * slotBytes, 0);
    BytesFromZZ(out.data(), m, out.size());                 // m < 2^(b*fgTrm), so nothing is lost

    ZZ_pX polyModProd;                                      // slot i of m is coefficient s_i
    polyModProd.SetLength(fgTrm);
    ZZ digit;
    for (long i = 0; i < fgTrm; ++i) {
        ZZFromBytes(digit, &out[i * slotBytes], slotBytes);
        conv(polyModProd.rep[i], digit);                    // reduce mod p
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
