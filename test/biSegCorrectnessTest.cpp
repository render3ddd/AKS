/*
    Correctness tests for the binary segmentation routines, checked against
    NTL's own (trusted) polynomial arithmetic.

    Build (from the repo root):
        clang++ -O2 -std=c++17 -pthread test/biSegCorrectnessTest.cpp -o test/biSegCorrectnessTest.out -lntl -lgmp -lm
    macOS / Homebrew: add  -I"$(brew --prefix)/include" -L"$(brew --prefix)/lib"

    Exit code is 0 only if every check passes.
*/

#include <cstdio>
#include <string>
#include <cstdlib>

#include "NTL/ZZ.h"
#include "NTL/ZZX.h"
#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
NTL_CLIENT

// The two headers each define overloads of getMaxCoeff/evaluate for different
// types (ZZX vs ZZ_pX), so they can be included together.
#include "../devVR/biSegMultiplyZZX.h"
#include "../devVR/biSegMultiplyZZpX.h"

static long failures = 0;
static long checks = 0;

static void check(bool ok, const std::string& what) {
    ++checks;
    if (!ok) {
        ++failures;
        if (failures <= 20) std::printf("FAIL: %s\n", what.c_str());
    }
}

static ZZX randomZZX(long terms, long bits, bool allowNegative) {
    ZZX f;
    for (long i = 0; i < terms; ++i) {
        ZZ c = RandomBits_ZZ(bits);
        if (allowNegative && RandomBnd(2) == 1) c = -c;
        SetCoeff(f, i, c);
    }
    return f;
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);   // keep output if a routine crashes
    SetSeed(ZZ(20211017));

    // ---------------------------------------------------------------
    // 1. polyMultiply (ZZX) -- signed integer coefficients
    // ---------------------------------------------------------------
    {
        // Hand-picked case from the dissertation (p.36): (x^2+2x+1)(x^2-x+1)
        ZZX f, g;
        SetCoeff(f, 2, 1); SetCoeff(f, 1, 2); SetCoeff(f, 0, 1);
        SetCoeff(g, 2, 1); SetCoeff(g, 1, -1); SetCoeff(g, 0, 1);
        check(polyMultiply(f, g) == f * g, "ZZX dissertation example (x^2+2x+1)(x^2-x+1)");

        // (x - 5)^2 : negative coefficient dominates
        ZZX h; SetCoeff(h, 1, 1); SetCoeff(h, 0, -5);
        check(polyMultiply(h, h) == h * h, "ZZX (x-5)^2");

        // all-negative polynomial
        ZZX k; SetCoeff(k, 0, -3); SetCoeff(k, 1, -7);
        check(polyMultiply(k, k) == k * k, "ZZX (-7x-3)^2");

        // constant * constant
        ZZX c3; SetCoeff(c3, 0, 3);
        check(polyMultiply(c3, c3) == c3 * c3, "ZZX 3*3");

        for (long t = 0; t < 300; ++t) {
            long tf = 1 + RandomBnd(40), tg = 1 + RandomBnd(40), bits = 1 + RandomBnd(80);
            ZZX a = randomZZX(tf, bits, true), b = randomZZX(tg, bits, true);
            check(polyMultiply(a, b) == a * b, "ZZX random signed #" + std::to_string(t));
        }
        for (long t = 0; t < 300; ++t) {
            long tf = 1 + RandomBnd(40), tg = 1 + RandomBnd(40), bits = 1 + RandomBnd(80);
            ZZX a = randomZZX(tf, bits, false), b = randomZZX(tg, bits, false);
            check(polyMultiply(a, b) == a * b, "ZZX random unsigned #" + std::to_string(t));
        }
    }

    // ---------------------------------------------------------------
    // 2. polyPowMod (ZZX) -- a(x)^n mod b(x) over Z  (Algorithm 5.1.3)
    //    Set AKS_SKIP_ZZX_POWMOD=1 to skip (the unfixed code can crash here).
    // ---------------------------------------------------------------
    if (!std::getenv("AKS_SKIP_ZZX_POWMOD")) {
        for (long t = 0; t < 40; ++t) {
            long r = 3 + RandomBnd(12);
            ZZX modulus; SetCoeff(modulus, r, 1); SetCoeff(modulus, 0, -1); // x^r - 1
            ZZX base; SetCoeff(base, 1, 1); SetCoeff(base, 0, -(1 + RandomBnd(5))); // x - a
            long e = 1 + RandomBnd(12);
            ZZX expect; SetCoeff(expect, 0, 1);
            for (long i = 0; i < e; ++i) expect = (expect * base) % modulus;  // NTL has no ZZX PowerMod
            ZZX got = polyPowMod(base, ZZ(e), modulus);
            check(got == expect, "ZZX polyPowMod (x-a)^" + std::to_string(e) + " mod x^" + std::to_string(r) + "-1");
        }
    }

    // ---------------------------------------------------------------
    // 3. ZZpXmultiply / ZZpPowMod -- the routines the Z version of AKS uses
    // ---------------------------------------------------------------
    {
        const char* moduli[] = {"7", "10007", "23456789", "1000000007",
                                "9007199254740997", "689960931088884849033689023336009222695077"};
        for (const char* ms : moduli) {
            ZZ p = conv<ZZ>(ms);
            ZZ_p::init(p);

            // constant * constant (U = V = 1): the case where X/2 < product <= bound
            ZZ_pX c; SetCoeff(c, 0, to_ZZ_p(p - 1));
            check(ZZpXmultiply(c, c) == c * c, std::string("ZZ_pX const*const mod ") + ms);

            for (long t = 0; t < 150; ++t) {
                long tf = 1 + RandomBnd(30), tg = 1 + RandomBnd(30);
                ZZ_pX a, b;
                random(a, tf); random(b, tg);
                // sometimes force small coefficients so that X can be smaller than 2p
                if (t % 3 == 0) { for (long i = 0; i <= deg(a); ++i) SetCoeff(a, i, to_ZZ_p(RandomBnd(4))); a.normalize(); }
                if (IsZero(a) || IsZero(b)) continue;      // zero polynomial is covered separately below
                check(ZZpXmultiply(a, b) == a * b, std::string("ZZ_pX random mod ") + ms + " #" + std::to_string(t));
            }

            // Coefficients of the product that lie in (X/2, p): e.g. (1 + x + ... + x^(k-1))^2
            // has middle coefficient k while the bound is max(U,V)*1*1 = k, so X < 2k.
            for (long k = 2; k <= 40; k += 7) {
                ZZ_pX ones;
                for (long i = 0; i < k; ++i) SetCoeff(ones, i, 1);
                check(ZZpXmultiply(ones, ones) == ones * ones,
                      std::string("ZZ_pX all-ones^2, k=") + std::to_string(k) + " mod " + ms);
            }

            long r = 7 + RandomBnd(40);
            ZZ_pX m = ZZ_pX(r, 1) - 1;
            for (long aa = 1; aa <= 4; ++aa) {
                ZZ_pX base = ZZ_pX(1, 1) - aa;
                check(ZZpPowMod(base, p, m) == PowerMod(base, p, m),
                      std::string("ZZpPowMod (x-") + std::to_string(aa) + ")^p mod x^r-1, p=" + ms);
            }
        }
    }

    // ---------------------------------------------------------------
    // 4. Zero polynomial (deg = -1). Set AKS_SKIP_ZERO=1 to skip
    //    (the unfixed code writes past a zero-length array here).
    // ---------------------------------------------------------------
    if (!std::getenv("AKS_SKIP_ZERO")) {
        ZZX z, one; SetCoeff(one, 0, 1);
        check(polyMultiply(z, one) == z, "ZZX 0*1");
        ZZ_p::init(ZZ(10007));
        ZZ_pX zp, onep; SetCoeff(onep, 0, 1);
        check(ZZpXmultiply(zp, onep) == zp, "ZZ_pX 0*1");
    }

    std::printf("%ld / %ld checks passed\n", checks - failures, checks);
    return failures == 0 ? 0 : 1;
}
