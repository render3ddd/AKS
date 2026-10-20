#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
// #include <NTL/BasicThreadPool.h>
NTL_CLIENT

// #include "biSegMultiplyZZX.h"
#include "biSegMultiplyZZpX.h"

unsigned int ncores = std::thread::hardware_concurrency(); // machine cores - may return 0 when not able to detect
// const auto SetNumThreads(ncores); // number of threads - should correspond to the number of available cores on your machine

long CongruenceZ(const ZZ& n, const ZZ& r, const ZZ& r2, const long& a) {
    // congruence test of polynomials in large integer form

    ZZ_p::init(n);                      // initialise mod n

    printf("\nCalculating x^r - 1 (mod n) ...\n");
    ZZ_pX b = ZZ_pX(to_long(r), 1) - 1; // b = x^r - 1 (mod n);
    printf("Done.\n");

    printf("Initialising x (mod n) ...\n");
    auto e = ZZ_pX(1, 1);              // e = x (mod n)
    printf("Done.\n");

    printf("Calculating x^n (mod (x^r - 1), n) ...\n");
    ZZ_pX d = ZZpPowMod(e, n, b);        // d = x^n (mod b, n)
    printf("Done.\n");

    printf("\nCommencing congruence tests:\n\n");
    for(long j = 1; j <= a; ++j){
        printf("a = %ld\n",j);

        ZZ_pX c = ZZ_pX(1, 1) - j ;         // c = x - a (mod n);
        ZZ_pX f = ZZpPowMod(c, n, b);        // f = (x - a)^n (mod b, n) - LHS
        ZZ_pX g = d - j;                    // g = x^n - a (mod b, n) - RHS

        if(f != g){
            return j; // n is not prime
        }
    }

    return 0; // n is prime
}