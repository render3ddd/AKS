
#include <vector>
#include <numeric>

#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
#include <NTL/BasicThreadPool.h>
NTL_CLIENT

// Per-a progress output (REFACTOR 2.6). Printing every a from the hot loop,
// and from every worker thread, costs I/O time and serialises the threads on
// stdout. It is off by default; build with -DAKS_VERBOSE=1 to restore it.
#ifndef AKS_VERBOSE
#define AKS_VERBOSE 0
#endif


unsigned int ncores = std::thread::hardware_concurrency(); // machine cores - may return 0 when not able to detect

long CongruenceZnx(const ZZ& n, const ZZ& r, const ZZ& r2, const long& a){
    // congruence test of polynomials in regular form
    SetNumThreads(ncores); // number of threads - should correspond to the number of available cores on your machine
    std::cout << AvailableThreads() << " of " << ncores << " threads available.\n";

    ZZ_p::init(n);                      // initialise mod n
    ZZ_pContext context;
    context.save();

    printf("\nCalculating x^r - 1 (mod n) ...\n");
    ZZ_pX b = ZZ_pX(to_long(r), 1) - 1; // b = x^r - 1 (mod n);
    printf("Done.\n");

    printf("Initialising x (mod n) ...\n");
    ZZ_pX e = ZZ_pX(1, 1);              // e = x (mod n)
    printf("Done.\n");

    printf("Calculating x^n (mod (x^r - 1), n) ...\n");
    ZZ_pX d = PowerMod(e, n, b);        // d = x^n (mod b, n)
    printf("Done.\n");

    vector <long> test;
    test.reserve(ncores);

    printf("\nCommencing congruence tests (1 <= a <= %ld):\n\n", a);
    long first, last;
    NTL_EXEC_RANGE(a,first,last)

        context.restore();

        for(long j = first; j < last; ++j){
            if (j < last) {                            // ensures each thread does not continue past j = a
                if (AKS_VERBOSE) {
                    printf("a = %ld\n",(j + 1));
                }

                ZZ_pX c = ZZ_pX(1, 1) - (j + 1);    // c = x - a (mod n);
                ZZ_pX f = PowerMod(c, n, b);        // f = (x - a)^n (mod b, n) - LHS
                ZZ_pX g = d - (j + 1);              // g = x^n - a (mod b, n) - RHS

                if(f != g){
                    test.push_back((j + 1)); // n is not prime
                    break;
                }
            }
        }

    NTL_EXEC_RANGE_END

    if (test.empty()) {
        return 0; // n is prime
    }
    else {
        long smallest = test.at(0);
        for (long k : test) {
            if (k < smallest) {
                smallest = k;
            }
        }

        return smallest;
    }
}