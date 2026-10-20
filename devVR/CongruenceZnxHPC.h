
#include <atomic>
#include <climits>
#include <thread>

#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
#include <NTL/BasicThreadPool.h>
NTL_CLIENT

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

    // Smallest failing a found so far (LONG_MAX = none). Shared between the
    // worker threads, so it is atomic; the old std::vector::push_back from
    // several threads at once was a data race (REFACTOR 1.6).
    std::atomic<long> smallestFail(LONG_MAX);

    long nthreads = AvailableThreads();
    if (nthreads > a) {
        nthreads = a;
    }

    printf("\nCommencing congruence tests:\n\n");
    if (a < 1) {
        return 0;
    }

    // Thread t tests a = t+1, t+1+T, t+1+2T, ... in increasing order and
    // stops as soon as its next a is not below the smallest failure already
    // found by any thread. Every a below the final minimum has therefore been
    // tested, so the result is the smallest failing a and does not depend on
    // scheduling, and threads stop early on composites.
    NTL_EXEC_INDEX(nthreads, t)

        context.restore();

        for (long j = t + 1; j <= a; j += nthreads) {
            if (j >= smallestFail.load(std::memory_order_relaxed)) {
                break;                              // a smaller failure is already known
            }
            printf("a = %ld\n", j);

            ZZ_pX c = ZZ_pX(1, 1) - j;              // c = x - a (mod n);
            ZZ_pX f = PowerMod(c, n, b);            // f = (x - a)^n (mod b, n) - LHS
            ZZ_pX g = d - j;                        // g = x^n - a (mod b, n) - RHS

            if (f != g) {                           // n is not prime: record min(smallestFail, j)
                long current = smallestFail.load(std::memory_order_relaxed);
                while (j < current && !smallestFail.compare_exchange_weak(current, j)) {
                }
                break;
            }
        }

    NTL_EXEC_INDEX_END

    long result = smallestFail.load();
    return result == LONG_MAX ? 0 : result;         // 0: n is prime
}