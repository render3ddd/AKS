# Known issues, inefficiencies and carried-over misconceptions

This review covers the C++ code in this repository (`devVR/` primarily) and the dissertation *Implementing alternative polynomial multiplication algorithms and the AKS primality test in parallel* (Renders, 2021). It records:

1. code errors,
2. inefficiencies,
3. assertions in the dissertation that are incorrect or unsupported *and* show up in the code,

apart from the binary segmentation bug the supervisor flagged, which is fixed on the local branch `fix/binary-segmentation-sign-bit` (§1.1).

Every claim marked **[verified]** was reproduced by compiling and running the code (Linux, GCC 13 / Clang, NTL 11.5.1, GMP 6.3) or checked against an independent reference (SymPy). Line numbers refer to `main` unless stated otherwise.

**Severity:** 🔴 wrong results or undefined behaviour · 🟠 weakens the correctness argument or the experimental conclusions · 🟡 efficiency or maintainability.

---

## Summary

| # | Issue | Where | Severity | Status |
|---|---|---|---|---|
| 1.1 | Binary segmentation bound / evaluation point / sign handling; `polyPowMod` missing `mod b` (supervisor, pp. 37–38) | `biSegMultiplyZZX.h`, `biSegMultiplyZZpX.h` | 🔴 | **Fixed on branch** |
| 1.2 | `r` does not satisfy ord_r(n) > ⌊log₂ n⌋²; the search uses ln, truncates, tests the wrong condition and overflows | `LenstraZnx.cpp`, `LenstraZ.cpp` | 🟠 | Open |
| 1.3 | Range of *a* (φ(r)−1 or λ(r)−1) is not the theorem's bound; the Carmichael speed-up is an artefact of the invalid `r` | as above, `Carmichael.h` | 🟠 | Open |
| 1.4 | AKS step "if n ≤ r then PRIME" missing; inverted in Algorithm 2.4.1 | as above | 🟡 | Open |
| 1.5 | `PerfectPower` misses every perfect power with base ≥ 2⁶³ and can loop forever | `PerfectPower.h` | 🔴 | Open |
| 1.6 | Data race on `std::vector::push_back` inside the thread pool | `CongruenceZnxHPC.h:51` | 🔴 | Open |
| 1.7 | Large *n* printed and logged via `to_long(n)` (truncated) | `Lenstra*.cpp` | 🟡 | Open |
| 1.8 | `Carmichael.h`: correct values, but O(r) trial division, fragile LCM indexing, uninitialised variable | `Carmichael.h` | 🟡 | Open |
| 1.9 | Headers define globals and functions, have no include guards and depend on include order | all headers | 🟡 | Open |
| 1.10 | Stale or non-portable test programs | `test/` | 🟡 | Open |
| 1.11 | Review notes on the 2024 WIP commit | branch commit `2e78c69` | 🟡 | Open |
| 2.1 | Binary segmentation is O(D²·b), not O(D); a packed version runs within ~3× of NTL | `biSegMultiplyZZpX.h` | 🟠 | Open |
| 2.2–2.8 | Other inefficiencies (bounds, reductions, `r` search, logging, partitioning) | various | 🟡 | Open |
| 3 | Dissertation assertions that carry into the code | see table | — | — |

---

## 1. Code errors

### 1.1 🔴 Binary segmentation and `polyPowMod`: fixed on `fix/binary-segmentation-sign-bit`

This is the bug in the supervisor's marker notes ("p. 37: we want 2^{b−1} > …, not 2^b − 1 (to allow for the sign bit), and you need to substitute in 2^b"; "p. 38: Line 6 also needs to be modulo b"). The branch fixes all of the following:

- **Bound and evaluation point.** Choose *b* with 2^(b−1) > max(U,V)·max|f_j|·max|g_k| and evaluate at X = 2^b. Every product coefficient then satisfies |s_i| < X/2, and the signed-digit recovery is exact. The original code used X = 2^b − 1 > bound. Commit `e173feb` used X = 2^(b−1) > bound, which still has no sign margin, and it changed only two of the four copies of the loop.
- **Comparing after reduction** (`ZZpXmultiply`). The code compared `rep(s_i)`, the digit *already reduced mod p*, against X/2. Any true coefficient in (X/2, p) was then corrupted. Example: (1 + x + … + x⁸)² mod 10007 is wrong **[verified]**.
- **Signed coefficients** (`getMaxCoeff(ZZX)`). The function returned the largest *signed* coefficient, so for x − 5 the "max" was 1, and (x − 5)², (−7x − 3)² and 3·3 were all wrong **[verified]**. It must use |f_j|. This is also an error in Algorithm 5.1.2 itself (§3).
- **`polyPowMod` (`ZZX`).** Line 6 of Algorithm 5.1.3 omitted `mod b`, so a(x) grew without bound. The wrong products eventually produced a zero polynomial, and the code then wrote past a zero-length array and **segfaulted** **[verified with ASan]**.
- **Zero operand.** `deg = −1` gave a zero-length VLA and an out-of-bounds write. Fixed with an early return.
- **Duplication and VLAs.** `polyMulMod` and `polyMulModN` now reuse `polyMultiply` instead of three diverging copies. Non-standard VLAs are replaced with `std::vector`.

`test/biSegCorrectnessTest.cpp` (on the branch) checks every routine against NTL's arithmetic: the old code failed 35 of 1,608 checks and crashed; the fixed code passes all 1,608, clean under ASan and UBSan. The Z version of AKS on the branch agrees with a reference primality test for every n from 2 to 699.

In practice the `ZZ_pX` bug rarely fired inside AKS, because coefficients mod a large *n* make X ≫ p. That is why the dissertation's tests passed, but the routine was not correct.

### 1.2 🟠 The chosen `r` does not satisfy the AKS/Lenstra condition

Every AKS-class proof needs ord_r(n) > ⌊log₂ n⌋². The search loop (`LenstraZnx.cpp` ≈ lines 160–180, identical in `LenstraZ.cpp`) has four separate problems:

1. `log(n)` in NTL is the **natural** log, but the dissertation states that log means log₂ (p. 7). ln n is about 0.69 log₂ n, so the threshold is about 2× too small.
2. `to_long(log(n))` truncates *before* squaring.
3. The test is inverted. It accepts the first r for which n^v ≡ 1 (mod r) for **some** v in [⌊ln n⌋², r]. That only shows the order *divides* some v ≥ the threshold, not that the order *exceeds* it.
4. `power_long(to_long(n), to_long(v))` computes n^v in a 64-bit `long`. It overflows for every realistic input, which is undefined behaviour, and `to_long(n)` truncates n > 2⁶³. The result of the search is effectively arbitrary.

**[verified]** Values chosen by the code compared with the smallest valid r (SymPy `n_order`):

| n | code's r | ord_r(n) | required ord > | smallest valid r |
|---|---|---|---|---|
| 10007 | 87 | 28 | 169 | 173 |
| 23456789 | 256 | 64 | 576 | 599 |
| 689960931088884849033689023336009222695077 (42 digits) | 9323 | 4661 | 19044 | 19051 |

The theorem therefore does not guarantee the program's PRIME answers. No false positive turned up in 409 composites tested (semiprimes, Carmichael numbers, prime powers), and the program agrees with a reference for n ≤ 1499, but that is empirical luck, not a proof. The benchmark timings in Chapters 5–6 were taken with r about half its valid size. With a valid r (≈2×) and a valid range of a (≈2×), the 42-digit run would take roughly 4× longer.

**Fix.** Compute ⌊log₂ n⌋ as `NumBits(n) − 1`. For each r coprime to n, compute ord_r(n) by repeated multiplication of `n % r` modulo r (a word-sized loop), stopping as soon as k exceeds the threshold.

### 1.3 🟠 Range of *a*, and the Carmichael λ(r) optimisation

The supervisor noted on p. 16 that the Lemma says "all such a" up to ⌈√r·log n⌉, but the algorithm tests a different range. The code tests 1 ≤ a ≤ φ(r) − 1 (Z version) or 1 ≤ a ≤ λ(r) − 1 (Zn[x] version). The standard bound is a ≤ ⌊√φ(r)·log₂ n⌋, capped below n.

- No theorem supports λ(r) − 1. It can fall far below the required bound. For example, λ(256) − 1 = 63, while √φ(256)·log₂ n ≈ 277 for n = 23456789.
- The large λ-versus-φ reductions in §5.2.4 and Appendix A.2 occur because the invalid r values are often **composite** (256 = 2⁸, 87 = 3·29). The smallest *valid* r is usually prime: 173, 599, 2441 and 19051 in the cases above. For prime r, λ(r) = φ(r), so the optimisation gives no gain.

**Fix.** Use the Lemma's bound. Keep λ/φ only as switchable experiment modes, clearly labelled as not certified.

### 1.4 🟡 Missing "n ≤ r ⇒ PRIME" step

The code relies on `while (r < n)` and never handles r ≥ n. Algorithm 2.4.1 in the dissertation has this step inverted: lines 8–9 say "if n ≤ r then return NOT PRIME". AKS says output PRIME, because every a < n has already passed the gcd test. If the search loop ever finishes without finding r, `r1` is left at 0 and is then passed to `Carmichael(0)`.

### 1.5 🔴 `PerfectPower` is wrong for large inputs, and can hang

`PerfectPower.h` runs Newton's method in `double` and assigns `a = to_long(e)`. **[verified]**

- It correctly detects aᵇ for bases up to 4294967311 (2³² + 15).
- It reports **"not a perfect power"** for every b ∈ {2,3,4,5} when a = 18446744073709551629 (> 2⁶⁴), and the same for a = 10²⁰ + 39 with b ∈ {2,3,4}.
- For (10²⁰ + 39)⁵ it **does not terminate**.

**Fix.** Use an exact integer k-th root, for example Newton iteration entirely in `ZZ`, or GMP's `mpz_perfect_power_p`. Only exponents b ≤ ⌊log₂ n⌋ need testing, and only prime b.

### 1.6 🔴 Data race in the parallel congruence test

In `CongruenceZnxHPC.h:51`, `test.push_back(j + 1)` runs concurrently from the NTL thread pool's workers on a shared `std::vector`, with no lock. That is undefined behaviour. **[verified]** ThreadSanitizer reports the race at that exact line for n = 23456611 and n = 100160063.

The loop has related problems:

- `printf` runs inside the parallel loop.
- `if (j < last)` is redundant.
- `NTL_EXEC_RANGE` gives each thread a fixed contiguous block of *a* values with no early exit, so on a composite every thread finishes its first failing *a* before the call returns (as noted in §6.2.1).

**Fix.** Use a per-thread result slot (`std::vector<long> firstFail(nthreads)`) or `std::atomic<long>` with a compare-and-swap minimum, plus a shared `std::atomic<bool>` that every iteration checks. Interleaving the *a* values (thread t handles a ≡ t mod T) finds the smallest failing *a* sooner than contiguous blocks do.

### 1.7 🟡 Output truncates large *n*

`std::printf("%ld …", to_long(n))` and `std::to_string(to_long(n))` truncate any n > 2⁶³. The console output, and the CSV `note` field for perfect powers and divisors, are therefore wrong for the 42-digit test primes. (The first CSV column is fine because it streams the `ZZ` with `<<`.) Print with `std::cout << n`.

Separately, `printf("%ld", time)` with a chrono `count()` is a `-Wformat` warning on macOS, where the type is `long long`. The WIP changes it to `%lld`, which then warns on Linux, where the type is `long`. Cast to `long long` and use `%lld`.

### 1.8 🟡 `Carmichael.h`

**[verified]** `Carmichael(r)` and `Euler(r)` return the correct value for every r from 2 to 3,000, with no out-of-bounds access under ASan and `_GLIBCXX_ASSERTIONS`. Even so:

- `isPrime` and `primeFactors` trial-divide up to **n/2** (see p. 4 in §3), and `primeFactors` pushes *every* prime ≤ n/2 with exponent 0. Both are O(n) work for each call.
- `LCM()` aligns factors by **position** in each factor list, not by prime. That only works because every λ(pᵉ) value except 1 and 2 is even, which fixes index 0 to the prime 2. It is easy to break.
- In `subLambda`, `lambda` is uninitialised on the (unreachable) `else` branch, which produces a compiler warning. `isPrime(1)` returns `true`.

**Fix.** Factor r once by trial division up to √r, and derive both φ(r) and λ(r) from that single factorisation.

### 1.9 🟡 Header hygiene

- No header has an include guard or `#pragma once`.
- Headers *define* functions and the global `unsigned int ncores`. As a result `CongruenceZnx.h` and `CongruenceZnxHPC.h` cannot both be included. Listing C.1 in the dissertation includes both and calls `CongruenceZnxHPC()`, so the printed listing does not compile and does not match the repository, where the function is called `CongruenceZnx`.
- Headers rely on whatever the `.cpp` included first. For example, `CongruenceZnxHPC.h` uses `std::thread`, `std::cout` and `printf` without including `<thread>`, `<iostream>` or `<cstdio>`. On the WIP branch, `<chrono>` is no longer included at all and only arrives through `<thread>` (which happens to work with both libstdc++ and libc++).
- `using namespace NTL` (via `NTL_CLIENT`) appears in every header.

**Fix.** Make headers self-contained, mark functions `inline` (or move them into `.cpp` files), and pass `ncores` as a parameter.

### 1.10 🟡 Tests

- `test/chronoTest.cpp` does not compile on any platform. It calls `getMaxCoeff`, `evaluate` and `polyMultiply` with `ZZX` signatures through the `ZZ_pX` header.
- `test/getHardwareTest.cpp` uses x86-only `<cpuid.h>`.
- `test/*_cli.h` and `perflogBiSegMultiplyZZpX.h` are older copies of the binary segmentation code that still contain the unfixed algorithm.
- Nothing runs automatically.

**Fix.** Delete the stale copies, keep `biSegCorrectnessTest.cpp`, add end-to-end AKS checks against NTL's `ProbPrime` for n < 10⁶, and wire them into a `Makefile` or CMake `ctest` target.

### 1.11 🟡 Review of the 2024 WIP commit (`2e78c69` on the branch)

The WIP was committed as-is and is not assumed correct. Observations:

- **Builds.** It compiles and runs with GCC, Clang/libstdc++ and Clang/libc++. The Zn[x] program matches a reference for every n from 2 to 1,499.
- **Include changes.** `<chrono>` was removed but `std::chrono` is still used (see §1.9). Moving the NTL includes after the project headers is only safe because `PerfectPower.h` now includes NTL itself and `Euler.h` qualifies `NTL::`.
- **`main()` reads n from stdin.** That is good for testing, but `nos[]`, `nosSize`, `nosEnd` and `prime` are now unused.
- **Deleted `.cu` files.** `.vscode/launch.json` still points at `devVR/LenstraZnx.cu`. The deleted kernels could never have worked, since NTL types cannot be used in `__global__` code (Ch. 3.4).
- **2^(b−1) comment edits.** These describe the `e173feb` change, which is superseded by the fix in §1.1.

---

## 2. Inefficiencies

### 2.1 🟠 Binary segmentation is quadratic as implemented

Section 5.1.2 states the routine is O(D_s). It is O(D²·b) as written:

- `evaluate()` builds `xpow = xpow*x` and `coeff*xpow`, which multiplies ever-larger integers, roughly D times.
- The extraction loop performs D full-width `m % X` and `m / X` operations on an integer of D·b bits.

That is why "construction" took about 90% of the run time (§5.2.2). **[verified]** Squaring a random `ZZ_pX` mod the 42-digit prime:

| degree D | `ZZpXmultiply` (fixed) | NTL `sqr` | ratio |
|---|---|---|---|
| 500 | 84 ms | 0.71 ms | 119× |
| 2000 | 975 ms | 0.94 ms | 1032× |
| 8000 | 15,963 ms | 4.99 ms | 3197× |

Time ×4 per doubling of D: quadratic.

With X = 2^b, substitution and extraction are just **bit packing**: Kronecker substitution, the O(M(D·b)) method Crandall & Pomerance describe. A 25-line prototype using byte-aligned slots, `ZZFromBytes` and `BytesFromZZ` gives the following timings **[verified, product identical to NTL]**:

| degree D | Kronecker prototype | NTL `sqr` | ratio |
|---|---|---|---|
| 500 | 0.55 ms | 0.43 ms | 1.3× |
| 2000 | 3.33 ms | 1.02 ms | 3.3× |
| 16000 | 37.0 ms | 12.8 ms | 2.9× |

That is the 2–4× overhead §4.3.2 originally anticipated. The conclusions that the Z version is "several orders of magnitude" slower, that the overhead is "a factor of 14", and that it is "redundant" (§5.2.2–5.2.3, §7.3–7.4) reflect the implementation, not the method.

```cpp
// Kronecker substitution for ZZ_pX: coefficients are in [0, p-1], so no sign bit is needed.
ZZ_pX kronMul(const ZZ_pX& f, const ZZ_pX& g) {
    if (IsZero(f) || IsZero(g)) return ZZ_pX();
    long U = deg(f) + 1, V = deg(g) + 1, T = U + V - 1;
    ZZ bound = ZZ(std::min(U, V)) * sqr(ZZ_p::modulus() - 1);   // exact a-priori bound
    long w = (NumBits(bound) + 7) / 8;                             // slot width in bytes
    auto pack = [&](const ZZ_pX& h, long n) {
        std::vector<unsigned char> buf(n * w, 0);
        for (long i = 0; i <= deg(h); ++i) BytesFromZZ(&buf[i * w], rep(h.rep[i]), w);
        ZZ z; ZZFromBytes(z, buf.data(), buf.size()); return z;
    };
    ZZ m = pack(f, U) * pack(g, V);                                // one big GMP multiply
    std::vector<unsigned char> out(T * w, 0);
    BytesFromZZ(out.data(), m, out.size());
    ZZ_pX s; s.SetLength(T); ZZ c;
    for (long i = 0; i < T; ++i) { ZZFromBytes(c, &out[i * w], w); conv(s.rep[i], c); }
    s.normalize(); return s;
}
```

This also disproves the §5.2.2 claim that construction is "a very serial process": with X a power of two and non-negative digits, each coefficient is an independent bit-slice.

### 2.2 🟡 `getMaxCoeff` is unnecessary for `ZZ_pX`

min(U,V)·(p−1)² is an exact a-priori bound, which is the "approximating" variant done correctly. Note min(U,V), not max(U,V).

### 2.3 🟡 Reduction mod xʳ − 1 and repeated modulus set-up

- `ZZpPowMod` uses general polynomial division (`%=`). Reducing modulo xʳ − 1 is a wrap-around: add coefficient i + r into coefficient i, which is O(r).
- In `CongruenceZnx*`, `PowerMod(c, n, b)` is called for every *a* with a plain `ZZ_pX` modulus. NTL's `PowerMod` only takes a `ZZ_pXModulus`, so a temporary one, with its pre-computation, is built on every call. Build `ZZ_pXModulus B(b)` once and pass `B`.
- The right-hand side xⁿ mod (xʳ − 1) is computed with a full `PowerMod` (`ZZpPowMod` in the Z version). It is simply the single term x^(n mod r).

### 2.4 🟡 The `r` search

Even after the fix in §1.2, the `r` search should work with word-sized `n % r` arithmetic. It currently calls `ZZ_p::init(r)` and does big-integer work inside the loop.

### 2.5 🟡 φ and λ

These are computed by O(r) trial division, separately from each other. See §1.8.

### 2.6 🟡 I/O inside the hot loop

The code prints a line with `printf` for every value of *a*, including from the worker threads.

### 2.7 🟡 Configuration drift between versions

`LenstraZ.cpp` uses φ(r) and `LenstraZnx.cpp` uses λ(r), so whole-program comparisons between the Z and Zn[x] versions are not like-for-like. Make the variant a runtime option.

### 2.8 🟡 Global state at static initialisation

The log file's `std::ofstream` is opened during static initialisation with a timestamped name, so a file is created even when the program exits early. The CSV header is written on every run, into a fresh file each time.

---

## 3. Dissertation assertions that carry into the code

| Where | Assertion | Problem | Code impact |
|---|---|---|---|
| p. 4, §2.1.1 | "the largest factor of n must be n/2", so search only to √n | The two halves are inconsistent. Trial division needs only primes ≤ √n. | `Carmichael.h` `isPrime` / `primeFactors` loop to n/2 (§1.8) |
| p. 6, Eq. 2.1.5; Alg. 5.2.3 lines 7, 9 | λ(p^α) = φ(**n**) or ½φ(**n**) | Should be φ(p^α). | Code is correct (`Phi(power_long(a,b))`), but the pseudocode does not match it |
| p. 7, §2.1.2 | "log n" denotes log₂ n | The code uses NTL `log`, the natural log. | `r` threshold about 2× too small (§1.2) |
| p. 16, Alg. 2.4.1 lines 8–9 | "if n ≤ r then return NOT PRIME" | Inverted: AKS outputs PRIME here. | Step missing from code (§1.4) |
| p. 16, Lemma 2.4.1 vs p. 17, Alg. 2.4.2 | a ≤ ⌈√r·log n⌉ vs a ∈ [1, φ(r) − 1] (supervisor, p. 16) | Inconsistent; the standard bound is ⌊√φ(r)·log₂ n⌋. | Code tests φ(r)−1 or λ(r)−1 (§1.3) |
| p. 17, §2.4.2; §5.2.4; App. A.2 | λ(r) − 1 iterations suffice (Egan, 2005) and give up to 36× speed-ups | Not supported by the theorem. The gains come from composite, invalid r; at a valid (typically prime) r, λ = φ. | `LenstraZnx.cpp` uses λ (§1.3) |
| p. 37, Alg. 5.1.2 | 2^b − 1 > max{U,V}·max{f_j}·max{g_k}; substitute 2^b − 1 (supervisor, p. 37) | No room for the sign bit; it must be 2^(b−1) > … and X = 2^b. | **Fixed on branch** (§1.1) |
| p. 37, Alg. 5.1.2 | max{f_j}, max{g_k} for signed coefficients | The bound needs max\|f_j\|. | **Fixed on branch** (§1.1) |
| p. 38, Alg. 5.1.3 line 6 | a(x) := a(x) × a(x) | Missing mod b (supervisor, p. 38). | **Fixed on branch** (§1.1) |
| p. 37, §5.1.2 | Binary segmentation is O(D_s) | As implemented it is O(D²·b); the packed version is O(M(D·b)). | §2.1 |
| §5.2.2 | Construction "is a very serial process" with "no data-level parallelism" | With X = 2^b each digit is an independent bit-slice. | §2.1 |
| §5.2.2–5.2.3, §7.3–7.4 | Overhead factor of 14; several orders of magnitude slower than NTL; Z version "redundant" | An artefact of the quadratic implementation; the packed version runs within 1.3–3.3× of NTL. | §2.1 |
| §5.2.1, Alg. 5.2.1 | Powermod runs in "O(log n) time" | That is O(log n) *multiplications*; time is O(log n · M(r·log n)). | — |
| §5.2.5 | Passing `const long&` avoids copying; prefix `++` is faster for `long`; `inline` is ignored for non-void return types; `printf` is 12× faster than `cout` | Scalars are cheaper to pass by value. `++i` and `i++` compile identically for built-in types. Return type has no effect on inlining. The `printf` gap comes from stream synchronisation and is irrelevant here. | Signatures such as `const long& a`, `const bool& PRIME`, `const unsigned int& cores`; `inline bool Lenstra` |
| §6.2.1 | A vector with `push_back` is used to collect results from the threads | Unsynchronised concurrent writes. | Data race (§1.6) |
| §4.3.5 | `high_resolution_clock` used for timing | Not guaranteed to be steady (on libstdc++ it is `system_clock`). | Use `steady_clock` for intervals |
| App. C.1 | Everything compiles with the given `g++ -std=c++11 … -march=native` command | Listing C.1 includes both congruence headers, which is a redefinition error. On Apple Silicon, `-march=native` and `/usr/local` paths do not apply. | See `README.md` |
| §7.3 | Both implementations are "correct" | The binary segmentation was not correct (§1.1), and the `r` chosen does not satisfy the theorem (§1.2), so PRIME results were not certified. | — |

---

## 4. Suggested order of work

1. Merge `fix/binary-segmentation-sign-bit` after reviewing the WIP commit (§1.11).
2. Make the AKS parameters correct (§1.2–1.4): log₂, true order test, Lemma bound for *a*, the n ≤ r step. Keep φ/λ as labelled experimental modes.
3. Replace `PerfectPower` with exact integer roots (§1.5).
4. Fix the data race and add early termination (§1.6).
5. Make headers self-contained, with a `Makefile` or CMake build and a `ctest` suite that runs `biSegCorrectnessTest` plus a primality sweep against `ProbPrime` (§1.9–1.10).
6. Replace the quadratic binary segmentation with Kronecker packing (§2.1) and re-run the Chapter 5 comparison. Plot time against log n on log–log axes, as the supervisor suggested (p. 45).
7. Remaining efficiency items (§2.2–2.8).

These map directly onto the Rust plan in the AKS2 repository.
