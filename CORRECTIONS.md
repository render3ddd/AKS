# CORRECTIONS.md

Changes on `fix/binary-segmentation-sign-bit` relative to `main` (`e173feb`). The branch contains:

1. **The main fix** for the binary segmentation bug the supervisor identified (`3242448`), below.
2. **Other code changes** from the 2024 work-in-progress snapshot (`2e78c69`) and tidy-ups made alongside the fix.
3. **Eighteen refactor pull requests** (#2–#22): one per issue in `REFACTOR.md`, plus two follow-ups (quick fixes and dedicated squaring), listed in [Refactor pull requests](#refactor-pull-requests). Together they change the programs' command line, defaults and output; see [Behaviour changes](#behaviour-changes-for-anyone-running-the-programs).

## Main fix: binary segmentation (supervisor's feedback, pp. 37–38)

The marker's notes on the dissertation identified two defects in the binary segmentation algorithms (Algorithms 5.1.2 and 5.1.3). Both were also present in `devVR/biSegMultiplyZZpX.h` and `devVR/biSegMultiplyZZX.h`.

> **p. 37.** "I think we want 2^{b−1} > …, not 2^b − 1 (to allow for the sign bit), and you need to substitute in 2^b, not 2^b − 1."
> **p. 38.** "Line 6 also needs to be modulo b, which seems to be missing here."

### Errors

1. **Bound and evaluation point.** The code chose *b* with 2^b − 1 > max(U,V)·max f_j·max g_k and substituted X = 2^b − 1. Signed-digit recovery treats a digit greater than X/2 as negative, so it is only exact when every product coefficient satisfies |s_i| < X/2. With X only just above the bound, legitimate positive coefficients in (X/2, bound] were turned into negatives. Initial work changed X to 2^(b−1), which still leaves no sign margin, and it updated only two of the four copies of the loop.
2. **Comparing after reduction** (`ZZpXmultiply`). Each digit was reduced mod p *before* it was compared with X/2. Any true coefficient in (X/2, p) was then corrupted. Example: (1 + x + … + x⁸)² mod 10007.
3. **Unsigned bound for signed coefficients** (`getMaxCoeff(ZZX)`). The function returned the largest signed coefficient, so for x − 5 the "maximum" was 1 and the bound was far too small. (x − 5)², (−7x − 3)² and even 3·3 came out wrong.
4. **Missing `mod b`** (`polyPowMod`, line 6 of Algorithm 5.1.3). a(x)² was never reduced, so a(x) grew without bound. The wrong products eventually produced a zero polynomial.
5. **Zero operand.** A zero polynomial has deg = −1, which gave a zero-length array, and the code then wrote past its end. This produced the segfault in item 4.

### Changes

| File | Change |
|---|---|
| `devVR/biSegMultiplyZZpX.h` | Choose *b* with 2^(b−1) > bound and evaluate at **X = 2^b**. Recover each signed digit on the integer value, then reduce mod p. Return 0 early for a zero operand. |
| `devVR/biSegMultiplyZZX.h` | Same bound, evaluation point and zero handling. `getMaxCoeff` now returns max \|f_i\|. `polyPowMod` reduces a(x) mod b(x) before the loop and after every squaring. |
| `test/biSegCorrectnessTest.cpp` (new) | 1,608 checks of `polyMultiply`, `polyPowMod`, `ZZpXmultiply` and `ZZpPowMod` against NTL's own arithmetic, including the edge cases above. |

### Verification

- **Old code:** fails 35 of the 1,561 checks it can run, and segfaults once the `polyPowMod` and zero-polynomial sections are included.
- **Fixed code:** passes **1,608 / 1,608** with GCC and Clang, clean under AddressSanitizer and UBSan.
- **End to end:** the Z version of AKS agrees with a reference primality test for every n from 2 to 699, and for n = 10007, 23456611 and 23456789.

```bash
clang++ -O2 -std=c++17 -pthread test/biSegCorrectnessTest.cpp -o test/biSegCorrectnessTest.out -lntl -lgmp -lm
./test/biSegCorrectnessTest.out      # 1608 / 1608 at the time of this fix; 2322 / 2322 since #22 added squaring checks
```

---

## Other code changes

These changes are not part of the bug fix. They come from the WIP commit (`2e78c69`), plus the tidy-ups made alongside the fix in `3242448`.

### Semantics: changes in behaviour or build structure

| File | Change | Note |
|---|---|---|
| `devVR/LenstraZnx.cpp` | `main()` reads *n* from standard input instead of looping over the hard-coded `nos[]` array | The leftover `nos`, `nosSize`, `nosEnd` and `prime` were removed in #12; #19 then replaced `main()` |
| `devVR/LenstraZ.cpp`, `LenstraZnx.cpp` | Standard-library includes trimmed: `<unistd.h>`, `<filesystem>`, `<sys/*.h>`, the duplicate C headers, and **`<chrono>`** removed. NTL includes moved after the project headers and reduced to `ZZ`, `ZZ_p` and `ZZX`. | `<chrono>` restored in #12 |
| `devVR/PerfectPower.h` | Now includes the NTL headers it needs; names qualified with `NTL::` | Lets the header be included before `NTL_CLIENT` |
| `devVR/Euler.h` | `ZZ` / `to_ZZ` qualified as `NTL::ZZ` / `NTL::to_ZZ` | Same reason |
| `devVR/LenstraZ.cu`, `devVR/LenstraZnx.cu` | **Deleted** | The CUDA kernels could not work with NTL types; the stale `launch.json` entry was removed in #12 |
| `devVR/LenstraZnx.cpp` | `n == 2` and `n == 3` branches merged and printed with `%d` / `to_int(n)` | Same output |
| `devVR/Lenstra*.cpp` | `filename` made `const` | — |
| `devVR/Lenstra*.cpp` | `printf("%ld", time)` changed to `%lld` | `static_cast<long long>` added in #8 so it is correct on both |
| `devVR/biSegMultiplyZZX.h` | `polyMulMod` and `polyMulModN` now call `polyMultiply` instead of carrying their own copies of the algorithm | Same results; removes three diverging copies |
| `devVR/biSegMultiplyZZX.h` | `polyPowMod(ZZX a, ZZ n, ZZX b)` now takes `const ZZX& b` | Source-compatible for callers |
| `devVR/biSegMultiplyZZ*X.h` | `std::vector` replaces variable-length arrays (`ZZ s[fgTrm]`) | VLAs are a non-standard compiler extension; Clang warns about them |

### Formatting and other stylistic changes

| File | Change |
|---|---|
| `devVR/Carmichael.h`, `CongruenceZnxHPC.h` | Index loops replaced by range-based `for`; `test.size() == 0` became `test.empty()` |
| `devVR/CongruenceZ.h`, `PerfectPower.h` | `return(j);` / `return(0);` became `return j;` / `return 0;`; one `ZZ_pX` declaration became `auto` |
| `devVR/biSegMultiplyZZpX.h` | `#include <math.h>` became `<cmath>`; `ZZ xpow` became `auto xpow` |
| `devVR/Lenstra*.cpp` | Unreachable `break;` after `return false;` removed, plus a commented-out `// break;` |
| `devVR/Lenstra*.cpp` | Commented-out `mkdir` / `create_directory` / alternative-folder lines removed |
| `devVR/Lenstra*.cpp` | Spaces after commas in `fileWrite(...)` and `printf(...)` arguments |
| `devVR/biSegMultiply*.h` | Comments updated to the corrected formulas (2^(b−1) > bound, X = 2^b) |

---

## Refactor pull requests

Each issue in `REFACTOR.md` was fixed on its own `refactor/<issue>-…` branch and merged here, in the order listed. Every merge was checked with `make test` (1,608 binary segmentation checks, 2,322 from #22 on, and an AKS sweep of n = 2…1000 against NTL's `ProbPrime`).

| PR | Issue | Change |
|---|---|---|
| #11 | 1.10 tests & build | `Makefile` (`make`, `make test`, Homebrew paths found automatically); `test/aksSweep.sh` + `test/probPrime.cpp`; `chronoTest.cpp` and `getHardwareTest.cpp` repaired; two stale `_cli.h` copies deleted. |
| #12 | 1.11 WIP review | `<chrono>` restored; unused variables removed; stale CUDA launch config removed. |
| #14 | 2.2 a-priori bound | Multiply size bound is min(U,V)·(p−1)², with no scan for the largest coefficient. |
| #18 | 2.6 quiet hot loops | Per-*a* progress output only with `-DAKS_VERBOSE=1`. |
| #6 | 1.5 exact perfect power | Integer k-th roots replace floating-point Newton; detects perfect powers of any size; no longer hangs. |
| #8 | 1.7 print large *n* | Large numbers printed and logged in full. |
| #9 | 1.8 Carmichael robustness | √n trial division, gcd-based LCM, edge cases. |
| #20 | 2.8 lazy log file | Log opened on first result; `logs/Z`, `logs/Znx` created if missing; header written once. |
| #7 | 1.6 thread-safe congruence | Data race removed (atomic minimum); interleaved *a* with early exit; reports the smallest failing *a*. |
| #2 | 1.2 valid r (and 1.4) | r chosen with ord_r(n) > ⌊log₂ n⌋², tested exactly in new `FindR.h`; "no r < n ⇒ PRIME" step added. |
| #19 | 2.7 runtime config | Same command line for both programs; *n* as argument or on stdin; exit status. |
| #4 | 1.3 a-range bound | New `ARange.h`; default range of *a* is the theorem's ⌊√φ(r)·log₂ n⌋; `--phi` / `--lambda` kept as "not certified" modes. |
| #17 | 2.5 single factorisation | φ(r) and λ(r) from one factorisation (`Totients.h`). |
| #13 | 2.1 Kronecker packing | Binary segmentation by bit packing: O(D²·b) → O(M(D·b)), within about 1–3× of NTL. |
| #16 | 2.3 cheap reductions | Reduce mod xʳ − 1 by folding; NTL modulus built once; xⁿ computed as x^(n mod r). |
| #10 | 1.9 header hygiene | `#pragma once`, self-contained includes, `inline`; `ncores` in new `Cores.h`. |
| #22 | §5 dedicated squaring | When both operands are the same polynomial, pack once and square the big integer; new `ZZpXsqr()` used by both powering loops. 714 new checks (2,322 in total). Squares 1.4–1.65× faster; about 8% end to end. |
| #21 | §5 quick fixes | `-Wall`-clean with GCC and Clang (cast in #2's guard; `unitTest.cpp` formats); `steady_clock`; scalars passed by value; `ioTest.cpp` gets `<cstdint>` and no longer writes past its arrays; `make` defaults to `c++`; `.vscode` set up for LLDB and Homebrew on Apple Silicon. |

Not merged: the 1.4 PR (already contained in #2) and the 2.4 PR (its search was replaced by #2).

### Behaviour changes for anyone running the programs

| Before | After |
|---|---|
| `LenstraZnx.out` looped over hard-coded 42-digit primes; `LenstraZ.out` tested a hard-coded n | Both: `<program> [--theorem \| --phi \| --lambda] [n]`; *n* is read from standard input if not given |
| Zn[x] tested a ≤ λ(r) − 1, Z tested a ≤ φ(r) − 1 | Both default to a ≤ ⌊√φ(r)·log₂ n⌋ (`--theorem`); `--lambda` / `--phi` reproduce the dissertation's runs and are marked "not certified" |
| r chosen by an invalid search (e.g. 256 for 23456789) | Smallest r with ord_r(n) > ⌊log₂ n⌋² (599 for 23456789), so runs do more work than the dissertation's but are backed by the theorem |
| Exit status always 0 | 0 if *n* is prime, 1 otherwise |
| A CSV with a header was created on every run, even with no result; silently nothing if `logs/…` was missing | CSV created on the first result, folder created if needed; notes record r, φ(r), λ(r) and the range mode |
| A line printed for every *a* tested | Only with `-DAKS_VERBOSE=1` |
| Z version on 23456789 at the valid r = 599: about 125 s with the quadratic multiply | About 2.7 s (#13 and #16) |
