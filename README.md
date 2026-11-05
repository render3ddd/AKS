# **University of Bath Computer Science MSc. Dissertation Project**

## **On the implementation of the AKS algorithm in CUDA**

An implementation of the Lenstra variant of the AKS primality test and a further, unimplemented variant suggestion from [Crandall & Papadopoulos (2003)](https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.95.58).

### Stages of Development

- [ ] Review and validate (respectively) previous implementations from:
  - [ ]  [Gallot](http://yves.gallot.pagesperso-orange.fr/src/)
  - [ ]  [Hua Li](https://researchportal.bath.ac.uk/en/publications/the-analysis-and-implementation-of-the-aks-algorithm-and-its-impr) - Lenstra/$Z_n(x)$ version
- [ ] Implement Lenstra/$Z_n(x)$ version as GPU accelerated versions in CUDA C++ for deployment on a GEFORCE RTX 3090 GPU cloud.
- [ ] Develop $Z_n$ version from [Crandall & Papadopoulos (2003)](https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.95.58)
- [ ] Optimise above implementation, using MIRACL or NTL

### Other information

Personal System Specifications for home benchmarking purposes:

- Eight-core AMD Ryzen 7 4700U, 8GB RAM
- Dual-core Intel(R) Pentium(R) CPU G860 @ 3.00GHz, 16GB DDR3 RAM @ 1333MHz, Radeon HD 7770 GPU

---

## Local Development Environment — Apple Silicon macOS

These steps recreate the development environment on an **Apple M3 Pro MacBook running macOS 27.0.1 "Golden Gate"**. They should also work on any Apple Silicon Mac on a recent macOS.

### Outputs

| Component | Version | Source |
| --- | --- | --- |
| Compiler | Apple clang (Xcode Command Line Tools) | `xcode-select --install` |
| Package manager | Homebrew, installed at `/opt/homebrew` | [brew.sh](https://brew.sh) |
| GMP | 6.3.0 | `brew install gmp` |
| NTL | 11.6.0, built against GMP with threads and thread boost on | `brew install ntl` |
| Editor (optional) | VS Code with the *C/C++* and *CodeLLDB* extensions | — |

> **Why not the dissertation's NTL 11.5.1 + GMP 6.2.0?** GMP 6.2.0 predates Apple Silicon support. GMP added Apple Arm support in 6.2.1, with a known register (`x18`) caveat, and the GMP maintainers recommend 6.3.0 on M-series Macs. Homebrew's NTL 11.6.0 is built against GMP 6.3.0, so the version clash described in Appendix C.1 of the dissertation does not apply. If you need the exact NTL 11.5.1 for reproducibility, see [Building NTL from source](#optional-building-ntl-from-source-pinned-version) below.

### 1. Install the Xcode Command Line Tools

```bash
xcode-select --install          # skip if already installed
clang++ --version               # should report "Apple clang" and "Target: arm64-apple-darwin..."
uname -m                        # must print arm64 (if it prints x86_64 your terminal is running under Rosetta)
```

On macOS, `g++` is a shim for Apple clang, so the `g++ ...` commands in the dissertation also run clang. This guide uses `clang++` explicitly.

### 2. Install Homebrew

Follow the one-line installer at [brew.sh](https://brew.sh), then make sure `brew` is on your `PATH` (the installer prints these lines):

```bash
echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zprofile
eval "$(/opt/homebrew/bin/brew shellenv)"
brew --prefix                   # should print /opt/homebrew
```

### 3. Install GMP and NTL

```bash
brew install gmp ntl
brew list --versions gmp ntl    # expect gmp 6.3.0 and ntl 11.6.0 (or newer)
```

Optionally, confirm that NTL was built with the GMP backend and thread boost, which `CongruenceZnxHPC.h` relies on for its parallel `NTL_EXEC_INDEX` loop:

```bash
cat > /tmp/ntlcheck.cpp <<'EOF'
#include <NTL/version.h>
#include <NTL/ZZ.h>
#include <cstdio>
int main() {
    std::printf("NTL %s\n", NTL_VERSION);
#ifdef NTL_GMP_LIP
    std::puts("GMP backend: on");
#else
    std::puts("GMP backend: OFF");
#endif
#ifdef NTL_THREAD_BOOST
    std::puts("Thread boost: on");
#else
    std::puts("Thread boost: OFF");
#endif
}
EOF
clang++ -std=c++17 -I"$(brew --prefix)/include" -L"$(brew --prefix)/lib" /tmp/ntlcheck.cpp -o /tmp/ntlcheck -lntl -lgmp && /tmp/ntlcheck
```

### 4. Point the compiler at Homebrew

Apple clang does not search `/opt/homebrew` by default. The `Makefile` (step 6) adds the Homebrew paths itself, so this step is only needed for building by hand. Either pass `-I`/`-L` flags on every build, or add these lines to `~/.zshrc` once so that plain `-lntl -lgmp` works:

```bash
export CPATH="$(brew --prefix)/include${CPATH:+:$CPATH}"
export LIBRARY_PATH="$(brew --prefix)/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"
```

Then open a new terminal or run `source ~/.zshrc`. The commands below assume you have done this.

### 5. Get the code

```bash
git clone https://github.com/render3d/AKS.git
cd AKS
```

### 6. Build

Run all commands from the **repository root**.

```bash
make            # builds devVR/LenstraZnx.out and devVR/LenstraZ.out
make test       # 2,322 binary segmentation checks + an AKS sweep of n = 2..1000 against NTL's ProbPrime
make clean
```

`make` compiles with `c++` (Apple clang on macOS); pass `CXX=clang++`, `CXX=g++` or any other compiler to choose. To build by hand instead (after step 4):

```bash
# Zn[x] version (NTL PowerMod, parallel congruence tests)
clang++ -g -O2 -std=c++17 -pthread devVR/LenstraZnx.cpp -o devVR/LenstraZnx.out -lntl -lgmp -lm

# Z version (binary segmentation)
clang++ -g -O2 -std=c++17 -pthread devVR/LenstraZ.cpp -o devVR/LenstraZ.out -lntl -lgmp -lm
```

Changes from the dissertation's command:

- `-march=native` is dropped. NTL and GMP already use Homebrew's arm64-tuned code, and the flag is aimed at x86 builds.
- `-std=c++17` instead of `-std=c++11`, matching the `MacOS` profile in `.vscode/c_cpp_properties.json`.
- Add `-DAKS_VERBOSE=1` to print a line for every value of *a* tested (off by default because it slows the parallel loop).

Both programs, and every program in `test/`, build without warnings with `-Wall` under Clang (on macOS and Linux) and GCC.

### 7. Run

Both programs take the same arguments:

```bash
./devVR/LenstraZnx.out [--theorem | --phi | --lambda] [n]
./devVR/LenstraZ.out   [--theorem | --phi | --lambda] [n]
```

- ***n*** is given as an argument, or read from standard input if omitted (`echo 23456789 | ./devVR/LenstraZnx.out`).
- **Range of *a*.** `--theorem` (the default) tests 1 ≤ a ≤ ⌊√φ(r)·log₂ n⌋, the bound the AKS/Lenstra proof requires. `--phi` (a ≤ φ(r) − 1) and `--lambda` (a ≤ λ(r) − 1) reproduce the dissertation's Z and Zn[x] runs respectively; they are marked "not certified" in the output because the proof does not cover them.
- **Exit status** is 0 if *n* is prime and 1 otherwise, so the programs can be used in scripts.
- **Logs.** Each run that reaches a result appends to a timestamped CSV in `logs/Znx/` or `logs/Z/` (created if missing), recording n, the core count, the verdict, the time in milliseconds, r, φ(r), λ(r) and the range mode.

Example:

```bash
./devVR/LenstraZnx.out 23456789
#   r = 599
#   Phi(599) = 598, Lambda(599) = 598
#   Testing 1 <= a <= 598 (theorem: floor(sqrt(phi(r)) log2 n))
#   ...
#   23456789 is prime.
```

Rough run times on two cores (an M3 Pro will be faster): about 2 s (Zn[x]) and 3 s (Z) for 23456789, and about 5 s and 7 s for 1000000007. A 16-digit prime such as 1003026954441971 takes about 100 s with Zn[x]. Time grows quickly with the number of digits: the 42-digit primes from the dissertation (listed in a comment in `LenstraZnx.cpp`) need hours. These times are longer than the dissertation's for the same *n*, because the r the code now chooses is valid (and larger); see `REFACTOR.md` §1.2.

The Zn[x] program runs its congruence tests on an NTL thread pool sized by `ncores` in `devVR/Cores.h` (`std::thread::hardware_concurrency()`, falling back to 1). On an M3 Pro that is 11 or 12, mixing performance and efficiency cores; lower it when timing runs, to keep benchmarks consistent.

### 8. Tests

`make test` runs:

- `test/biSegCorrectnessTest.cpp`: every binary segmentation routine checked against NTL's own arithmetic (expect `2322 / 2322 checks passed`).
- `test/aksSweep.sh`: `LenstraZnx.out` for every n from 2 to 1000, compared with NTL's `ProbPrime` (expect `aksSweep: 2..1000 all correct`). For a longer sweep, run `./test/aksSweep.sh 5000` after `make test` (which builds the reference checker it needs).

### VS Code Setup (optional)

The repository's `.vscode/` folder is already set up for Apple Silicon:

- `c_cpp_properties.json`: the `MacOS` profile looks for headers in `/opt/homebrew/include` and uses `/usr/bin/clang++`.
- `launch.json`: `(lldb)` entries for macOS and `(gdb)` entries for Linux (GDB does not run on Apple Silicon). Each runs from the repository root, so the relative `logs/` paths resolve; the AKS programs are started with `23456789` as a sample argument.

To use them:

1. Install the **C/C++** (`ms-vscode.cpptools`) and **CodeLLDB** (`vadimcn.vscode-lldb`) extensions.
2. Build first with `make` (and `make test/unitTest.out` for the unit-test entry); the launch entries run the built `.out` files.
3. Pick a configuration such as "LenstraZnx (lldb)" in the Run and Debug view.

### What does not build or run on a Mac

- **CUDA.** NVIDIA's CUDA toolkit has not supported macOS since CUDA 10.2, and Apple Silicon has no NVIDIA GPU, so `CUDA/` and the `.cu` files are reference only. GPU work needs a Linux machine with an NVIDIA card. On Apple Silicon the native GPU route would be Metal, which is a separate port.
- **`test/getHardwareTest.cpp`** builds, but the CPU vendor query uses x86-only CPUID, so on Apple Silicon it reports only the thread count and working directory.
- **`GallotY/AKSgallotMIRACL.cpp`** needs MIRACL, which is not packaged by Homebrew. The file is also in `.gitignore`.

### Troubleshooting

| Symptom | Fix |
| --- | --- |
| `fatal error: 'NTL/ZZ.h' file not found` | Step 4 (`CPATH`) not applied, or pass `-I"$(brew --prefix)/include"`. |
| `ld: library 'ntl' not found` | Step 4 (`LIBRARY_PATH`) not applied, or pass `-L"$(brew --prefix)/lib"`. |
| `ld: warning: ignoring file ... built for macOS-x86_64` | Mixed Intel/arm64 toolchains: `uname -m` must be `arm64` and `brew --prefix` must be `/opt/homebrew`. |
| Program runs but no CSV appears | Logs are written relative to the current folder: run from the repo root. A warning on stderr means the log file could not be created there. |
| "1 of 1 threads available" on a multi-core Mac | `hardware_concurrency()` returned 0, so `Cores.h` fell back to 1; set `ncores` manually in `devVR/Cores.h`. |

### Optional: building NTL from source (pinned version)

Only needed to reproduce the dissertation's NTL 11.5.1 exactly. GMP 6.3.0 from Homebrew is used as the backend.

```bash
brew install gmp
mkdir -p ~/src && cd ~/src
curl -LO https://libntl.org/ntl-11.5.1.tar.gz
tar xzf ntl-11.5.1.tar.gz && cd ntl-11.5.1/src
./configure PREFIX="$HOME/opt/ntl-11.5.1" GMP_PREFIX="$(brew --prefix gmp)" \
            CXX=clang++ NTL_THREADS=on NTL_THREAD_BOOST=on
make -j"$(sysctl -n hw.ncpu)" && make check && make install
```

Then build with `-I"$HOME/opt/ntl-11.5.1/include" -L"$HOME/opt/ntl-11.5.1/lib"` in place of the Homebrew NTL paths. Keep `$(brew --prefix)/lib` on the link path for GMP.

### Linux

On Debian/Ubuntu, `sudo apt install build-essential libgmp-dev libntl-dev` provides everything. The commands in steps 6–8 work unchanged, with `g++` or `clang++`.
