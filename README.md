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

Optionally, confirm that NTL was built with the GMP backend and thread boost, which `CongruenceZnxHPC.h` relies on for its parallel `NTL_EXEC_RANGE` loop:

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

Apple clang does not search `/opt/homebrew` by default. Either pass `-I`/`-L` flags on every build, or add these lines to `~/.zshrc` once so that plain `-lntl -lgmp` works:

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

Run all commands from the **repository root**. The programs write their CSV logs to the relative paths `logs/Z/` and `logs/Znx/`.

```bash
# Zn[x] version (NTL PowerMod, parallel congruence tests)
clang++ -g -O2 -std=c++17 -pthread devVR/LenstraZnx.cpp -o devVR/LenstraZnx.out -lntl -lgmp -lm

# Z version (binary segmentation)
clang++ -g -O2 -std=c++17 -pthread devVR/LenstraZ.cpp -o devVR/LenstraZ.out -lntl -lgmp -lm
```

Changes from the dissertation's command:

- `-march=native` is dropped. NTL and GMP already use Homebrew's arm64-tuned code, and the flag is aimed at x86 builds.
- `-std=c++17` instead of `-std=c++11`. Both compile, but C++17 matches the `MacOS` profile in `.vscode/c_cpp_properties.json`.

Expected compiler warnings: `-Wformat` (`%ld` with a `long long` chrono count), unused variables, and on `main` a `-Wvla-cxx-extension` warning from `biSegMultiplyZZpX.h`. They are harmless for building. The fix branch removes the VLA.

### 7. Run

```bash
./devVR/LenstraZnx.out
./devVR/LenstraZ.out
```

Check what each `main()` does before you run it:

- On `main`, `LenstraZnx.cpp` loops over the four hard-coded 42-digit primes in `nos[]`. The dissertation reports about 3 h 45 min for **each one** in serial. For a quick smoke test, edit `nos[]` or use the fix branch.
- On the fix branch, `LenstraZnx.out` reads *n* from standard input:

  ```bash
  echo 23456789 | ./devVR/LenstraZnx.out      # prints "23456789 is prime." in a few seconds
  ```

- `LenstraZ.cpp` (both branches) tests the hard-coded `n = 1003026954441971`. The binary segmentation version is very slow, so change `p` in `main()` to something small such as `23456789` for a first run.

`CongruenceZnxHPC.h` sets the NTL thread pool to `std::thread::hardware_concurrency()`. On an M3 Pro that is 11 or 12, depending on the model, mixing performance and efficiency cores. Lower `ncores` when timing runs, to keep benchmarks consistent.

### Running the correctness test

```bash
clang++ -O2 -std=c++17 -pthread test/biSegCorrectnessTest.cpp -o test/biSegCorrectnessTest.out -lntl -lgmp -lm
./test/biSegCorrectnessTest.out             # expect: "1608 / 1608 checks passed", exit code 0
```

This checks every binary segmentation routine against NTL's own arithmetic.

### VS Code Setup (optional)

1. Install the **C/C++** (`ms-vscode.cpptools`) and **CodeLLDB** (`vadimcn.vscode-lldb`) extensions.
2. In `.vscode/c_cpp_properties.json`, `MacOS` configuration: change `"/usr/local/include/**"` to `"/opt/homebrew/include/**"` (the Intel-Mac Homebrew path will not find NTL on Apple Silicon), and set `"compilerPath": "/usr/bin/clang++"`.
3. `.vscode/launch.json` uses `"MIMode": "gdb"`, but GDB does not run on Apple Silicon. Use LLDB instead, for example:

    ```json
    {
        "name": "LenstraZnx (lldb)",
        "type": "lldb",
        "request": "launch",
        "program": "${workspaceFolder}/devVR/LenstraZnx.out",
        "cwd": "${workspaceFolder}"
    }
    ```

    Set `cwd` to the repository root (not `${fileDirname}`) so the relative `logs/` paths resolve.

### What does not build or run on a Mac

- **CUDA.** NVIDIA's CUDA toolkit has not supported macOS since CUDA 10.2, and Apple Silicon has no NVIDIA GPU, so `CUDA/` and the `.cu` files are reference only. GPU work needs a Linux machine with an NVIDIA card. On Apple Silicon the native GPU route would be Metal, which is a separate port.
- **`test/getHardwareTest.cpp`** includes the x86-only `<cpuid.h>`.
- **`test/chronoTest.cpp`** is stale: it calls functions with signatures that no longer exist, and fails to compile on any platform.
- **`GallotY/AKSgallotMIRACL.cpp`** needs MIRACL, which is not packaged by Homebrew. The file is also in `.gitignore`.

### Troubleshooting

| Symptom | Fix |
| --- | --- |
| `fatal error: 'NTL/ZZ.h' file not found` | Step 4 (`CPATH`) not applied, or pass `-I"$(brew --prefix)/include"`. |
| `ld: library 'ntl' not found` | Step 4 (`LIBRARY_PATH`) not applied, or pass `-L"$(brew --prefix)/lib"`. |
| `ld: warning: ignoring file ... built for macOS-x86_64` | Mixed Intel/arm64 toolchains: `uname -m` must be `arm64` and `brew --prefix` must be `/opt/homebrew`. |
| Program runs but no CSV appears | Run from the repo root; `logs/Z` and `logs/Znx` must exist. |
| "0 of 0 threads available" | `hardware_concurrency()` returned 0; set `ncores` manually in `CongruenceZnxHPC.h`. |

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

On Debian/Ubuntu, `sudo apt install build-essential libgmp-dev libntl-dev` provides everything. The commands in steps 6–8 work unchanged with `g++` or `clang++`.
