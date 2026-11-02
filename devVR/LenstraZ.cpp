/*
    A C++ implementation of Lenstra's algorithm, adapted from that from Hua Li:
        https://researchportal.bath.ac.uk/en/publications/the-analysis-and-implementation-of-the-aks-algorithm-and-its-impr
    Compile with:
        $ g++ -g -O2 -std=c++11 -pthread -march=native devVR/LenstraZ.cpp -o devVR/LenstraZ.out -lntl -lgmp -lm
*/

#include <cmath> // standard libraries
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <ctime>
#include <chrono>
#include <string>
#include <thread>
#include <array>
#include <sys/stat.h>
#include <cstring>

#include "PerfectPower.h" //Each Independent Test
#include "Euler.h"
#include "FindR.h"
#include "ARange.h"
#include "Carmichael.h"
#include "CongruenceZ.h"

#include "NTL/ZZ.h" // NTL Libraries
#include "NTL/ZZ_p.h"
#include "NTL/ZZX.h"
NTL_CLIENT

// Decimal string of an arbitrary-size integer. (std::to_string(to_long(x))
// truncated anything above 2^63, e.g. the 42-digit test primes.)
std::string ZZtoString(const ZZ& x) {
    std::ostringstream oss;
    oss << x;
    return oss.str();
}

std::string getDate() {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d");
    auto date = oss.str();

    return date;
}

std::string getTime() {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%H-%M-%S");
    auto time = oss.str();

    return time;
}

std::string getDateTime() {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d-%H-%M-%S");
    auto datetime = oss.str();

    return datetime;
}

std::string getFilename() {
    std::string fldr = "logs/Z/";
    std::string prfx = "LnstrZ-";
    std::string sffx = getDateTime();
    std::string extn = ".csv";

    std::string filename = fldr + prfx + sffx + extn;

    return filename;
}

// REFACTOR 2.7 / 1.3: how far the congruence tests run is a runtime option,
// shared by the Z and Zn[x] programs so that whole-program comparisons are
// like-for-like (see ARange.h):
//   default / --theorem : 1 <= a <= floor(sqrt(phi(r)) log2 n)   (AKS/Lenstra bound)
//   --phi               : 1 <= a <= phi(r) - 1      (Algorithm 2.4.2; not certified)
//   --lambda            : 1 <= a <= lambda(r) - 1   (Egan, 2005; not certified)
// This program's dissertation runs used --phi.
ARange aRange = ARange::Theorem;

// REFACTOR 2.8: the log used to be a global std::ofstream opened during
// static initialisation, so every run created a timestamped CSV (with just a
// header) even when it stopped early, and failed silently if logs/Z/ was
// missing. It is now opened on the first result: the folder is created if
// needed, a failure is reported, and the header is written only to a new file.
std::ofstream& perfLog() {
    static std::ofstream log;
    if (!log.is_open()) {
        mkdir("logs", 0777);                    // no-op if it already exists
        mkdir("logs/Z", 0777);
        const std::string filename = getFilename();
        const bool isNew = !std::ifstream(filename).good();
        log.open(filename, std::ios::app);
        if (!log) {
            std::fprintf(stderr, "warning: cannot open log file %s\n", filename.c_str());
        }
        else if (isNew) {
            log << "Int, Cores, Prime (T/F), Time (milliseconds), Comments\n";
        }
    }
    return log;
}

inline void fileWrite(const ZZ& n, const unsigned int& cores, const bool& PRIME, const long& time, const std::string& other) {
    perfLog() << n << "," << cores << "," << PRIME  << "," << time << "," << other << "\n";
}

inline bool Lenstra (const ZZ& n) {
    if(n < 1){
        std::printf("Integer n needs to be positive.\n");
        return false;
    }
    else if(n == 1){
        std::printf("1 is neither prime or composite.\n");
        return false;
    }
    else if(n == 2){
        std::printf("2 is prime.\n");
        return true;
    }
    else if(n == 3){
        std::printf("3 is prime.\n");
        return true;
    }

    std::cout << "n = " << n << "\n\n";

    // start timing
    auto start = std::chrono::high_resolution_clock::now();

    // Test if n is a perfect power
    int PP = PerfectPower(n);

    // returns 1 if n is a perfect power, 0 otherwise;
    if(PP == 1){
        auto finish = std::chrono::high_resolution_clock::now();
        auto duration = finish - start;
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();

        std::cout << n << " is not prime.\n";
        std::cout << n << " is a perfect power.\n";
        std::printf("Time taken: %lld milliseconds\n", static_cast<long long>(time));

        std::string note = ZZtoString(n) + " is a perfect power";
        fileWrite(n, ncores, false, time, note);

        return false;
    }

    // Find a suitable r: the smallest r with ord_r(n) > floor(log2 n)^2.
    // (Previously: natural log, truncated before squaring, an inverted order
    // test, and n^v computed in a 64-bit long, which overflowed.)
    const long log2n = NumBits(n) - 1;          // floor(log2 n)
    const long maxOrder = log2n * log2n;        // need ord_r(n) > floor(log2 n)^2
    ZZ r = to_ZZ(2);
    ZZ r1;                                       // 0 until a suitable r is found

    while(r < n){
        ZZ R = GCD(r, n);
        if(R != 1 ){
            auto finish = std::chrono::high_resolution_clock::now();
            auto duration = finish - start;
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();

            std::cout << n << " is not prime.\n";
            std::cout << R << " is a divisor.\n";
            std::printf("Time taken: %lld milliseconds\n", static_cast<long long>(time));

            std::string note = ZZtoString(R) + " is a divisor";
            fileWrite(n, ncores, false, time, note);

            return false;
        }
        else if (OrderExceeds(rem(n, to_long(r)), to_long(r), maxOrder)) {
            r1 = r;                              // store value of r
            break;
        }
        r = r + 1;
    }

    if (IsZero(r1)) {
        // No suitable r below n, and gcd(a, n) = 1 for every 2 <= a < n,
        // so n is prime (AKS step: "if n <= r, output PRIME"). Without this
        // guard r would be 0 below.
        auto finish = std::chrono::high_resolution_clock::now();
        auto duration = finish - start;
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();

        std::printf("%ld is prime.\n",to_long(n));
        std::printf("Time taken: %lld milliseconds\n", time);

        fileWrite(n, ncores, true, time, "no r < n with ord_r(n) > log2(n)^2; all a < n coprime to n");

        return true;
    }

    r = r1;
    std::printf("r = %ld\n",to_long(r));

    ZZ r2 = Euler(to_long(r));                                  // phi(r)
    long lambda = to_long(Carmichael(to_long(r)));              // lambda(r)
    std::printf("Phi(%ld) = %ld, Lambda(%ld) = %ld\n", to_long(r), to_long(r2), to_long(r), lambda);

    long a = ABound(aRange, n, to_long(r2), lambda);
    std::printf("Testing 1 <= a <= %ld (%s)\n", a, ARangeName(aRange));
    long f = CongruenceZ(n,r,r2,a);

    if(f == 0){
        auto finish = std::chrono::high_resolution_clock::now();
        auto duration = finish - start;
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();

        std::cout << n << " is prime.\n";
        std::printf("Time taken: %lld milliseconds\n", static_cast<long long>(time));

        std::string note = "a = " + std::to_string(a) + " (" + ARangeName(aRange) + "); End: a = " + std::to_string(f) + "; r = " + std::to_string(to_long(r)) + "; phi(r) = " + std::to_string(to_long(r2)) + "; lambda(r) = " + std::to_string(lambda);
        fileWrite(n,ncores,true, time, note);

        return true;
    }
    else {
        auto finish = std::chrono::high_resolution_clock::now();
        auto duration = finish - start;
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();

        std::cout << n << " is not prime.\n";
        std::printf("The a which fails is %ld\n",f);
        std::printf("Time taken: %lld milliseconds\n", static_cast<long long>(time));

        std::string note = "a = " + std::to_string(a) + " (" + ARangeName(aRange) + "); End: a = " + std::to_string(f) + "; r = " + std::to_string(to_long(r)) + "; phi(r) = " + std::to_string(to_long(r2)) + "; lambda(r) = " + std::to_string(lambda);
        fileWrite(n, ncores, false, time, note);

        return false;
    }

}

int main (int argc, char * argv[]) {
    // usage: <program> [--theorem | --phi | --lambda] [n]
    //   n is read from standard input when it is not given as an argument.
    ZZ n;
    bool haveN = false;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--theorem") == 0) {
            aRange = ARange::Theorem;
        }
        else if (std::strcmp(argv[i], "--phi") == 0) {
            aRange = ARange::PhiMinusOne;
        }
        else if (std::strcmp(argv[i], "--lambda") == 0) {
            aRange = ARange::LambdaMinusOne;
        }
        else {
            n = conv<ZZ>(argv[i]);
            haveN = true;
        }
    }

    if (!haveN) {
        std::printf("Enter a positive integer number n you want to be tested:\n");
        std::cin >> n;
    }

    return Lenstra(n) ? 0 : 1;
}
