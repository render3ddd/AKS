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
#include <string>
#include <thread>
#include <array>
#include <cstring>

#include "PerfectPower.h" //Each Independent Test
#include "Euler.h"
#include "Carmichael.h"
#include "CongruenceZ.h"

#include "NTL/ZZ.h" // NTL Libraries
#include "NTL/ZZ_p.h"
#include "NTL/ZZX.h"
NTL_CLIENT

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

// REFACTOR 2.7: which totient bounds the congruence tests is a runtime option,
// shared by the Z and Zn[x] programs, so that whole-program comparisons are
// like-for-like (previously Z used phi(r) and Zn[x] lambda(r), hard-coded).
//   default / --phi : 1 <= a <= phi(r) - 1      (Algorithm 2.4.2)
//   --lambda        : 1 <= a <= lambda(r) - 1   (Egan, 2005)
bool useLambda = false;

const std::string filename = getFilename();
std::ofstream perflog(filename, std::ios::app); // output result into file

inline void fileWrite(const ZZ& n, const unsigned int& cores, const bool& PRIME, const long& time, const std::string& other) {
    perflog << n << "," << cores << "," << PRIME  << "," << time << "," << other << "\n";
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

        std::printf("%ld is not prime.\n",to_long(n));
        std::printf("%ld is a perfect power.\n",to_long(n));
        std::printf("Time taken: %lld milliseconds\n", time);

        std::string note = std::to_string(to_long(n)) + " is a perfect power";
        fileWrite(n, ncores, false, time, note);

        return false;
    }

    // Find a suitable r
    ZZ r = to_ZZ(2);
    ZZ R;
    ZZ r1;

    while(r < n){
        ZZ R = GCD(r, n);
        if(R != 1 ){
            auto finish = std::chrono::high_resolution_clock::now();
            auto duration = finish - start;
            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();

            std::printf("%ld is not prime.\n",to_long(n));
            std::printf("%ld is a divisor.\n",to_long(R));
            std::printf("Time taken: %lld milliseconds\n", time);

            std::string note = std::to_string(to_long(R)) + " is a divisor";
            fileWrite(n, ncores,false, time, note);

            return false;
        }
        else {
            ZZ v = to_ZZ(floor(power_long(to_long(log(n)), 2)));

            // order of n mod r is bigger than v;
            int p = 0;
            ZZ_p::init(r); // calculate mod r

            while(v <= r){
                ZZ x = to_ZZ(power_long(to_long(n), to_long(v))); // calculates x = n^v
                ZZ_p z = to_ZZ_p(x);
                if(z == to_ZZ_p(1)){
                    r1 = r; // store value of r
                    r = n + 1;
                    break;
                }
                else{
                    v = v + 1;
                }
            }
        }
        r = r + 1;
    }

    r = r1;
    std::printf("r = %ld\n",to_long(r));

    ZZ r2 = useLambda ? Carmichael(to_long(r)) : Euler(to_long(r));
    const std::string r2name = useLambda ? "lambda" : "phi";
    std::printf("%s(%ld) = %ld\n", r2name.c_str(), to_long(r), to_long(r2));

    long a = to_long(r2 - 1);
    long f = CongruenceZ(n,r,r2,a);

    if(f == 0){
        auto finish = std::chrono::high_resolution_clock::now();
        auto duration = finish - start;
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();

        std::printf("%ld is prime.\n",to_long(n));
        std::printf("Time taken: %lld milliseconds\n", time);

        std::string note = "a = " + std::to_string(a) + "; End: a = " + std::to_string(f) + "; r = " + std::to_string(to_long(r)) + "; " + r2name + "(r) = " + std::to_string(to_long(r2));
        fileWrite(n,ncores,true, time, note);

        return true;
    }
    else {
        auto finish = std::chrono::high_resolution_clock::now();
        auto duration = finish - start;
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();

        std::printf("%ld is not prime.\n",to_long(n));
        std::printf("The a which fails is %ld\n",f);
        std::printf("Time taken: %lld milliseconds\n", time);

        std::string note = "a = " + std::to_string(a) + "; End: a = " + std::to_string(f) + "; r = " + std::to_string(to_long(r)) + "; " + r2name + "(r) = " + std::to_string(to_long(r2));
        fileWrite(n, ncores, false, time, note);

        return false;
    }

}

int main (int argc, char * argv[]) {
    // usage: <program> [--phi | --lambda] [n]
    //   n is read from standard input when it is not given as an argument.
    ZZ n;
    bool haveN = false;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--lambda") == 0) {
            useLambda = true;
        }
        else if (std::strcmp(argv[i], "--phi") == 0) {
            useLambda = false;
        }
        else {
            n = conv<ZZ>(argv[i]);
            haveN = true;
        }
    }

    perflog << "Int, Cores, Prime (T/F), Time (milliseconds), Comments\n";

    if (!haveN) {
        std::printf("Enter a positive integer number n you want to be tested:\n");
        std::cin >> n;
    }

    return Lenstra(n) ? 0 : 1;
}
