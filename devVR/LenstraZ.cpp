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

#include "PerfectPower.h" //Each Independent Test
#include "Euler.h"
// #include "Carmichael.h"
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
            long v = to_long(to_ZZ(floor(power_long(to_long(log(n)), 2))));

            // order of n mod r is bigger than v;
            //
            // REFACTOR 2.4: same acceptance condition as before, but evaluated
            // with word-sized modular arithmetic on n mod r. Previously
            // n^v was computed in a 64-bit long (power_long), which overflows
            // for any realistic n and v, and ZZ_p::init(r) ran for every r.
            // (Whether this condition is the right one is REFACTOR 1.2.)
            const long rl = to_long(r);
            const long nr = rem(n, rl);                 // n mod r
            long x = PowerMod(nr, v, rl);               // n^v mod r

            while (v <= rl) {
                if (x == 1 % rl) {
                    r1 = r;                             // store value of r
                    r = n + 1;
                    break;
                }
                else {
                    x = MulMod(x, nr, rl);              // n^(v+1) mod r
                    v = v + 1;
                }
            }
        }
        r = r + 1;
    }

    r = r1;
    std::printf("r = %ld\n",to_long(r));

    ZZ r2 = Euler(to_long(r));
    std::printf("Phi(%ld) = %ld\n",to_long(r),to_long(r2));
    // ZZ r2 = Carmichael(to_long(r));
    // std::printf("Lambda(%ld) = %ld\n",to_long(r),to_long(r2));

    long a = to_long(r2 - 1);
    long f = CongruenceZ(n,r,r2,a);

    if(f == 0){
        auto finish = std::chrono::high_resolution_clock::now();
        auto duration = finish - start;
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();

        std::printf("%ld is prime.\n",to_long(n));
        std::printf("Time taken: %lld milliseconds\n", time);

        std::string note = "a = " + std::to_string(a) + "; End: a = " + std::to_string(f) + "; r = " + std::to_string(to_long(r)) + "; phi(r) = " + std::to_string(to_long(r2));
        // std::string note = "a = " + std::to_string(a) + "; End: a = " + std::to_string(f) + "; r = " + std::to_string(to_long(r)) + "; lambda(r) = " + std::to_string(to_long(r2));
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

        std::string note = "a = " + std::to_string(a) + "; End: a = " + std::to_string(f) + "; r = " + std::to_string(to_long(r)) + "; phi(r) = " + std::to_string(to_long(r2));
        // std::string note = "a = " + std::to_string(a) + "; End: a = " + std::to_string(f) + "; r = " + std::to_string(to_long(r)) + "; lambda(r) = " + std::to_string(to_long(r2));
        fileWrite(n, ncores, false, time, note);

        return false;
    }

}

int main (int argc, char * argv[]) {

    perflog << "Int, Cores, Prime (T/F), Time (milliseconds), Comments\n";

    bool prime;

    // ZZ p = conv<ZZ>("11663");
    // ZZ p = conv<ZZ>("11639");
    // ZZ p = conv<ZZ>("160427");
    // ZZ p = conv<ZZ>("7740919");
    // ZZ p = conv<ZZ>("23456611");
    ZZ p = conv<ZZ>("1003026954441971");

    prime = Lenstra(p);

    // // int nos[] = {137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347};
    // // int nos[] = {11491, 11497, 11503, 11519, 11527, 11549, 11551, 11579, 11587, 11593, 11597, 11617, 11621, 11633, 11657, 11677, 11681, 11689, 11699, 11701};
    // // ZZ nos[] = {689960931088884849033689023336009222694927, 689960931088884849033689023336009222694971, 689960931088884849033689023336009222695053, 689960931088884849033689023336009222695077};
    // long nos[] = {23456597, 23456603, 23456627, 23456669, 23456681, 23456683, 23456717, 23456723, 23456743, 23456747, 23456749, 23456761, 23456789};
    // // long nos[] = {1003026954441829, 1003026954441857, 1003026954441859, 1003026954441947, 1003026954441961, 1003026954441971};
    // // unsigned long long nos[] = {9007199254740677, 9007199254740727, 9007199254740761, 9007199254740847, 9007199254740881, 9007199254740997};
    // // unsigned long long nos[] = {10012141, 10012157, 10012181, 10012193, 10012213, 10012217, 10012229, 10012253, 10012271, 10012313, 10012333};
    // int nosSize = sizeof(nos)/sizeof(*nos);
    // int nosEnd = (sizeof(nos)/sizeof(*nos)) - 1;

    // for (int i = 0; i < nosSize; ++i) {
    // // for (unsigned long long i = nos[0]; i < nos[nosEnd] + 1; ++i) {
    // // // for (int i = 5; i < 506; ++i) {
    //     // ZZ n;
    //     // n = 0;

    //     // std::printf("Enter a positive integer number n you want to be tested:\n");
    //     // std::cin >> n;

    //     prime = Lenstra(to_ZZ(nos[i]));
    //     // prime = Lenstra(ZZ(i));
    // }

}