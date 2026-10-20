/*
    A C++ implementation of Lenstra's algorithm, adapted from that from Hua Li:
        https://researchportal.bath.ac.uk/en/publications/the-analysis-and-implementation-of-the-aks-algorithm-and-its-impr
    Compile with:
        $ g++ -g -O2 -std=c++11 -pthread -march=native devVR/LenstraZnx.cpp -o devVR/LenstraZnx.out -lntl -lgmp -lm
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
#include "FindR.h"
#include "Carmichael.h"
// #include "CongruenceZnx.h"
#include "CongruenceZnxHPC.h"

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
    std::string fldr = "logs/Znx/";
    std::string prfx = "LnstrZnx-";
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
    else if(n == 2 || n == 3){
        std::printf("%d is prime.\n", to_int(n));
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

            std::printf("%ld is not prime.\n",to_long(n));
            std::printf("%ld is a divisor.\n",to_long(R));
            std::printf("Time taken: %lld milliseconds\n", time);

            std::string note = std::to_string(to_long(R)) + " is a divisor";
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

    // ZZ r2 = Euler(to_long(r));
    // std::printf("Phi(%ld) = %ld\n",to_long(r),to_long(r2));
    ZZ r2 = Carmichael(to_long(r));
    std::printf("Lambda(%ld) = %ld\n",to_long(r),to_long(r2));

    long a = to_long(r2 - 1);
    long f = CongruenceZnx(n,r,r2,a);

    if(f == 0){
        auto finish = std::chrono::high_resolution_clock::now();
        auto duration = finish - start;
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();

        std::printf("%ld is prime.\n",to_long(n));
        std::printf("Time taken: %lld milliseconds\n", time);

        // std::string note = "a = " + std::to_string(a) + "; End: a = " + std::to_string(f) + "; r = " + std::to_string(to_long(r)) + "; phi(r) = " + std::to_string(to_long(r2));
        std::string note = "a = " + std::to_string(a) + "; End: a = " + std::to_string(f) + "; r = " + std::to_string(to_long(r)) + "; lambda(r) = " + std::to_string(to_long(r2));
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

    // // ZZ p = conv<ZZ>("11663");
    // // ZZ p = conv<ZZ>("11639");
    // ZZ p = conv<ZZ>("23456611"); // not prime
    // ZZ p = conv<ZZ>("4467165232203241"); // not prime
    // ZZ p = conv<ZZ>("1003026954441971");
    // ZZ p = conv<ZZ>("4467165232203221");
    // ZZ p = conv<ZZ>("689960931088884849033689023336009222695077");

    // prime = Lenstra(p);

    // unsigned long long nos[] = {137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347};

    // unsigned long long nos[] = {9923, 9929, 9931, 9941, 9949, 9967, 9973, 10007, 10009, 10037, 10039, 10061, 10067, 10069, 10079, 10091, 10093, 10099, 10103, 10111, 10133};

    // unsigned long long nos[] = {11491, 11497, 11503, 11519, 11527, 11549, 11551, 11579, 11587, 11593, 11597, 11617, 11621, 11633, 11657, 11677, 11681, 11689, 11699, 11701};
    // unsigned long long nos[] = {32353, 32359, 32363, 32369, 32371, 32377, 32381, 32401, 32411, 32413, 32423, 32429, 32441, 32443, 32467, 32479, 32491, 32497, 32503, 32507, 32531, 32533, 32537, 32561};
    // unsigned long long nos[] = {88721, 88729, 88741, 88747, 88771, 88789, 88793, 88799, 88801, 88807, 88811, 88813, 88817, 88819, 88843, 88853, 88861, 88867, 88873, 88883, 88897, 88903};
    // unsigned long long nos[] = {99191, 99223, 99233, 99241, 99251, 99257, 99259, 99277, 99289, 99317, 99347, 99349, 99367, 99371, 99377, 99391, 99397};

    // unsigned long long nos[] = {104527, 104537, 104543, 104549, 104551, 104561, 104579, 104593, 104597, 104623, 104639, 104651, 104659, 104677, 104681, 104683, 104693, 104701, 104707, 104711, 104717, 104723, 104729};
    // unsigned long long nos[] = {107137, 107171, 107183, 107197, 107201, 107209, 107227, 107243, 107251, 107269, 107273, 107279, 107309, 107323, 107339};
    // unsigned long long nos[] = {160373, 160387, 160397, 160403, 160409, 160423, 160441, 160453, 160481, 160483, 160499, 160507, 160541, 160553, 160579, 160583};
    // unsigned long long nos[] = {267863, 267877, 267887, 267893, 267899, 267901, 267907, 267913, 267929, 267941, 267959, 267961, 268003, 268013, 268043, 268049, 268063, 268069};
    // unsigned long long nos[] = {777419, 777421, 777431, 777433, 777437, 777451, 777463, 777473, 777479, 777541, 777551, 777571, 777583, 777589, 777617, 777619};
    // unsigned long long nos[] = {999907, 999917, 999931, 999953, 999959, 999961, 999979, 999983, 1000003, 1000033, 1000037, 1000039, 1000081, 1000099};

    // unsigned long long nos[] = {2968813, 2968831, 2968841, 2968859, 2968871, 2968877, 2968891, 2968907, 2968913, 2968937, 2968961, 2968967, 2969003, 2969009, 2969011, 2969017, 2969023};
    // unsigned long long nos[] = {4740413, 4740419, 4740469, 4740499, 4740509, 4740511, 4740521, 4740523, 4740559, 4740583, 4740623};
    // unsigned long long nos[] = {7740833, 7740851, 7740871, 7740907, 7740913, 7740983, 7740991, 7741007, 7741009};

    // unsigned long long nos[] = {10012141, 10012157, 10012181, 10012193, 10012213, 10012217, 10012229, 10012253, 10012271, 10012313, 10012333};
    // unsigned long long nos[] = {12434629, 12434647, 12434657, 12434659, 12434687, 12434711, 12434717, 12434729, 12434753, 12434761, 12434771, 12434777, 12434831, 12434839};
    // unsigned long long nos[] = {23456597, 23456603, 23456627, 23456669, 23456681, 23456683, 23456717, 23456723, 23456743, 23456747, 23456749, 23456761, 23456789};
    // unsigned long long nos[] = {26933407, 26933429, 26933449, 26933453, 26933461, 26933471, 26933479, 26933497, 26933503, 26933509, 26933537, 26933551, 26933563, 26933579, 26933581, 26933587, 26933591, 26933603, 26933611};
    // unsigned long long nos[] = {86456873, 86456891, 86456899, 86456911, 86456921, 86456933, 86456941, 86456957, 86456959, 86456969, 86456999, 86457031, 86457071, 86457073, 86457079};

    // unsigned long long nos[] = {100041493, 100041509, 100041527, 100041533, 100041551, 100041553, 100041569, 100041589, 100041649, 100041659, 100041661, 100041679, 100041703};
    // unsigned long long nos[] = {100982747, 100982759, 100982767, 100982809, 100982821, 100982851, 100982867, 100982887, 100982897, 100982909, 100982927, 100982933};
    // unsigned long long nos[] = {103736093, 103736131, 103736153, 103736189, 103736231, 103736249, 103736287, 103736291, 103736293};
    // unsigned long long nos[] = {199999831, 199999841, 199999853, 199999889, 199999901, 199999903, 199999931, 199999949, 199999957, 199999963, 199999991, 200000033};
    // unsigned long long nos[] = {879672929, 879672931, 879672947, 879672953, 879672961, 879672977, 879673049, 879673051, 879673073, 879673087, 879673099, 879673117};
    // unsigned long long nos[] = {999999797, 999999883, 999999893, 999999929, 999999937, 1000000007};

    // unsigned long long nos[] = {1001771887, 1001771893, 1001771933, 1001771957, 1001771963, 1001771977, 1001772011, 1001772019, 1001772041, 1001772047, 1001772089, 1001772091};
    // unsigned long long nos[] = {2492813317, 2492813353, 2492813363, 2492813399, 2492813413, 2492813423, 2492813429, 2492813437, 2492813441, 2492813459, 2492813461, 2492813483, 2492813497}
    // unsigned long long nos[] = {2958346843, 2958346861, 2958346879, 2958346903, 2958346913, 2958346921, 2958346939, 2958346943, 2958346961, 2958346967, 2958347027, 2958347047};
    // unsigned long long nos[] = {9999954821, 9999954851, 9999954857, 9999954889, 9999954911, 9999954919, 9999954929, 9999954941, 9999954947, 9999954967, 9999954997};

    // unsigned long long nos[] = {10015571473, 10015571509, 10015571557, 10015571561, 10015571563, 10015571569, 10015571617, 10015571623, 10015571633, 10015571677};
    // unsigned long long nos[] = {44426254963, 44426254967, 44426254973, 44426255053, 44426255143};

    // unsigned long long nos[] = {100020817121, 100020817159, 100020817181, 100020817201, 100020817219, 100020817277, 100020817327, 100020817331};
    // unsigned long long nos[] = {333267326563, 333267326573, 333267326597, 333267326603, 333267326627, 333267326683, 333267326699, 333267326701, 333267326719, 333267326731, 333267326747, 333267326767};
    // unsigned long long nos[] = {435465768541, 435465768547, 435465768553, 435465768557, 435465768559, 435465768611, 435465768613, 435465768643, 435465768659, 435465768671, 435465768733};

    // unsigned long long nos[] = {1000528294751, 1000528294771, 1000528294793, 1000528294837, 1000528294849, 1000528294861, 1000528294919, 1000528294937, 1000528294943};
    // unsigned long long nos[] = {3222583708391, 3222583708499, 3222583708507, 3222583708567};
    // unsigned long long nos[] = {6999999999827, 6999999999857, 6999999999911, 6999999999931, 6999999999973, 7000000000009};
    // unsigned long long nos[] = {9999999999863, 9999999999971, 10000000000037};

    // unsigned long long nos[] = {10083087720577, 10083087720593, 10083087720661, 10083087720679, 10083087720701, 10083087720709, 10083087720737, 10083087720779};
    // unsigned long long nos[] = {35466059872447, 35466059872613, 35466059872651};
    // unsigned long long nos[] = {99999999999821, 99999999999829, 99999999999853, 99999999999923, 99999999999929, 99999999999931, 99999999999959, 99999999999971, 99999999999973, 100000000000031}

    // unsigned long long nos[] = {112272535095109, 112272535095113, 112272535095131, 112272535095187, 112272535095199, 112272535095233, 112272535095239, 112272535095293};
    // unsigned long long nos[] = {281702565145993, 281702565146017, 281702565146047, 281702565146083, 281702565146089, 281702565146117, 281702565146153, 281702565146161, 281702565146179};

    // unsigned long long nos[] = {1003026954441829, 1003026954441857, 1003026954441859, 1003026954441947, 1003026954441961, 1003026954441971};
    // unsigned long long nos[] = {4467165232203221, 4467165232203239, 4467165232203269, 4467165232203283, 4467165232203337, 4467165232203343, 4467165232203361, 4467165232203367, 4467165232203371, 4467165232203397};
    // unsigned long long nos[] = {9007199254740677, 9007199254740727, 9007199254740761, 9007199254740847, 9007199254740881, 9007199254740997};

    // unsigned long long nos[] = {10022390619214619, 10022390619214661, 10022390619214693, 10022390619214777, 10022390619214801, 10022390619214807};

    // unsigned long long nos[] = {100055128505715869, 100055128505715871, 100055128505715913, 100055128505715941, 100055128505715961, 100055128505715983, 100055128505716009};

    // unsigned long long nos[] = {1083717775299973637, 1083717775299973651, 1083717775299973657, 1083717775299973661, 1083717775299973673, 1083717775299973771};

    // unsigned long long nos[] = {29546363270378696821, 29546363270378696837, 29546363270378696953, 29546363270378696999, 29546363270378697007};

    // unsigned long long nos[] = {326070784035774767779, 326070784035774767789, 326070784035774767903, 326070784035774767933, 326070784035774767971};

    // unsigned long long nos[] = {4973004941902396102547, 4973004941902396102559, 4973004941902396102577, 4973004941902396102607, 4973004941902396102619, 4973004941902396102723, 4973004941902396102727};

    // unsigned long long nos[] = {51005856776120585676917, 51005856776120585676973, 51005856776120585676983, 51005856776120585677021, 51005856776120585677099, 51005856776120585677103};

    // unsigned long long nos[] = {614783152143098270145193, 614783152143098270145199, 614783152143098270145221, 614783152143098270145287, 614783152143098270145299, 614783152143098270145349, 614783152143098270145367, 614783152143098270145397};

    // unsigned long long nos[] = {7666569009190249923345199, 7666569009190249923345239, 7666569009190249923345271, 7666569009190249923345277, 7666569009190249923345281};

    // unsigned long long nos[] = {86084043198752959566539087, 86084043198752959566539147, 86084043198752959566539209};

    // unsigned long long nos[] = {912613825844053990694090939, 912613825844053990694091053, 912613825844053990694091143};

    // unsigned long long nos[] = {1000474617637553175973957531, 1000474617637553175973957657, 1000474617637553175973957663};

    // unsigned long long nos[] = {20476096752860587951845236729, 20476096752860587951845236731, 20476096752860587951845236737, 20476096752860587951845236741, 20476096752860587951845236747, 20476096752860587951845236753, 20476096752860587951845236761, 20476096752860587951845236783, 20476096752860587951845236797, 20476096752860587951845236809, 20476096752860587951845236837, 20476096752860587951845236857, 20476096752860587951845236863, 20476096752860587951845236879, 20476096752860587951845236899, 20476096752860587951845236929};

    // unsigned long long nos[] = {387121083116233373653498534693, 387121083116233373653498534753, 387121083116233373653498534811, 387121083116233373653498534849};

    // unsigned long long nos[] = {4313339400115792413779939217997, 4313339400115792413779939218099};

    // ZZ nos[] = {conv<ZZ>("53437079999999999999999994656083"),conv<ZZ>("53437079999999999999999994656113"),conv<ZZ>("53437079999999999999999994656137"),conv<ZZ>("53437079999999999999999994656153"),conv<ZZ>("53437079999999999999999994656293")};



    ZZ nos[] = {conv<ZZ>("689960931088884849033689023336009222694927"),conv<ZZ>("689960931088884849033689023336009222694971"),conv<ZZ>("689960931088884849033689023336009222695053"),conv<ZZ>("689960931088884849033689023336009222695077")};

    int nosSize = sizeof(nos)/sizeof(*nos);
    int nosEnd = (sizeof(nos)/sizeof(*nos)) - 1;

//    for (int i = 0; i < nosSize; ++i) {
    // for (unsigned long long i = nos[0]; i < nos[nosEnd] + 1; ++i) {
    // for (ZZ i = ZZ(nos[0]); i < ZZ(nos[nosEnd] + 1); ++i) {
    // // for (int i = 5; i < 506; ++i) {
        ZZ n;
        // n = 0;

        std::printf("Enter a positive integer number n you want to be tested:\n");
        std::cin >> n;
        Lenstra(n);

        // prime = Lenstra(ZZ(nos[i]));
//        Lenstra(nos[i]);
        // prime = Lenstra(ZZ(i));
//    }

}