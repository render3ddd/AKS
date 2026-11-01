// type test for chrono library steady_clock
// cp test/chronoTest.cpp test/HPCchronoTest.cpp
// $ g++ -g -O2 -std=c++11 -pthread -march=native test/chronoTest.cpp -o test/chronoTest.out -lntl -lgmp -lm

#include <typeinfo>

#include <math.h> // standard libraries
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
// #include <windows.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
// #include <mmsystem.h>
#include <time.h>
#include <ctime>
#include <chrono>
#include <string>

#include "NTL/ZZ.h" // NTL Libraries
#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
#include "NTL/ZZX.h"
#include "NTL/vec_ZZ.h"
#include "NTL/ZZXFactoring.h"
NTL_CLIENT

#include "../devVR/biSegMultiplyZZX.h" // personal header (ZZX getMaxCoeff/evaluate/polyMultiply)

template <typename T> std::string type_name();

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

int main() {
    // std::cout << getTime() << "\n";
    // std::cout << getDate();

    // int a = 850, n = 520, b = 961;
    // printf("\na^n (mod b) = %d^%d (mod %d)\n\n",a,n,b);

    ZZX f;

    cin >> f;
    // test cases
    // [2 10 14 6] = 2 + 10x + 14x^2 +6x^3
    // 
    // 

    auto now = std::chrono::high_resolution_clock::now();
    // auto now = std::chrono::steady_clock::now();
    printf("\nStart\n");

    // // *** DO STUFF HERE ***

    // unsigned int microsecond = 1000000;
    // usleep(3 * microsecond);//sleeps for 3 second

    // for (int i = 10; --i; ) {
    // for (int i = 0; i < 10; ++i) {
    // for (int i = 0; i < 10; --i) {
    // for (int i = 10-1; i >= 0; --i) {
    //     printf("%d\n",i);
    // }

    // long tntl = PowerMod(a,n,b);
    // long test = powMod(a,n,b);

    // ZZ x = conv<ZZ>("12345678901234567890000");
    // std::cout << "x = " << x << "\n";
    // ZZ y = conv<ZZ>("98765432109876543210");
    // std::cout << "y = " << y << "\n\n";

    // ZZ q = x/y;
    // ZZ q = to_ZZ(x/y);

    // printf("Bitwise AND: %d", n & 1);

    ZZ coeffMaxF = getMaxCoeff(f);
    std::cout << "Max coeff of f(x) = " << coeffMaxF << "\n";

    ZZ x = ZZ(6);
    ZZ resEval = evaluate(f,x);
    std::cout << "When x = " << x << ", f(x) = " << resEval << "\n";

    ZZX g;
    cin >> g;

    ZZX fgBin = polyMultiply(f,g);
    ZZX fgNTL = f * g;

    std::cout << "Binary Segmentation: f(x) * g(x) = " << fgBin << "\n";
    std::cout << "NTL Polynomial Prod: f(x) * g(x) = " << fgNTL << "\n";


    auto then = std::chrono::high_resolution_clock::now();
    auto duration = then - now;

    // printf("PowerMod = %ld\n", tntl);
    // printf("powMod = %ld\n", test);

    // printf("\nInteger (ZZ) division: %ld\n", to_long(q));

    std::cout << "\nTime: " << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() << " milliseconds\n\n";
}
