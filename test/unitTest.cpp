// type test for chrono library high_resolution_clock
// cp test/chronoTest.cpp test/HPCchronoTest.cpp
// $ g++ -g -O2 -std=c++11 -pthread -march=native test/unitTest.cpp -o test/unitTest.out -lntl -lgmp -lm
// $ g++ -g -O2 -std=c++11 -pthread -march=native test/unitTest.cpp -o test/unitTest.out -lntl -lgmp -lm && ./test/unitTest.out

#include <typeinfo>
#include <algorithm>

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

// #include "../devVR/biSegMultiplyZZX.h" // personal headers
#include "../devVR/biSegMultiplyZZpX.h"
// #include "../archive/biSegMultiplyZZpXpreJHD.h"
// #include "../test/perflogBiSegMultiplyZZpX.h"

#include "../devVR/Carmichael.h"
// #include "../BathAKS/Carmichael.h"
#include "../devVR/Euler.h"

std::string compfile = "logs/other/BSegSquareVsNTL.csv";
std::ofstream complog(compfile, std::ios::app); // output result into file

// ZZ evaluate(const ZZ_pX& f, const ZZ& x) {
//     /*
//         Returns the integer result when the polynomial f(y) (mod p)
//         is evaluated for y = x in O(D) time where D is the number of
//         terms in f(x) (mod p), including those with coefficients
//         equal to zero.
//     */

//     long fDeg = deg(f);
//     ZZ ans = rep(ConstTerm(f));

//     for (long i = 1; i <= fDeg; ++i) {
//         ans += rep(coeff(f,i))*power(x,i);
//     }

//     return ans;
// }

// ZZ fastEval(const ZZ_pX& f, const ZZ& x) {
//     /*
//         Returns the integer result when the polynomial f(y) (mod p)
//         is evaluated for y = x in O(D) time where D is the number of
//         terms in f(x) (mod p), including those with coefficients
//         equal to zero.
//     */

//     long fDeg = deg(f);
//     ZZ ans = rep(ConstTerm(f));

//     ZZ xpow = ZZ(1);
//     for (long i = 1; i <= fDeg; ++i) {
//         xpow = xpow*x;
//         ans += rep(coeff(f,i))*xpow;
//     }

//     return ans;
// }

int main(int argc, char * argv[]) {

    // test cases
    // [2 10 14 6] = 2 + 10x + 14x^2 + 6x^3
    // [9 16 7 1]  = 9 + 16x + 07x^2 + 1x^3
    // [14 10 6 2]
    // [16 9 7 1]
    // [2 -1 1]
    // [2 -1 -1]

    printf("\nStart\n");

    // // *** DO STUFF HERE ***

    // ############### Multiply Tests ############################################

    // ZZX f;
    // std::cout << "Enter a polynomial using its coefficients in the place of the numbers";
    // std::cout << " in the list [0,1,2,...,n] where n is its exponent (e.g. [2 -1 1] ";
    // std::cout << " = x^2 - x + 2):\n";
    // cin >> f;

    // ZZX g;
    // std::cout << "Enter the polynomial you wish to multiply it by:\n";
    // cin >> g;

    // ZZX h;
    // std::cout << "And the polynomial you wish to mod the product by:\n";
    // cin >> h;

    // double then = GetTime();

    // ZZX fgNTL = f * g;
    // double lap1 = GetTime();

    // ZZX fgBin = polyMultiply(f,g);
    // double lap2 = GetTime();

    // double tNTL = lap1-then;
    // std::cout << "NTL Polynomial Prod: f(x) * g(x) = " << fgNTL << " (" << tNTL*1000 << " milliseconds)\n";

    // double tBin = lap2-lap1;
    // std::cout << "Binary Segmentation: f(x) * g(x) = " << fgBin << " (" << tBin*1000 << " milliseconds)\n";

    // if (fgNTL == fgBin) {
    //     printf("\nBinary Segmentation Multiply Successful.\n");
    // }
    // else {
    //     printf("\nBinary Segmentation Multiply Failed.\n");
    // }

    // ############### MulMod Tests ############################################

    // ZZX f;
    // std::cout << "Enter a polynomial using its coefficients in the place of the numbers";
    // std::cout << " in the list [0,1,2,...,n] where n is its exponent (e.g. [2 -1 1] ";
    // std::cout << " = x^2 - x + 2):\n";
    // cin >> f;

    // ZZX g;
    // std::cout << "Enter the polynomial you wish to multiply it by:\n";
    // cin >> g;

    // ZZX h;
    // std::cout << "And the polynomial you wish to mod the product by:\n";
    // cin >> h;

    // double then = GetTime();

    // ZZX fgnNTL = MulMod(f,g,h);
    // double lap1 = GetTime();

    // ZZX fgnBin = polyMulMod(f,g,h);
    // double lap2 = GetTime();

    // double tNTL = lap1-then;
    // std::cout << "NTL Polynomial Prod: f(x) * g(x) (mod h(x)) = " << fgnNTL << " (" << tNTL*1000 << " milliseconds)\n";

    // double tBin = lap2-lap1;
    // std::cout << "Binary Segmentation: f(x) * g(x) (mod h(x)) = " << fgnBin << " (" << tBin*1000 << " milliseconds)\n";

    // if (fgnNTL == fgnBin) {
    //     printf("\nBinary Segmentation MulMod Successful.\n");
    // }
    // else {
    //     printf("\nBinary Segmentation MulMod Failed.\n");
    // }

    // ############### PowMod Tests ############################################

    // ZZX f;
    // std::cout << "Enter a polynomial using its coefficients in the place of the numbers";
    // std::cout << " in the list [0,1,2,...,n] where n is its exponent (e.g. [2 -1 1] ";
    // std::cout << " = x^2 - x + 2):\n";
    // cin >> f;

    // ZZ e;
    // std::cout << "Enter the exponent you wish to raise it to:\n";
    // cin >> e;

    // ZZX g;
    // std::cout << "And the polynomial you wish to mod the product by:\n";
    // cin >> g;

    // double then = GetTime();
    // ZZX fegBin = polyPowMod(f,e,g);
    // double lap1 = GetTime();

    // // ZZ_pX f_p = to_ZZ_pX(f);
    // // ZZ_pX g_p = to_ZZ_pX(g);

    // double lap2 = GetTime();
    // // ZZ_pX fegNTL = PowerMod(f_p,e,g_p);
    // ZZX fegNTL;
    // SetCoeff(fegNTL,0,1);
    // while (e > 0) {
    //     if (e % 2 == 1) {               // if (n is odd) then
    //         fegNTL = fegNTL*f % g;
    //     }
    //     f = f*f % g;    // a = a^2 (mod b)
    //     e /= 2;                         // n = n/2
    // }
    // double lap3 = GetTime();

    // double tNTL = lap3-lap2;
    // std::cout << "NTL Polynomial PowerMod: f(x)^e (mod h(x)) = " << fegNTL << " (" << tNTL*1000 << " milliseconds)\n";

    // double tBin = lap1-then;
    // std::cout << "Binary Sgmtn polyPowMod: f(x)^e (mod h(x)) = " << fegBin << " (" << tBin*1000 << " milliseconds)\n";

    // if (fegNTL == fegBin) {
    //     printf("\nBinary Segmentation PowMod Successful.\n");
    // }
    // else {
    //     printf("\nBinary Segmentation PowMod Failed.\n");
    // }

    // ############### PowMod ZZ_p Tests ############################################

    // ZZ_p::init(ZZ(7));

    // ZZ mdls = ZZ_p::modulus();
    // std::cout << "\nModulus is: " << mdls << "\n\n";

    // ZZ_pX f;
    // std::cout << "Enter a polynomial using its coefficients in the place of the numbers";
    // std::cout << " in the list [0,1,2,...,n] where n is its exponent (e.g. [2 -1 1] ";
    // std::cout << " = x^2 - x + 2):\n";
    // cin >> f;
    // std::cout << "\nf(x) = " << f << "\n\n";

    // ZZ e;
    // std::cout << "Enter the exponent you wish to raise it to:\n";
    // cin >> e;
    // std::cout << "\ne = " << e << "\n\n";

    // ZZ_pX g;
    // std::cout << "And the polynomial you wish to mod the product by:\n";
    // cin >> g;
    // std::cout << "\ng(x) = " << g << "\n\n";

    // std::cout << "\nOperation " << f << "^" << e << " mod " << g << "\n";

    // ZZ_pX fegBin = ZZpPowMod(f,e,g);
    // ZZ_pX fegNTL = PowerMod(f,e,g);

    // auto start = std::chrono::high_resolution_clock::now();
    // for (int i = 0; i < 100000; ++i) {
    //     PowerMod(f,e,g);
    // }
    // auto lap1 = std::chrono::high_resolution_clock::now();
    // std::cout << "\nNTL: f(x)^e (mod p, g(x)) = ";
    // std::cout << "PowerMod(" << f << ", " << e << ", " << g << ") = " << fegNTL << "\n";

    // auto duration1 = lap1 - start;
    // auto tNTL = std::chrono::duration_cast<std::chrono::milliseconds>(duration1).count();
    // std::printf("Time taken: %ld milliseconds\n",tNTL);

    // auto lap2 = std::chrono::high_resolution_clock::now();
    // for (int i = 0; i < 100000; ++i) {
    //     ZZpPowMod(f,e,g);
    // }
    // auto lap3 = std::chrono::high_resolution_clock::now();
    // std::cout << "\nBin: f(x)^e (mod p, g(x)) = ";
    // std::cout << "PowerMod(" << f << ", " << e << ", " << g << ") = " << fegBin << "\n";

    // auto duration2 = lap3 - lap2;
    // auto tBin = std::chrono::duration_cast<std::chrono::milliseconds>(duration2).count();
    // std::printf("Time taken: %ld milliseconds\n",tBin);

    // long double tNTLd = tNTL;
    // long double tBind = tBin;

    // long double diff = ((tBind - tNTLd)/tNTLd)*100;

    // if (fegNTL == fegBin) {
    //     printf("\nBinary Segmentation PowMod Successful.\n");
    //     if (diff < 0) {
    //         std::printf("New function is %LG percent faster.\n",diff);
    //     }
    //     else {
    //         std::printf("New function is %LG percent slower.\n",diff);
    //     }
    // }
    // else {
    //     printf("\nBinary Segmentation PowMod Failed.\n");
    // }

    // ##################### TEST CASES ##############################################
    // [2 10 14 6] = 2 + 10x + 14x^2 + 6x^3
    // [9 16 7 1]  = 9 + 16x + 07x^2 + 1x^3
    // [14 10 6 2]
    // [16 9 7 1]
    // [2 -1 1]
    // [2 -1 -1]

    // Fix: [1 0 0 1] ^ 8 mod [-1 1 1 0 -1 0 1]
    // Fix: [1 1 1 1 1] ^ 2 mod [-1 1 1 0 -1 0 1]
    // Fix: [1 1 1 1 1] ^ 25 mod [-1 1 1 0 -1 0 1]
    // [133 197 27 69 13 53 17 34 407 190 82 39 48 3]
    // [133 197 27 69 13 53 17 34 407 190 82 39 48 14 10 6 2]
    // 4885403941208064

    // ############### PowModNpX Tests ############################################

    // ZZ n = conv<ZZ>(2);
    // ZZ_p::init(n);

    // ZZX f;
    // std::cout << "Enter a polynomial using its coefficients in the place of the numbers";
    // std::cout << " in the list [0,1,2,...,n] where n is its exponent (e.g. [2 -1 1] ";
    // std::cout << " = x^2 - x + 2):\n";
    // cin >> f;

    // ZZ_pX f_p;
    // long fDeg = deg(f);
    // for (long i = 0; i <= fDeg; ++i) {
    //     SetCoeff(f_p,i,to_ZZ_p(coeff(f,i)));
    // }
    // std::cout << "\nf(x) = " << f << "\n";
    // std::cout << "f(x) = " << f_p << "\n\n";

    // ZZ e;
    // std::cout << "Enter the exponent you wish to raise it to:\n";
    // cin >> e;
    // std::cout << "\ne = " << e << "\n\n";

    // ZZX g;
    // std::cout << "And the polynomial you wish to mod the product by:\n";
    // cin >> g;

    // ZZ_pX g_p;
    // long gDeg = deg(g);
    // for (long i = 0; i <= gDeg; ++i) {
    //     SetCoeff(g_p,i,to_ZZ_p(coeff(g,i)));
    // }
    // std::cout << "\ng(x) = " << g << "\n";
    // std::cout << "g(x) = " << g_p << "\n\n";

    // std::cout << "\nOperation " << f << "^" << e << " mod " << g << "\n";

    // double then = GetTime();
    // ZZX fegBin = PowModNpX(f,e,g,n);
    // double lap1 = GetTime();

    // double lap2 = GetTime();
    // ZZ_pX fegNTL = PowerMod(f_p,e,g_p);
    // double lap3 = GetTime();

    // double tNTL = lap3-lap2;
    // std::cout << "NTL Polynomial PowerMod: f(x)^e (mod h(x)) = " << fegNTL << " (" << tNTL*1000 << " milliseconds)\n";

    // double tBin = lap1-then;
    // std::cout << "Binary Sgmtn polyPowMod: f(x)^e (mod h(x)) = " << fegBin << " (" << tBin*1000 << " milliseconds)\n";

    // if (fegNTL == fegBin) {
    //     printf("\nBinary Segmentation PowMod Successful.\n");
    // }
    // else {
    //     printf("\nBinary Segmentation PowMod Failed.\n");
    // }

    // ############### fastEval Tests ############################################

    // ZZ_p::init(ZZ(11));

    // ZZ mdls = ZZ_p::modulus();

    // std::cout << "\nModulus is: " << mdls << "\n\n";

    // ZZ_pX f;
    // std::cout << "Enter a polynomial using its coefficients in the place of the numbers";
    // std::cout << " in the list [0,1,2,...,n] where n is its exponent (e.g. [2 -1 1] ";
    // std::cout << " = x^2 - x + 2):\n";
    // cin >> f;
    // std::cout << "\nf(x) = " << f << "\n\n";

    // ZZ x;
    // std::cout << "Enter the value you wish to solve for:\n";
    // cin >> x;
    // std::cout << "\ne = " << x << "\n\n";

    // auto start = std::chrono::high_resolution_clock::now();
    // for (int i = 0; i < 10000000; ++i) {
    //     slowEval(f,x);
    // }
    // std::cout << "Original: ";
    // std::cout << "f(" << x << ") = " << slowEval(f,x) << "\n";
    // auto lap1 = std::chrono::high_resolution_clock::now();

    // auto duration1 = lap1 - start;
    // auto eval = std::chrono::duration_cast<std::chrono::milliseconds>(duration1).count();
    // std::printf("Time taken: %ld milliseconds\n",eval);

    // auto lap2 = std::chrono::high_resolution_clock::now();
    // for (int i = 0; i < 10000000; ++i) {
    //     fastEval(f,x);
    // }
    // auto lap3 = std::chrono::high_resolution_clock::now();
    // std::cout << "Faster: ";
    // std::cout << "f(" << x << ") = " << fastEval(f,x) << "\n";

    // auto duration2 = lap3 - lap2;
    // auto fast = std::chrono::duration_cast<std::chrono::milliseconds>(duration2).count();
    // std::printf("Time taken: %ld milliseconds\n",fast);

    // long double evald = eval;
    // long double fastd = fast;

    // long double diff = ((fastd - evald)/evald)*100;

    // if (fastEval(f,x) == evaluate(f,x)) {
    //     std::printf("Faster evaluation successful.\n");
    //     if (diff < 0) {
    //         std::printf("New function is %LG percent faster.\n",diff);
    //     }
    //     else {
    //         std::printf("New function is %LG percent slower.\n",diff);
    //     }
    // }
    // else {
    //     std::printf("Faster evaluation failed.\n");
    // }

    // ##################### TEST CASES ##############################################
    // [2 10 14 6] = 2 + 10x + 14x^2 + 6x^3
    // [9 16 7 1]  = 9 + 16x + 07x^2 + 1x^3
    // [14 10 6 2]
    // [16 9 7 1]
    // [2 -1 1]
    // [2 -1 -1]

    // Fix: [1 0 0 1] ^ 8 mod [-1 1 1 0 -1 0 1]
    // Fix: [1 1 1 1 1] ^ 2 mod [-1 1 1 0 -1 0 1]
    // [133 197 27 69 13 53 17 34 407 190 82 39 48 3]
    // 4885403941208064

    // // ############### Carmichael Tests ############################################

    // // unsigned long long test = 241;
    // // unsigned long long ans = 240;

    // // unsigned long long test = 256;
    // // unsigned long long ans = 64;

    // // unsigned long long test = 360;
    // // unsigned long long ans = 12;

    // // unsigned long long test = 1190;
    // // unsigned long long ans = 48;

    // // unsigned long long test = 1312;
    // // unsigned long long ans = 40;         // ans/5

    // // unsigned long long test = 1314;
    // // unsigned long long ans = 72;

    // // unsigned long long test = 1331;
    // // unsigned long long ans = 1210;

    // // unsigned long long test = 1348;
    // // unsigned long long ans = 336;        // ans/2

    // // unsigned long long test = 1356;
    // // unsigned long long ans = 112;

    // // unsigned long long test = 9323;
    // // unsigned long long ans = 9322;

    // // unsigned long long test = 49392;
    // // unsigned long long ans = 588;

    // // unsigned long long test = 12545938;
    // // unsigned long long ans = 72240;

    // // unsigned long long test = 0;
    // // unsigned long long ans = 0;

    // unsigned long long test = 1;
    // unsigned long long ans = 1;

    // auto start = std::chrono::high_resolution_clock::now();
    // ZZ phi = Euler(test);
    // auto lap1 = std::chrono::high_resolution_clock::now();
    // ZZ lmd = Carmichael(test);
    // auto lap2 = std::chrono::high_resolution_clock::now();

    // auto timePhi = std::chrono::duration_cast<std::chrono::milliseconds>(lap1 - start).count();
    // auto timeLmd = std::chrono::duration_cast<std::chrono::milliseconds>(lap2 - lap1).count();

    // std::cout << "Euler Totient Function: phi(" << test << ") = " << phi << "(" << timePhi <<" milliseconds)" << "\n";
    // std::cout << "Carmichael Function: lambda(" << test << ") = " << lmd << "(" << timeLmd <<" milliseconds)" << "\n";

    // if (lmd == ans) {
    //     std::printf("Function computed successfully.\n");
    // }
    // else {
    //     std::cout << "Function failed - lambda(" << test << ") != " << ans << "\n";
    // }

    // // ########################### MAIN LAMBDA HELPERS ##########################################

    // long n = 1348;
    // std::cout << "\n n = "<< n << "\n";
    // std::cout << "\n n/2 = "<< n/2 << "\n";

    // vector<long> p, e;

    // long rootn = sqrt(n);
    // // std::cout << "sqrt(" << n << ") = " << rootn << "\n";
    // // long i = 2;
    // // while (i <= sqrt(2)*sqrt(n)) {
    // //     std::cout << "i = " << i << "\n";
    // //     i += 1;
    // //     i = NextPrime(i);
    // // }

    // // std::cout << isPrime(n) << "\n";

    // primeFactors(n,p,e);

    // // e.erase(std::remove(e.begin(), e.end(), 0), e.end());

    // std::cout << "Prime factors:\t";
    // for (int i = 0; i < p.size(); i++)
    //     std::cout << p[i] << "\t";

    // std::cout << "\n\nExponents:\t";
    // for (int i = 0; i < p.size(); i++)
    //     std::cout << e[i] << "\t";
    // std::cout << "\n";

    // std::vector<long> comp;

    // for (long i = 0; i < p.size(); ++i) {
    //     comp.push_back(subLambda(p[i],e[i]));
    // }

    // long ones = 1;
    // // comp.erase(std::remove(comp.begin(), comp.end(), ones), comp.end());

    // std::cout << "\nLambda Components:\t";
    // for (int i = 0; i < p.size(); i++)
    //     std::cout << comp[i] << "\t";
    // std::cout << "\n\n";

    // long testLCM = LCM(comp);
    // std::cout << "\nLambda(" << n << ") = " << testLCM << "\n";

    // // ############### PowMod ZZ_p Profiling Tests ############################################

    // ZZ_p::init(ZZ(11));

    // ZZ mdls = ZZ_p::modulus();
    // std::cout << "\nModulus is: " << mdls << "\n\n";

    // ZZ_pX f;
    // std::cout << "Enter a polynomial using its coefficients in the place of the numbers";
    // std::cout << " in the list [0,1,2,...,n] where n is its exponent (e.g. [2 -1 1] ";
    // std::cout << " = x^2 - x + 2):\n";
    // cin >> f;
    // std::cout << "\nf(x) = " << f << "\n\n";

    // ZZ e;
    // std::cout << "Enter the exponent you wish to raise it to:\n";
    // cin >> e;
    // std::cout << "\ne = " << e << "\n\n";

    // ZZ_pX g;
    // std::cout << "And the polynomial you wish to mod the product by:\n";
    // cin >> g;
    // std::cout << "\ng(x) = " << g << "\n\n";

    // std::cout << "\nOperation " << f << "^" << e << " mod " << g << "\n";

    // ZZ_pX fegBin = ZZpPowMod(f,e,g);
    // ZZ_pX fegNTL = PowerMod(f,e,g);

    // auto start = std::chrono::high_resolution_clock::now();
    // for (int i = 0; i < 100000; ++i) {
    //     PowerMod(f,e,g);
    // }
    // auto lap1 = std::chrono::high_resolution_clock::now();
    // std::cout << "\nNTL: f(x)^e (mod p, g(x)) = ";
    // std::cout << "PowerMod(" << f << ", " << e << ", " << g << ") = " << fegNTL << "\n";

    // auto duration1 = lap1 - start;
    // auto tNTL = std::chrono::duration_cast<std::chrono::milliseconds>(duration1).count();
    // std::printf("Time taken: %ld milliseconds\n",tNTL);

    // auto lap2 = std::chrono::high_resolution_clock::now();
    // for (int i = 0; i < 100000; ++i) {
    //     ZZpPowMod(f,e,g);
    // }
    // auto lap3 = std::chrono::high_resolution_clock::now();
    // std::cout << "\nBin: f(x)^e (mod p, g(x)) = ";
    // std::cout << "PowerMod(" << f << ", " << e << ", " << g << ") = " << fegBin << "\n";

    // auto duration2 = lap3 - lap2;
    // auto tBin = std::chrono::duration_cast<std::chrono::milliseconds>(duration2).count();
    // std::printf("Time taken: %ld milliseconds\n",tBin);

    // long double tNTLd = tNTL;
    // long double tBind = tBin;

    // long double diff = ((tBind - tNTLd)/tNTLd)*100;

    // if (fegNTL == fegBin) {
    //     printf("\nBinary Segmentation PowMod Successful.\n");
    //     if (diff < 0) {
    //         std::printf("New function is %LG percent faster.\n",diff);
    //     }
    //     else {
    //         std::printf("New function is %LG percent slower.\n",diff);
    //     }
    // }
    // else {
    //     printf("\nBinary Segmentation PowMod Failed.\n");
    // }

    // // ##################### TEST CASES ##############################################
    // // [2 10 14 6] = 2 + 10x + 14x^2 + 6x^3
    // // [9 16 7 1]  = 9 + 16x + 07x^2 + 1x^3
    // // [14 10 6 2]
    // // [16 9 7 1]
    // // [2 -1 1]
    // // [2 -1 -1]

    // // Fix: [1 0 0 1] ^ 8 mod [-1 1 1 0 -1 0 1]
    // // Fix: [1 1 1 1 1] ^ 2 mod [-1 1 1 0 -1 0 1]
    // // Fix: [1 1 1 1 1] ^ 25 mod [-1 1 1 0 -1 0 1]
    // // [133 197 27 69 13 53 17 34 407 190 82 39 48 3]
    // // [133 197 27 69 13 53 17 34 407 190 82 39 48 14 10 6 2]
    // // 4885403941208064

    // ############### PowMod ZZ_p Profiling Tests 2 ############################################

    // ZZ nos[] = {conv<ZZ>("11"),conv<ZZ>("347"),conv<ZZ>("9929"),conv<ZZ>("32353"),conv<ZZ>("104527"),conv<ZZ>("7741009"),conv<ZZ>("2968813"),conv<ZZ>("86456873"),conv<ZZ>("200000033"),conv<ZZ>("9999954997"),conv<ZZ>("10015571473"),conv<ZZ>("435465768733"),conv<ZZ>("1000528294751"),conv<ZZ>("99999999999973"),conv<ZZ>("1003026954441829"),conv<ZZ>("100055128505715869"),conv<ZZ>("29546363270378696821"),conv<ZZ>("4973004941902396102547"),conv<ZZ>("614783152143098270145193"),conv<ZZ>("86084043198752959566539087"),conv<ZZ>("1000474617637553175973957531"),conv<ZZ>("387121083116233373653498534693"),conv<ZZ>("53437079999999999999999994656293"),conv<ZZ>("689960931088884849033689023336009222695077")};

    // int nosSize = sizeof(nos)/sizeof(*nos);
    // int nosEnd = (sizeof(nos)/sizeof(*nos)) - 1;

    // for (int i = 0; i < nosSize; ++i) {
    //     ZZ n = nos[i];

    //     // ZZ mdls = ZZ_p::modulus();
    //     // std::cout << "\nModulus is: " << mdls << "\n\n";
    //     std::cout << "\nn = " << n << "\n";
    //     std::printf("Running %d of %d\n",i+1,nosSize+1);

    //     // Find a suitable r
    //     ZZ r = to_ZZ(2);
    //     ZZ R;
    //     ZZ r1;

    //     while(r < n){
    //         ZZ R = GCD(r, n);
    //         if(R != 1 ){
    //             std::printf("%ld is not prime.\n",to_long(n));
    //             std::printf("%ld is a divisor.\n",to_long(R));
    //             break;
    //         }
    //         else {
    //             ZZ v = to_ZZ(floor(power_long(to_long(log(n)), 2)));

    //             // order of n mod r is bigger than v;
    //             int p = 0;
    //             ZZ_p::init(r); // calculate mod r

    //             while(v <= r){
    //                 ZZ x = to_ZZ(power_long(to_long(n), to_long(v))); // calculates x = n^v
    //                 ZZ_p z = to_ZZ_p(x);
    //                 if(z == to_ZZ_p(1)){
    //                     r1 = r; // store value of r
    //                     r = n + 1;
    //                     break;
    //                 }
    //                 else{
    //                     v = v + 1;
    //                 }
    //             }
    //         }
    //         r = r + 1;
    //     }

    //     r = r1;

    //     ZZ r2 = Euler(to_long(r));
    //     std::printf("Phi(%ld) = %ld\n",to_long(r),to_long(r2));

    //     long a = to_long(r2 - 1);

    //     ZZ_p::init(n);                      // initialise mod n

    //     ZZ_pX b = ZZ_pX(to_long(r), 1) - 1; // b = x^r - 1 (mod n);
    //     ZZ_pX e = ZZ_pX(1, 1);              // e = x (mod n)
    //     ZZ_pX d = PowerMod(e, n, b);        // d = x^n (mod b, n)

    //     std::cout << "a = " << a << "\n";
    //     ZZ_pX c = ZZ_pX(1, 1) - a;    // c = x - a (mod n);
    //     ZZ_pX fBSg = ZZpPowMod(c, n, b);        // f = (x - a)^n (mod b, n) - LHS

    //     if (fBSg == PowerMod(c, n, b)) {
    //         printf("\nBinary Segmentation PowMod Successful.\n");
    //     }
    //     else {
    //         printf("\nBinary Segmentation PowMod Failed.\n");
    //     }
    // }

    // // ##################### TEST CASES ##############################################
    // // [2 10 14 6] = 2 + 10x + 14x^2 + 6x^3
    // // [9 16 7 1]  = 9 + 16x + 07x^2 + 1x^3
    // // [14 10 6 2]
    // // [16 9 7 1]
    // // [2 -1 1]
    // // [2 -1 -1]

    // // Fix: [1 0 0 1] ^ 8 mod [-1 1 1 0 -1 0 1]
    // // Fix: [1 1 1 1 1] ^ 2 mod [-1 1 1 0 -1 0 1]
    // // Fix: [1 1 1 1 1] ^ 25 mod [-1 1 1 0 -1 0 1]
    // // [133 197 27 69 13 53 17 34 407 190 82 39 48 3]
    // // [133 197 27 69 13 53 17 34 407 190 82 39 48 14 10 6 2]
    // // 4885403941208064

    // ############### Multiply ZZ_pX Profiling Tests 3 ############################################

    ZZ nos[] = {conv<ZZ>("11"),conv<ZZ>("347"),conv<ZZ>("9929"),conv<ZZ>("32353"),conv<ZZ>("104527"),conv<ZZ>("7741009"),conv<ZZ>("2968813"),conv<ZZ>("86456873"),conv<ZZ>("200000033"),conv<ZZ>("9999954997"),conv<ZZ>("10015571473"),conv<ZZ>("435465768733"),conv<ZZ>("1000528294751"),conv<ZZ>("99999999999973"),conv<ZZ>("1003026954441829"),conv<ZZ>("100055128505715869"),conv<ZZ>("29546363270378696821"),conv<ZZ>("4973004941902396102547"),conv<ZZ>("614783152143098270145193"),conv<ZZ>("86084043198752959566539087"),conv<ZZ>("1000474617637553175973957531"),conv<ZZ>("387121083116233373653498534693"),conv<ZZ>("53437079999999999999999994656293")};
    // ZZ nos[] = {conv<ZZ>("11"),conv<ZZ>("347")};

    int nosSize = sizeof(nos)/sizeof(*nos);
    int nosEnd = (sizeof(nos)/sizeof(*nos)) - 1;

    complog << "n,r,phi(r),a,tBSg,tNTL,CORRECT\n";

    if (complog.is_open()) {
        printf("Writing to logs...");
        for (int i = 0; i < nosSize; ++i) {
            ZZ n = nos[i];

            // ZZ mdls = ZZ_p::modulus();
            // std::cout << "\nModulus is: " << mdls << "\n\n";
            std::cout << "\nn = " << n << "\n";
            std::printf("Running %d of %d\n",i+1,nosSize+1);

            // Find a suitable r
            ZZ r = to_ZZ(2);
            ZZ R;
            ZZ r1;

            while(r < n){
                ZZ R = GCD(r, n);
                if(R != 1 ){
                    std::printf("%ld is not prime.\n",to_long(n));
                    std::printf("%ld is a divisor.\n",to_long(R));
                    break;
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

            ZZ r2 = Euler(to_long(r));
            std::printf("Phi(%ld) = %ld\n",to_long(r),to_long(r2));

            long a = to_long(r2 - 1);

            ZZ_p::init(n);                      // initialise mod n

            ZZ_pX b = ZZ_pX(to_long(r), 1) - 1; // b = x^r - 1 (mod n);
            ZZ_pX e = ZZ_pX(1, 1);              // e = x (mod n)
            ZZ_pX d = PowerMod(e, n, b);        // d = x^n (mod b, n)

            std::cout << "a = " << a << "\n";
            ZZ_pX c = ZZ_pX(1, 1) - a;    // c = x - a (mod n);

            ZZ_pX pmBsg;
            ZZ_pX pmNTL;

            auto startNTL = std::chrono::high_resolution_clock::now();

            pmNTL = PowerMod(c,n,b);

            auto stopNTL = std::chrono::high_resolution_clock::now();
            auto durationNTL = stopNTL - startNTL;
            auto tNTL = std::chrono::duration_cast<std::chrono::milliseconds>(durationNTL).count();
            std::printf("NTL Time:\t%lld milliseconds\n",static_cast<long long>(tNTL));

            auto startBsg = std::chrono::high_resolution_clock::now();

            pmBsg = ZZpPowMod(c,n,b);

            auto stopBsg = std::chrono::high_resolution_clock::now();
            auto durationBsg = stopBsg - startBsg;
            auto tBsg = std::chrono::duration_cast<std::chrono::milliseconds>(durationBsg).count();
            std::printf("BSg Time:\t%lld milliseconds\n",static_cast<long long>(tBsg));

            bool equal;
            if (pmBsg == pmNTL) {
                equal = true;
                printf("\nBinary Segmentation PowMod Successful.\n");
            }
            else {
                equal = false;
                printf("\nBinary Segmentation PowMod Failed.\n");
            }

            complog << n << "," << r << "," << r2 << "," << a <<  "," << tBsg << "," << tNTL <<  "," << equal << "\n";
        }
    }

    printf("Stop\n");

}
