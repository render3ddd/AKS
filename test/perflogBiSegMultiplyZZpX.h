// NOTE (REFACTOR 1.10): instrumented copy used to produce the Chapter 5 profiling
// figures. It still contains the ORIGINAL, unfixed binary segmentation (2^b - 1
// evaluation point, digit compared after reduction mod p). Use
// devVR/biSegMultiplyZZpX.h for anything other than reproducing those figures.


#include <math.h>
#include <vector>
#include <algorithm>

#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
NTL_CLIENT

// std::string filename = "logs/other/BSegVsNTL.csv";
std::string filename = "logs/other/BSegProfile2.csv";
std::ofstream perflog(filename, std::ios::app); // output result into file

ZZ_p getMaxCoeff(const ZZ_pX& f) {
    /*
        Returns the largest integer coefficient, f_i, of the input
        polynomial, f(x) (mod p), in O(D) time where D is the number
        of terms in f(x) (mod p), including those with coefficients
        equal to zero.
    */

    ZZ_p f_i = ConstTerm(f);
    ZZ comp = rep(ConstTerm(f));

    for (long i = 1; i <= deg(f); ++i) {
        if (comp < rep(coeff(f,i))) {
            f_i = coeff(f,i);
            comp = rep(coeff(f,i));
        }
        else {
            continue;
        }
    }

    return f_i;
}

ZZ evaluate(const ZZ_pX& f, const ZZ& x) {
    /*
        Returns the integer result when the polynomial f(y) (mod p)
        is evaluated for y = x in O(D) time where D is the number of
        terms in f(x) (mod p), including those with coefficients
        equal to zero.
    */

    long fDeg = deg(f);
    ZZ ans = rep(ConstTerm(f));

    ZZ xpow = ZZ(1);
    for (long i = 1; i <= fDeg; ++i) {
        xpow = xpow*x;
        ans += rep(coeff(f,i))*xpow;
    }

    return ans;
}

ZZ_pX ZZpXmultiply(const ZZ_pX& f, const ZZ_pX& g) {
    /*
        Returns the polynomial product s(x) =  f(x) * g(x) (mod p),
        using binary segmentation, in O(D) time where D is the number
        of terms in s(x), including those with coefficients equal to zero.
    */

    ZZ mdls = ZZ_p::modulus();

    std::ofstream perflog;
    perflog.open(filename, std::ios::app);

    if (perflog.is_open()) {
        // ################ INITIALISE ########################################

        auto start1 = std::chrono::high_resolution_clock::now();

        // ZZ termsF = to_ZZ(deg(f)+1);
        // ZZ termsG = to_ZZ(deg(g)+1);
        long termsF = deg(f)+1;
        long termsG = deg(g)+1;

        ZZ maxTermFG = to_ZZ(std::max(termsF,termsG));
        ZZ_p maxCoeffF = getMaxCoeff(f);
        ZZ_p maxCoeffG = getMaxCoeff(g);
        // ZZ maxCoeff = ZZ_p::modulus();

        ZZ rhs = maxTermFG * rep(maxCoeffF) * rep(maxCoeffG);   // rhs = max(U,V) * max(f_j) * max(g_k)
        // ZZ rhs = maxTermFG * sqr(maxCoeff);                     // rhs = max(U,V) * max(f_j) * max(g_k)

        long b = 1;
        ZZ lhs = (power2_ZZ(b)) - 1;                            // lhs = 2^b-1
        while (lhs <= rhs) {                                    // Choose b such that 2^b − 1 > max(U,V) * max(f_j) * max(g_k)
            b = b + 1;
            lhs = (power2_ZZ(b)) - 1;
        }

        for (int i = 0; i < 4; ++i) {
            // ZZ termsF = to_ZZ(deg(f)+1);
            // ZZ termsG = to_ZZ(deg(g)+1);
            long termsF = deg(f)+1;
            long termsG = deg(g)+1;

            ZZ maxTermFG = to_ZZ(std::max(termsF,termsG));
            ZZ_p maxCoeffF = getMaxCoeff(f);
            ZZ_p maxCoeffG = getMaxCoeff(g);
            // ZZ maxCoeff = ZZ_p::modulus();

            ZZ rhs = maxTermFG * rep(maxCoeffF) * rep(maxCoeffG);   // rhs = max(U,V) * max(f_j) * max(g_k)
            // ZZ rhs = maxTermFG * sqr(maxCoeff);                     // rhs = max(U,V) * max(f_j) * max(g_k)

            long b = 1;
            ZZ lhs = (power2_ZZ(b)) - 1;                            // lhs = 2^b-1
            while (lhs <= rhs) {                                    // Choose b such that 2^b − 1 > max(U,V) * max(f_j) * max(g_k)
                b = b + 1;
                lhs = (power2_ZZ(b)) - 1;
            }
        }

        auto stop1 = std::chrono::high_resolution_clock::now();

        auto duration1 = stop1 - start1;
        auto time1 = std::chrono::duration_cast<std::chrono::milliseconds>(duration1).count();
        // std::printf("Initialise Time:\t%ld\tmilliseconds\n",time1);

        // ################ EVALUATE ########################################

        auto start2 = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < 4; ++i) {

            // ZZ X = lhs;
            ZZ F = evaluate(f,lhs);                                   // evaluate f(x) for x = 2^b - 1
            ZZ G = evaluate(g,lhs);                                   // evaluate g(x) for x = 2^b - 1

            // ZZ m = F * G;                                           // Integer multiply

        }

        ZZ X = lhs;
        ZZ F = evaluate(f,X);                                   // evaluate f(x) for x = 2^b - 1
        ZZ G = evaluate(g,X);                                   // evaluate g(x) for x = 2^b - 1

        // ZZ m = F * G;                                           // Integer multiply

        auto stop2 = std::chrono::high_resolution_clock::now();

        auto duration2 = stop2 - start2;
        auto time2 = std::chrono::duration_cast<std::chrono::milliseconds>(duration2).count();
        // std::printf("Evaluation Time:\t%ld\tmilliseconds\n",time2);

        // ############## INTEGER MULTIPLY ####################################

        auto start5 = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < 4; ++i) {
            ZZ m = F * G;                                           // Integer multiply
        }

        ZZ m = F * G;                                           // Integer multiply

        auto stop5 = std::chrono::high_resolution_clock::now();

        auto duration5 = stop5 - start5;
        auto time5 = std::chrono::duration_cast<std::chrono::milliseconds>(duration5).count();
        // std::printf("Multiplication Time:\t%ld\tmilliseconds\n",time5);

        // ################ CONSTRUCT ########################################

        auto start3 = std::chrono::high_resolution_clock::now();

        // long fgDeg = deg(f) + deg(g);                           // Degree of polynomial product
        // long fgTrm = deg(f) + deg(g) + 1;                       // Terms in polynomial product
        long fgTrm = termsF + termsG - 1;                       // Terms in polynomial product
        ZZ_p s[fgTrm];                                          // Store coefficients and constant in an array
        // vector<ZZ_p> s;                                          // Store coefficients and constant in a vector
        // s.reserve(fgTrm);
        // std::cout << "Terms in s = " << s.capacity() << " = " << fgTrm << "?\n";

        s[0] = to_ZZ_p(m % lhs);                                // Reassemble coefficients into signal
        if (rep(s[0]) > lhs/2) {                                // Extract next b bits: s_i = floor( m/(2b−1)^i ) mod 2^b − 1
            s[0] = s[0] - to_ZZ_p(lhs);
            m = m + lhs;
        }
        for (long i = 1; i < fgTrm; ++i) {
            m = m / lhs;                                        // N.B. -- NTL "/" operator floors result by default
            s[i] = to_ZZ_p(m % lhs);

            if (rep(s[i]) > lhs/2) {
                s[i] = s[i] - to_ZZ_p(lhs);
                m = m + lhs;
            }
        }

        // s.push_back(to_ZZ_p(m % lhs));                                // Reassemble coefficients into signal
        // if (rep(s[0]) > lhs/2) {                                // Extract next b bits: s_i = floor( m/(2b−1)^i ) mod 2^b − 1
        //     s[0] = s[0] - to_ZZ_p(lhs);
        //     m = m + lhs;
        // }
        // for (long i = 1; i < fgTrm; ++i) {
        //     m = m / lhs;                                        // N.B. -- NTL "/" operator floors result by default
        //     s.push_back(to_ZZ_p(m % lhs));

        //     if (rep(s[i]) > lhs/2) {
        //         s[i] = s[i] - to_ZZ_p(lhs);
        //         m = m + lhs;
        //     }
        // }

        for (int i = 0; i < 4; ++i) {
            m = F * G;                                              // Integer multiply

            long fgTrm = termsF + termsG - 1;                       // Terms in polynomial product
            ZZ_p s[fgTrm];                                          // Store coefficients and constant in an array

            s[0] = to_ZZ_p(m % lhs);                                // Reassemble coefficients into signal
            if (rep(s[0]) > lhs/2) {                                // Extract next b bits: s_i = floor( m/(2b−1)^i ) mod 2^b − 1
                s[0] = s[0] - to_ZZ_p(lhs);
                m = m + lhs;
            }
            for (long i = 1; i < fgTrm; ++i) {
                m = m / lhs;                                        // N.B. -- NTL "/" operator floors result by default
                s[i] = to_ZZ_p(m % lhs);

                if (rep(s[i]) > lhs/2) {
                    s[i] = s[i] - to_ZZ_p(lhs);
                    m = m + lhs;
                }
            }

            // m = F * G;                                              // Integer multiply

            // long fgTrm = termsF + termsG - 1;                       // Terms in polynomial product
            // // ZZ_p s[fgTrm];                                          // Store coefficients and constant in an array
            // vector<ZZ_p> s;                                          // Store coefficients and constant in a vector
            // s.reserve(fgTrm);

            // s.push_back(to_ZZ_p(m % lhs));                                // Reassemble coefficients into signal
            // if (rep(s[0]) > lhs/2) {                                // Extract next b bits: s_i = floor( m/(2b−1)^i ) mod 2^b − 1
            //     s[0] = s[0] - to_ZZ_p(lhs);
            //     m = m + lhs;
            // }
            // for (long i = 1; i < fgTrm; ++i) {
            //     m = m / lhs;                                        // N.B. -- NTL "/" operator floors result by default
            //     s.push_back(to_ZZ_p(m % lhs));

            //     if (rep(s[i]) > lhs/2) {
            //         s[i] = s[i] - to_ZZ_p(lhs);
            //         m = m + lhs;
            //     }
            // }
        }

        auto stop3 = std::chrono::high_resolution_clock::now();

        auto duration3 = stop3 - start3;
        auto time3 = std::chrono::duration_cast<std::chrono::milliseconds>(duration3).count();
        // std::printf("Construction Time:\t%ld\tmilliseconds\n",time3);

        // ################ RETURN ########################################

        auto start4 = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < 4; ++i) {
            ZZ_pX polyModProd;                                      // Base-b digits of m are desired coefficients

            polyModProd.SetLength(fgTrm);                           // set the length of the underlying coefficient vector to number of Terms in polynomial product

            for (long j = 0; j < fgTrm; ++j) {
                SetCoeff(polyModProd,j,s[j]);
            }

            polyModProd.normalize();                                // remove leading zeros on coefficients

            // ZZ_pX polyModProd;                                      // Base-b digits of m are desired coefficients

            // polyModProd.SetLength(fgTrm);                           // set the length of the underlying coefficient vector to number of Terms in polynomial product

            // for (long j = 0; j < fgTrm; ++j) {
            //     SetCoeff(polyModProd,j,s.at(j));
            // }

            // polyModProd.normalize();                                // remove leading zeros on coefficients
        }

        ZZ_pX polyModProd;                                      // Base-b digits of m are desired coefficients

        polyModProd.SetLength(fgTrm);                           // set the length of the underlying coefficient vector to number of Terms in polynomial product

        for (long j = 0; j < fgTrm; ++j) {
            // SetCoeff(polyModProd,j,s.at(j));
            SetCoeff(polyModProd,j,s[j]);
        }

        polyModProd.normalize();                                // remove leading zeros on coefficients

        auto stop4 = std::chrono::high_resolution_clock::now();

        auto duration4 = stop4 - start4;
        auto time4 = std::chrono::duration_cast<std::chrono::milliseconds>(duration4).count();
        // std::printf("Return Time:\t\t%ld\tmilliseconds\n",time4);

        auto durationT = stop4 - start1;
        auto timeT = std::chrono::duration_cast<std::chrono::milliseconds>(durationT).count();
        // std::printf("Total Time:\t\t%ld\tmilliseconds\n",timeT);

        ZZ_pX fg;
        NTL::mul(fg,f,g);

        bool correct;
        if (polyModProd == fg) {
            correct = true;
        }
        else {
            correct = false;
        }

        perflog << mdls << "," << deg(f) << "," << deg(g) << "," << max(deg(f),deg(g)) << "," << deg(polyModProd) <<  "," << time1 << "," << time2  << "," << time5 << "," << time3 << "," << time4 << "," << timeT <<  "," << correct << "\n";
        return polyModProd;
    }
    else {
        printf("\n###### FILE FAILED TO OPEN ######\n");
        ZZ_pX fg;
        NTL::mul(fg,f,g);
        return fg;
    }
}

ZZ_pX ZZpPowMod(ZZ_pX a, ZZ n, const ZZ_pX& b) {
    /*
        Calculates a(x)^n (mod b(x), p) in O(log n) time.
    */

    ZZ_pX ans;                          // Initialise answer
    SetCoeff(ans,0,1);

    while (n > 0) {
        if (n % 2 == 1) {               // if (n is odd) then
            // std::printf("\nMultiply:\n");
            ans = ZZpXmultiply(ans,a);
            ans %= b;
        }
        // std::printf("\nSquaring:\n");
        a = ZZpXmultiply(a,a);          // a = a^2 (mod b)
        // a = sqr(a);                     // a = a^2 (mod b)
        a %= b;

        n /= 2;                         // n = n/2
    }

    return ans;
}
