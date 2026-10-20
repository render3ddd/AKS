
#include <vector>
#include <algorithm>

#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"
NTL_CLIENT

bool isPrime(long n) {

    long i = 2;
    while (i <= n/2) {

        if (n % i == 0) {
            return false;
        }

        i += 1;
        i = NextPrime(i);
    }

    return true;
}

void primeFactors(long n, vector<long>& p, vector<long>& e) {

    // std::cout << "\nChecking if " << n << " is prime...\n";

    if (isPrime(n)) {
        p.push_back(n);
        e.push_back(1);
        // std::cout << n << " is PRIME.\n";
        return;
    }

    // std::cout << n << " is NOT PRIME.\nCalculating prime factors and exponents...\n\n";

    long idx = 0;
    long i = 2, end = n/2;
    while (i <= end) {
        // std::cout << "i = " << i << "\n";
        // std::cout << "n = " << n << "\n";

        p.push_back(i);
        e.push_back(0);

        while (n % i == 0) {
            n /= i;
            e[idx] += 1;
        }

        idx +=1;

        i += 1;
        i = NextPrime(i);
    }

}

unsigned long long LCM(vector<long>& n) {

    unsigned long long lcm = 1;
    std::vector<long> p;
    std::vector<long> e;

    std::vector<long> b;
    std::vector<long> t;

    b.resize(n.size());
    t.resize(n.size());

    for (long i : n) {
        primeFactors(i,p,e);

        // std::cout << "Prime factors:\t";
        // for (int i = 0; i < p.size(); i++)
        //     std::cout << p[i] << "\t";

        // std::cout << "\nExponents:\t";
        // for (int i = 0; i < p.size(); i++)
        //     std::cout << e[i] << "\t";
        // std::cout << "\n";

        for (long j = 0; j < p.size(); ++j) {
            if (t[j] < e[j]) {
                t[j] = e[j];
            }
            if (b[j] < p[j]) {
                b[j] = p[j];
            }
        }

        p.clear();
        e.clear();
    }

    // std::cout << "\nLCM Prime factors:\t";
    // for (int i = 0; i < b.size(); i++)
    //     std::cout << b[i] << "\t";

    // std::cout << "\nLCM Exponents:\t\t";
    // for (int i = 0; i < t.size(); i++)
    //     std::cout << t[i] << "\t";
    // std::cout << "\n";

    for (long k = 0; k < t.size(); ++k) {
        lcm *= power_long(b[k],t[k]);
    }

    return lcm;
}

unsigned long long Phi(unsigned long long r){
    unsigned long long eu = 1;

    for (unsigned long long p = 2; p * p <= r; p += 2) {
        if (r % p == 0) {
            eu *= p - 1;
            r /= p;

            while (r % p == 0) {
                eu *= p;
                r /= p;
            }
        }

        if(p == 2) {
            --p;
        }
    }

    unsigned long long eu1 = eu;

    // now r is prime or 1
    if (r == 1) {
        return eu1;
    }
    else {
        return eu1 * (r - 1) ;
    }
}

long subLambda(const long& a, const long& b) {

    unsigned long long lambda;

    if (b <= 2 || a >= 3) {
        lambda = Phi(power_long(a,b));
        // std::cout << "\nsubLambda = phi(a^b) = " << a << "^" << b << " = " << lambda << "\n";
    }
    else if (b >= 3 && a == 2) {
        lambda = Phi(power_long(a,b))/2;
        // std::cout << "\nFinal: lambda = phi(a^b)/2 = " << a << "^" << b << "/2 = " << lambda << "\n";
    }
    else {
        std::printf("\nError: Out of condition bounds.\n");
    }

    return lambda;
}

// long mainLambda(const long& n) {

//     unsigned long long lambda;

//     std::vector<long> p;
//     std::vector<long> e;

//     primeFactors(n,p,e);

//     std::vector<long> comp;

//     for (long i = 0; i < p.size(); ++i) {
//         comp.push_back(subLambda(p[i],e[i]));
//     }

//     lambda = LCM(comp);

//     return lambda;
// }

// // making use of Fermats Theorem
// // a^lambda(r)=1 mod r for gcd(a,r) = 1;

// ZZ CarmichaelEgan(const long& r) {
//     unsigned long long a = 2;

//     if (isPrime(r)) {
//         return ZZ(Phi(r));
//     }

//     // find an a coprime to r
//     while (GCD(ZZ(r),ZZ(a)) != 1) {
//         // std::cout << "\na = " << a << "\n";
//         // std::cout << "GCD(" << r << ", " << a << ") = " << GCD(ZZ(r),ZZ(a)) <<"\n";
//         a += 1;
//     }
//     // std::cout << "\nGCD(" << r << ", " << a << ") = 1\n";
//     // std::cout << "\na = " << a << "\n";

//     // a <= r by construction
//     unsigned long long b = 1; // b is the order of a mod r.
//     ZZ g1; // g1 is a^b
//     // ZZ g2; // g2 is a^b mod r

//     for (;;) {
//         std::cout << "\nb = " << b << "\n";

//         g1 = PowerMod(a,b,r);
//         // g1 = power_ZZ(a,b);
//         // std::cout << "a^b = " << a << "^" << b << " = " << g1 << "\n";
//         // g2 = to_ZZ_p(g1);
//         // g2 = g1 % r;
//         std::cout << "a^b mod r = " << a << "^" << b << " mod " << r << " = " << g1 << "\n";

//         if (g1 == 1) { // if a^b (mod r) = 1
//             break;
//         }
//         else { // else increase value of b.
//             b += 1;
//         }
//     }

//     std::cout << "\nFinal: b = " << b << "\n";

//     long long lambda_r = mainLambda(r);

//     // return ZZ(b);
//     return ZZ(lambda_r);
// }

ZZ Carmichael(const long& n) {

    if (isPrime(n)) {
        return ZZ(Phi(n));
    }

    unsigned long long lambda;

    std::vector<long> p;
    std::vector<long> e;

    primeFactors(n,p,e);

    std::vector<long> comp;

    for (long i = 0; i < p.size(); ++i) {
        comp.push_back(subLambda(p[i],e[i]));
    }

    lambda = LCM(comp);

    return ZZ(lambda);
}