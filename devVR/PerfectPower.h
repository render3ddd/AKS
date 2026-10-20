#pragma once

#include <NTL/RR.h> // NTL Libraries
#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"
#include "NTL/ZZX.h"

// function to calculate if n = a^b
// takes input ZZ n and returns 1 if n is a perfect power, 0 otherwise

inline int PerfectPower(const NTL::ZZ& n){
    long b = 2;
    NTL::RR k = NTL::to_RR(log(n) / log(2));
    NTL::ZZ a;

    while(b <= to_long(k)){
        // b cannot be bigger than k
        long c = NTL::to_long(ceil((log(n) / log(2)) / NTL::to_long(b)));
        a = pow(2, c); // assign guess value for a

        while(power(a, b) > n){
            double d = to_double(((b - 1) * a + n / power(a, (b - 1))) / b); // Apply Integer Newton's Method
            NTL::ZZ e = NTL::to_ZZ(floor(d));
            a = to_long(e); // adjust a
        }

        if(n == power(a, b)){
            // if n is a perfect power
            printf("n = a^b\n");
            printf("b = %ld\n",b);
            printf("a = %ld\n\n",to_long(a));
            return 1;
        }
        else{
            b = b + 1;
        }
    }

    if(n != power(a,b)){
        return 0; // n is not a perfect power
    }

    return 0;
}