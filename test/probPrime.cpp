// Reference oracle for the AKS sweep: prints 1 if n is (probably) prime, else 0.
// NTL's ProbPrime runs Miller-Rabin; with 40 rounds the error probability is
// below 4^-40, and for the small n swept here it is effectively exact.
#include <iostream>
#include "NTL/ZZ.h"

int main() {
    NTL::ZZ n;
    while (std::cin >> n) {
        std::cout << (n > 1 && NTL::ProbPrime(n, 40) ? 1 : 0) << "\n";
    }
    return 0;
}
