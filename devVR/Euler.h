// Euler's phi function
//
// REFACTOR 2.5: computed with Carmichael's lambda from one factorisation in
// Totients.h (trial division to sqrt(r)).

#include "Totients.h"

NTL::ZZ Euler(long r){
    return NTL::to_ZZ(PhiAndLambda(r).phi);
}
