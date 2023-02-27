/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2014, 2015 Dana Jacobsen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int n_is_prime(mp_limb_t n)
{
    /* flint's "BPSW" checked against Feitsma and Galway's database [1, 2] 
       up to 2^64 by Dana Jacobsen.
       [1]  http://www.janfeitsma.nl/math/psp2/database
       [2]  http://www.cecm.sfu.ca/Pseudoprimes/index-2-to-64.html
    */

    if (n < 11) {
        if (n == 2 || n == 3 || n == 5 || n == 7)   return 1;
        else                                        return 0;
    }
    if (!(n%2) || !(n%3) || !(n%5) || !(n%7))       return 0;
    if (n <  121) /* 11*11 */                       return 1;
    if (!(n%11) || !(n%13) || !(n%17) || !(n%19) ||
        !(n%23) || !(n%29) || !(n%31) || !(n%37) ||
        !(n%41) || !(n%43) || !(n%47) || !(n%53))   return 0;
    if (n < 3481) /* 59*59 */                       return 1;
    if (n > 1000000 &&
        (!(n% 59) || !(n% 61) || !(n% 67) || !(n% 71) || !(n% 73) ||
         !(n% 79) || !(n% 83) || !(n% 89) || !(n% 97) || !(n%101) ||
         !(n%103) || !(n%107) || !(n%109) || !(n%113) || !(n%127) ||
         !(n%131) || !(n%137) || !(n%139) || !(n%149)))  return 0;

    return n_is_probabprime(n);
}
