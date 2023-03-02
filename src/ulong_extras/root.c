/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2015 Fredrik Johansson
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#define ulong ulongxx /* interferes with system includes */
#include <math.h>
#undef ulong
#define ulong mp_limb_t
#include "flint.h"
#include "ulong_extras.h"
    
/* A table of precomputed inverses of values from 1 to 64  
   inv_table[n] = i/n for all n in range[1, 64]  
   inv_table[0] is set to 0, albeit it will never be called */

static const double inv_table[] = { 
                       0.000000000000000, 1.000000000000000, 0.500000000000000, 
                       0.333333333333333, 0.250000000000000, 0.200000000000000, 
                       0.166666666666667, 0.142857142857143, 0.125000000000000, 
                       0.111111111111111, 0.100000000000000, 0.090909090909091, 
                       0.083333333333333, 0.076923076923077, 0.071428571428571, 
                       0.066666666666667, 0.062500000000000, 0.058823529411765, 
                       0.055555555555556, 0.052631578947368, 0.050000000000000, 
                       0.047619047619048, 0.045454545454545, 0.043478260869565, 
                       0.041666666666667, 0.040000000000000, 0.038461538461538, 
                       0.037037037037037, 0.035714285714286, 0.034482758620690, 
                       0.033333333333333, 0.032258064516129, 0.031250000000000, 
                       0.030303030303030, 0.029411764705882, 0.028571428571429, 
                       0.027777777777778, 0.027027027027027, 0.026315789473684, 
                       0.025641025641026, 0.025000000000000, 0.024390243902439, 
                       0.023809523809524, 0.023255813953488, 0.022727272727273, 
                       0.022222222222222, 0.021739130434783, 0.021276595744681, 
                       0.020833333333333, 0.020408163265306, 0.020000000000000, 
                       0.019607843137255, 0.019230769230769, 0.018867924528302, 
                       0.018518518518519, 0.018181818181818, 0.017857142857143, 
                       0.017543859649123, 0.017241379310345, 0.016949152542373, 
                       0.016666666666667, 0.016393442622951, 0.016129032258065, 
                       0.015873015873016, 0.015625000000000 };

                       /* this table consists 65 values */

/* This table has the max possible base for a given root. 
   max_base[n] = UWORD_MAX^(1/n) for n in range [1, FLINT_BITS] 
   max_base[0] is set to 0, although it will never be called */

static const mp_limb_t max_base[] = { 
#ifdef FLINT64
                        UWORD(0), UWORD_MAX, UWORD(4294967296), UWORD(2642245), UWORD(65536),
                        UWORD(7131), UWORD(1625), UWORD(565), UWORD(256), UWORD(138), UWORD(84),
                        UWORD(56), UWORD(40), UWORD(30), UWORD(23), UWORD(19), UWORD(16),
                        UWORD(13), UWORD(11), UWORD(10), UWORD(9), UWORD(8), UWORD(7),
                        UWORD(6), UWORD(6), UWORD(5), UWORD(5), UWORD(5), UWORD(4), UWORD(4),
                        UWORD(4), UWORD(4), UWORD(4), UWORD(3), UWORD(3), UWORD(3), UWORD(3),
                        UWORD(3), UWORD(3), UWORD(3), UWORD(3), UWORD(2), UWORD(2), UWORD(2),
                        UWORD(2), UWORD(2), UWORD(2), UWORD(2), UWORD(2), UWORD(2), UWORD(2),
                        UWORD(2), UWORD(2), UWORD(2), UWORD(2), UWORD(2), UWORD(2), UWORD(2),
                        UWORD(2), UWORD(2), UWORD(2), UWORD(2), UWORD(2), UWORD(2), UWORD(2)
#else
                        UWORD(0), UWORD(4294967295),  UWORD(65535), UWORD(1625), UWORD(255), 
                        UWORD(84),  UWORD(40), UWORD(23), UWORD(15), UWORD(11), UWORD(9),
                        UWORD(7), UWORD(6), UWORD(5), UWORD(4), UWORD(4), UWORD(3), UWORD(3),
                        UWORD(3), UWORD(3), UWORD(3), UWORD(2), UWORD(2), UWORD(2), UWORD(2),
                        UWORD(2), UWORD(2), UWORD(2), UWORD(2), UWORD(2), UWORD(2), UWORD(2)
#endif
                                                        };

                        /* this table consists of 65 values in case of FLINT64,
                           otherwise 33 */

mp_limb_t
n_root(mp_limb_t n, mp_limb_t root)
{
    mp_limb_t x, currval, base, upper_limit; 
    double dx;

    if (!n || !root)
        return 0;

    if (root == 1)
        return n;

    if (root == 2)
        return n_sqrt(n);

    if (root == 3)
        return n_cbrt(n);

    if (root >= FLINT_BITS || (UWORD(1) << root) > n)
        return 1;

    upper_limit = max_base[root];      /* n <= upper_limit^root */

    x = n_root_estimate((double)n, root);

    /* one round of newton iteration */

    currval = n_pow(x, root-1);
    dx = n / currval;
    dx -= x;
    dx *= inv_table[root];
    dx = floor(dx);
    x += dx;

    base = x;

    if (base >= upper_limit)
        base = upper_limit - 1;

    currval = n_pow(base, root);
    if (currval == n)
        goto final;

    while (currval <= n)
    {
        (base) += 1;
        currval = n_pow(base, root);
        if (base == upper_limit)
            break;
    }

    while (currval > n)
    {
        (base) -= 1;
        currval = n_pow(base, root);
    }

    final:
    return base;
}
