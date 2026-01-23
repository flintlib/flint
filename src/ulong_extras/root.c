/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2015, 2025 Fredrik Johansson
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include <math.h>
#include "ulong_extras.h"

static const float inv_table[] = {
   0.200000000000000,
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
   0.025641025641026, 0.025000000000000,
};

/* This table has the max possible base for a given root. For n >= 4,
   max_base[n-4] = floor(UWORD_MAX^(1/n)).*/
static const uint16_t max_base[] = {
#if FLINT_BITS == 64
    65535, 7131, 1625, 565, 255, 138, 84, 56, 40, 30, 23, 19, 15, 13, 11, 10,
    9, 8, 7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2
#else
    255, 84, 40, 23, 15, 11, 9, 7, 6, 5, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2
#endif
};

ulong
n_root(ulong n, ulong root)
{
    ulong x, currval, base, upper_limit;

    if (n == 0 || root == 0)
        return 0;

    if (root == 1)
        return n;

    if (root == 2)
        return n_sqrt(n);

    if (root == 3)
        return n_cbrt(n);

    if (root >= FLINT_BITS || n < (UWORD(1) << root))
        return 1;

    /* n <= upper_limit^root */
    upper_limit = max_base[root - 4];

    /* We have already established above that n >= 2^root */
    if (upper_limit == 2)
        return upper_limit;

    /* upper_limit = 2 for root >= 41 */
    FLINT_ASSERT(root <= 40);

    if (root == 4)
        x = sqrt(sqrt(n));
    else
        x = expf(inv_table[root - 5] * logf(n));

    base = x;

    if (base >= upper_limit)
        base = upper_limit - 1;

    currval = n_pow(base, root);
    if (currval == n)
        return base;

    while (currval <= n)
    {
        base++;
        currval = n_pow(base, root);
        if (base == upper_limit)
            break;
    }

    while (currval > n)
    {
        base--;
        currval = n_pow(base, root);
    }

    return base;
}

ulong
n_rootrem(ulong* remainder, ulong n, ulong root)
{
    ulong x, currval, base, upper_limit;

    if (!root)
        return 0;

    if (n == 0 || root == 1)
    {
        *remainder = 0;
        return n;
    }

    if (root == 2)
        return n_sqrtrem(remainder, n);

    if (root == 3)
        return n_cbrtrem(remainder, n);

    if (root >= FLINT_BITS || n < (UWORD(1) << root))
    {
        *remainder = n - 1;
        return 1;
    }

    /* n <= upper_limit^root */
    upper_limit = max_base[root - 4];

    /* We have already established above that n >= 2^root */
    if (upper_limit == 2)
    {
        *remainder = n - (UWORD(1) << root);
        return upper_limit;
    }

    FLINT_ASSERT(root <= 40);

    if (root == 4)
        x = sqrt(sqrt(n));
    else
        x = expf(inv_table[root - 5] * logf(n));

    base = x;

    if (base >= upper_limit)
        base = upper_limit - 1;

    currval = n_pow(base, root);
    if (currval == n)
    {
        *remainder = 0;
        return base;
    }

    while (currval <= n)
    {
        base++;
        currval = n_pow(base, root);
        if (base == upper_limit)
            break;
    }

    while (currval > n)
    {
        base--;
        currval = n_pow(base, root);
    }

    *remainder = base;
    *remainder = n_pow(*remainder, root);
    *remainder = n - *remainder;
    return base;
}
