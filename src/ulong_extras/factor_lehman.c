/*
    Copyright (C) 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "flint.h"
#include "ulong_extras.h"

ulong n_factor_lehman(ulong n)
{
    double limit;
    ulong cuberoot, k;
    n_factor_t factors;
    slong bound;

#if FLINT64 /* cannot compute enough primes */
    if (n > UWORD(10000000000000000)) return n;
#endif

    if ((n & 1) == 0) return 2;

    limit = pow(n, 1.0/3.0);
    cuberoot = (ulong) ceil(limit);
    bound = n_prime_pi(cuberoot);

    n_factor_init(&factors);
    if (n_factor_trial_range(&factors, n, 0, bound) != n)
        return factors.p[0];

    if ((factors.p[0] = n_factor_one_line(n, FLINT_FACTOR_ONE_LINE_ITERS)))
        if (factors.p[0] != n)
            return factors.p[0];

    for (k = 1; k <= cuberoot + 1; k++)
    {
        double low = 2.0*sqrt((double) k)*sqrt((double) n);
        ulong x = (ulong) ceil(low - 0.0001);
        ulong end = (ulong) floor(0.0001 + low + pow(n, 1.0/6.0)/((double) 4.0*sqrt((double) k)));
        ulong sub = k*n*4;

        for ( ; x <= end; x++)
        {
            ulong p, sq = x*x - sub;
            if (n_is_square(sq))
            {
                sq = sqrt((double) sq);
                p = n_gcd(n, x - sq);
                if (p != 1)
                    return p;
            }
        }
    }

    return n;
}
