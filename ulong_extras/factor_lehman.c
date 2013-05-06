/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 William Hart
    
******************************************************************************/

#undef ulong /* prevent clash with standard library */
#include <math.h>
#define ulong unsigned long
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t n_factor_lehman(mp_limb_t n)
{
    double limit;
    mp_limb_t cuberoot, k;
    n_factor_t factors;
    long bound;

#if FLINT64 /* cannot compute enough primes */
    if (n > 10000000000000000UL) return n;
#endif

    if ((n & 1) == 0) return 2;

    limit = pow(n, 1.0/3.0);
    cuberoot = (mp_limb_t) ceil(limit);
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
        mp_limb_t x = (mp_limb_t) ceil(low - 0.0001);
        mp_limb_t end = (mp_limb_t) floor(0.0001 + low + pow(n, 1.0/6.0)/((double) 4.0*sqrt((double) k)));
        mp_limb_t sub = k*n*4;

        for ( ; x <= end; x++)
        {
            mp_limb_t p, sq = x*x - sub;
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
