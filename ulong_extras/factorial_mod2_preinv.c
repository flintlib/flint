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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"


static const mp_limb_t small_factorials[] =
{
    1UL, 1UL, 2UL, 6UL, 24UL, 120UL, 720UL, 5040UL, 40320UL, 362880UL,
    3628800UL, 39916800UL, 479001600UL,
#if FLINT64
    6227020800UL, 87178291200UL,
    1307674368000UL, 20922789888000UL, 355687428096000UL, 6402373705728000UL,
    121645100408832000UL, 2432902008176640000UL,
#endif
};

#if FLINT64
#define MAX_SMALL_FACTORIAL 20
#else
#define MAX_SMALL_FACTORIAL 12
#endif


mp_limb_t n_factorial_mod2_preinv(ulong n, mp_limb_t p, mp_limb_t pinv)
{
    mp_limb_t prod, hi, lo;

    if (n <= MAX_SMALL_FACTORIAL)
        return n_mod2_preinv(small_factorials[n], p, pinv);

    if (n >= p)
        return 0UL;

    if (n >= 1000000UL)
        return n_factorial_fast_mod2_preinv(n, p, pinv);

    prod = small_factorials[MAX_SMALL_FACTORIAL];
    lo = n;
    n--;

    /* TODO: speedup for n in the range of sqrt(ULONG_MAX) */
    while (n > MAX_SMALL_FACTORIAL)
    {
        umul_ppmm(hi, lo, lo, n);

        if (hi)
        {
            lo = n_ll_mod_preinv(hi, lo, p, pinv);
            prod = n_mulmod2_preinv(prod, lo, p, pinv);
            lo = 1UL;
        }

        n--;
    }

    return n_mulmod2_preinv(prod, lo, p, pinv);
}
