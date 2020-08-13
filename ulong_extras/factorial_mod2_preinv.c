/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"


static const mp_limb_t small_factorials[] =
{
    UWORD(1), UWORD(1), UWORD(2), UWORD(6), UWORD(24), UWORD(120), UWORD(720), UWORD(5040), UWORD(40320), UWORD(362880),
    UWORD(3628800), UWORD(39916800), UWORD(479001600),
#if FLINT64
    UWORD(6227020800), UWORD(87178291200),
    UWORD(1307674368000), UWORD(20922789888000), UWORD(355687428096000), UWORD(6402373705728000),
    UWORD(121645100408832000), UWORD(2432902008176640000),
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
        return UWORD(0);

    if (n >= UWORD(1000000))
        return n_factorial_fast_mod2_preinv(n, p, pinv);

    prod = small_factorials[MAX_SMALL_FACTORIAL];
    lo = n;
    n--;

    /* TODO: speedup for n in the range of sqrt(UWORD_MAX) */
    while (n > MAX_SMALL_FACTORIAL)
    {
        umul_ppmm(hi, lo, lo, n);

        if (hi)
        {
            lo = n_ll_mod_preinv(hi, lo, p, pinv);
            prod = n_mulmod2_preinv(prod, lo, p, pinv);
            lo = UWORD(1);
        }

        n--;
    }

    return n_mulmod2_preinv(prod, lo, p, pinv);
}
