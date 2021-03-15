/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "mpn_extras.h"

#if FLINT_BITS == 64

/* Above NEWTON_CUTOFF, we need one Newton iteration. Shifting down
   the numerator and denominator by MAGIC_BITS allows doing the
   division using 32-bit unsigned ints. This constant can be changed
   slightly, but not too much, since we have to avoid overflow in the
   numerator and we simultaneously want to guarantee sufficiently
   many accurate quotient bits. */
#define NEWTON_CUTOFF (2 * 53)
#define NEED_NEWTON_ABOVE_HIGH_LIMB (UWORD(1) << (NEWTON_CUTOFF - 64))
#define MAGIC_BITS 45
#define TWO_EXP_64 18446744073709551616.0
#define MAX_LIMBIFYABLE_DOUBLE 1.844674407370955e+19
#define UU_LE(ahi, alo, bhi, blo) ((ahi) < (bhi) || ((ahi) == (bhi) && (alo <= blo)))

mp_size_t
flint_mpn_sqrtrem2(mp_ptr sp, mp_ptr rp, mp_srcptr np)
{
    double x;
    mp_limb_t s, rhi, rlo, thi, tlo, hi, lo;
    unsigned int rm, sm;

    hi = np[1];
    lo = np[0];

    FLINT_ASSERT(hi != 0);

    if (hi > NEED_NEWTON_ABOVE_HIGH_LIMB)
    {
        x = sqrt(hi * TWO_EXP_64 + lo);
        x = FLINT_MIN(x, MAX_LIMBIFYABLE_DOUBLE);
        s = (mp_limb_t) x;

        umul_ppmm(rhi, rlo, s, s);
        if (UU_LE(rhi, rlo, hi, lo))
        {
            /* r = n - x^2 */
            sub_ddmmss(rhi, rlo, hi, lo, rhi, rlo);
            /* rm = r >> MAGIC_BITS */
            rm = (rhi << (FLINT_BITS - MAGIC_BITS)) | (rlo >> MAGIC_BITS);
            sm = s >> (MAGIC_BITS - 1);
            s += rm / sm;
            s -= (s == 0);  /* fix possible wraparound */
        }
        else
        {
            /* r = x^2 - n */
            sub_ddmmss(rhi, rlo, rhi, rlo, hi, lo);
            /* as above, but round up to get a better approximation of the integer square root */
            rm = (rhi << (FLINT_BITS - MAGIC_BITS)) | (rlo >> MAGIC_BITS);
            sm = s >> (MAGIC_BITS - 1);
            s -= (rm + sm - 1) / sm;
        }
    }
    else
    {
        s = (mp_limb_t) sqrt(hi * TWO_EXP_64 + lo);
    }

    /* The integer square root satisfies s^2 <= n < (s+1)^2,
       or 0 <= r <= 2*s  where r = n - s^2. */

    /* If s^2 > n, we need to decrement s. */
    umul_ppmm(rhi, rlo, s, s);
    while (!UU_LE(rhi, rlo, hi, lo))
    {
        /* s^2 > n; set s' = s - 1 and s'^2 = s^2 - (2s-1)  */
        sub_ddmmss(thi, tlo, (s >> (FLINT_BITS - 1)), (s << 1), 0, 1);
        sub_ddmmss(rhi, rlo, rhi, rlo, thi, tlo);
        s -= 1;
    }

    /* rhi, rlo = n - s^2 */
    sub_ddmmss(rhi, rlo, hi, lo, rhi, rlo);
    /* thi, tlo = 2s */
    thi = s >> (FLINT_BITS - 1);
    tlo = s << 1;
    while (!UU_LE(rhi, rlo, thi, tlo))
    {
        /* r > 2*s; set s' = s + 1 and r' = r - (2s+1)  */
        add_ssaaaa(thi, tlo, thi, tlo, 0, 1);
        sub_ddmmss(rhi, rlo, rhi, rlo, thi, tlo);
        s += 1;
        thi = s >> (FLINT_BITS - 1);
        tlo = s << 1;
    }

    sp[0] = s;

    if (rp != NULL)
    {
        rp[0] = rlo;
        rp[1] = rhi;
    }

    if (rhi != 0)
        return 2;
    return (rlo != 0);
}

#else

/* todo: could implement the algorithm of n_sqrtrem if we cared
   about optimizing for 32-bit machines */
mp_size_t
flint_mpn_sqrtrem2(mp_ptr sp, mp_ptr rp, mp_srcptr np)
{
    return mpn_sqrtrem(sp, rp, np, 2);
}

#endif

