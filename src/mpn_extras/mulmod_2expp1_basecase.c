/*
    Copyright 2009 Jason Moxham

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpn_extras.h"

/* ret + (xp,n) = (yp,n)*(zp,n) % 2^b+1
   needs (tp,2n) temp space, everything reduced mod 2^b
   inputs, outputs are fully reduced
   NOTE: 2n is not the same as 2b rounded up to nearest limb
*/
static inline int
flint_mpn_mulmod_2expp1_internal(mp_ptr xp, mp_srcptr yp, mp_srcptr zp,
    flint_bitcnt_t b, mp_ptr tp)
{
    mp_size_t n, k;
    mp_limb_t c;

    n = BITS_TO_LIMBS(b);
    k = GMP_NUMB_BITS * n - b;

#if 0
    flint_mpn_mul_large(tp, yp, n, zp, n);
#else
    if (yp == zp)
        flint_mpn_sqr(tp, yp, n);
    else
        flint_mpn_mul_n(tp, yp, zp, n);
#endif

    if (k == 0)
    {
        c = mpn_sub_n(xp, tp, tp + n, n);
        return mpn_add_1 (xp, xp, n, c);
    }

    c = tp[n - 1];
    tp[n - 1] &= GMP_NUMB_MASK >> k;

#if HAVE_NATIVE_mpn_sublsh_nc
    c = mpn_sublsh_nc (xp, tp, tp + n, n, k, c);
#else
    {
        mp_limb_t c1;
        c1 = mpn_lshift (tp + n, tp + n, n, k);
        tp[n] |= c >> (GMP_NUMB_BITS - k);
        c = mpn_sub_n (xp, tp, tp + n, n) + c1;
    }
#endif
    c = mpn_add_1 (xp, xp, n, c);
    xp[n - 1] &= GMP_NUMB_MASK >> k;
    return c;
}

/* c is the top bits of the inputs, must be fully reduced */
int
flint_mpn_mulmod_2expp1_basecase (mp_ptr xp, mp_srcptr yp, mp_srcptr zp, int c,
    flint_bitcnt_t b, mp_ptr tp)
{
    int cy, cz;
    mp_size_t n, k;

    cy = c & 2;
    cz = c & 1;
    n = BITS_TO_LIMBS(b);
    k = GMP_NUMB_BITS * n - b;

    if (cy == 0)
    {
        if (cz == 0)
        {
            c = flint_mpn_mulmod_2expp1_internal(xp, yp, zp, b, tp);
        }
        else
        {
            c = mpn_neg(xp, yp, n);
            c = mpn_add_1 (xp, xp, n, c);
            xp[n - 1] &= GMP_NUMB_MASK >> k;
        }
    }
    else
    {
        if (cz == 0)
	     {
            c = mpn_neg(xp, zp, n);
            c = mpn_add_1(xp, xp, n, c);
            xp[n - 1] &= GMP_NUMB_MASK >> k;
        }
        else
        {
            c = 0;
            xp[0] = 1;
            flint_mpn_zero(xp + 1, n - 1);
        }
    }

    return c;
}

