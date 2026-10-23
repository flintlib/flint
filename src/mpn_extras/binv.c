/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/* Up to this many limbs the Newton steps are done with full low products
   (flint_mpn_mullow_n on zero-padded operands) rather than with middle
   products; the redundant low limbs are cheaper than the dispatch. */
#ifndef FLINT_MPN_BINV_MULLOW_CUTOFF
#define FLINT_MPN_BINV_MULLOW_CUTOFF 10
#endif

/*
    res = x^(-1) mod B^n, for odd x. Port of radix_invmod_bn with B = 2^64.

    Newton iteration: if y = x^(-1) mod B^m then, writing x y = 1 + B^m h,

        y' = y (2 - x y) = y - B^m (y h)   (mod B^(2m)),

    so the new limbs [m, n) of the inverse are -(y h) mod B^(n-m), where
    h = (x y)[m, n) is a high product whose low limbs are known to be
    00...01. For large sizes h comes from a guarded middle product
    (_flint_mpn_mulhigh_known_low); in the basecase range full low products
    are used.
*/

/* inverse modulo B^2 */
static void
_binv_2(mp_ptr res, mp_limb_t x0, mp_limb_t x1)
{
    mp_limb_t r0, h, t;

    r0 = n_binvert(x0);
    /* (x r0)[1] = mulhi(x0, r0) + x1 r0 */
    umul_ppmm(h, t, x0, r0);
    FLINT_ASSERT(t == 1);
    h += x1 * r0;
    res[0] = r0;
    res[1] = -(r0 * h);
}

/* inverse modulo B^3 from an inverse modulo B^2 */
static void
_binv_2_to_3(mp_ptr res, mp_limb_t x0, mp_limb_t x1, mp_limb_t x2)
{
    mp_limb_t r0 = res[0], r1 = res[1];
    mp_limb_t h, hi, lo, cy, t1;

    /* (x r)[2] where x r == 1 (mod B^2): limb 2 of
       x0 r0 + B (x0 r1 + x1 r0) + B^2 (x0 r2 ... ) */
    umul_ppmm(hi, lo, x0, r0);          /* lo == 1 */
    (void) lo;
    t1 = hi;
    umul_ppmm(hi, lo, x0, r1);
    add_ssaaaa(cy, t1, 0, t1, 0, lo);
    h = hi + cy;
    umul_ppmm(hi, lo, x1, r0);
    add_ssaaaa(cy, t1, 0, t1, 0, lo);
    h += hi + cy;
    FLINT_ASSERT(t1 == 0);
    h += x1 * r1 + x2 * r0;
    res[2] = -(r0 * h);
}

void
flint_mpn_binv(mp_ptr res, mp_srcptr x, mp_size_t xn, mp_size_t n)
{
    mp_size_t a[FLINT_BITS];
    mp_size_t i, m;
    mp_limb_t x1, x2;
    mp_ptr u, scratch;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(xn >= 1);
    FLINT_ASSERT(x[0] & 1);

    if (n == 1)
    {
        res[0] = n_binvert(x[0]);
        return;
    }

    x1 = (xn > 1) ? x[1] : 0;
    x2 = (xn > 2) ? x[2] : 0;

    if (n <= 3)
    {
        _binv_2(res, x[0], x1);
        if (n == 3)
            _binv_2_to_3(res, x[0], x1, x2);
        return;
    }

    /* schedule n -> ceil(n/2) -> ... -> 3 */
    a[i = 0] = n;
    while (a[i] > 3)
    {
        a[i + 1] = (a[i] + 1) / 2;
        i++;
    }

    _binv_2(res, x[0], x1);
    _binv_2_to_3(res, x[0], x1, x2);

    TMP_START;
    u = TMP_ALLOC(2 * n * sizeof(mp_limb_t));
    scratch = u + n;

    for (i--; i >= 0; i--)
    {
        mp_size_t rxn;

        m = a[i + 1];           /* current precision */
        n = a[i];               /* new precision, m < n <= 2m */

        /* the product res * x has m + xn limbs; limbs >= n are not needed */
        rxn = FLINT_MIN(n, m + xn);

        if (n <= FLINT_MPN_BINV_MULLOW_CUTOFF)
        {
            /* u = x * res mod B^n with both operands padded to n limbs */
            mp_limb_t xp[FLINT_MPN_BINV_MULLOW_CUTOFF];
            mp_srcptr xx = x;

            if (xn < n)
            {
                flint_mpn_copyi(xp, x, xn);
                flint_mpn_zero(xp + xn, n - xn);
                xx = xp;
            }

            flint_mpn_zero(res + m, n - m);
            flint_mpn_mullow_n(u, xx, res, n);
            FLINT_ASSERT(u[0] == 1 && flint_mpn_zero_p(u + 1, m - 1));

            /* res[m, n) = -(res * u[m, n)) mod B^(n-m); only the low
               n - m <= m limbs of res enter */
            flint_mpn_mullow_n(scratch, res, u + m, n - m);
            mpn_neg(res + m, scratch, n - m);
        }
        else
        {
            mp_limb_t one = 1;

            /* u = (res * x)[m, rxn) */
            _flint_mpn_mulhigh_known_low(u, res, m, x, FLINT_MIN(n, xn),
                &one, 1, m, rxn, scratch);
            /* res[m, n) = -(u * res) mod B^(n-m) */
            flint_mpn_mulmid(res + m, u, rxn - m, res, n - m, 0, n - m);
            mpn_neg(res + m, res + m, n - m);
        }
    }

    TMP_END;
}
