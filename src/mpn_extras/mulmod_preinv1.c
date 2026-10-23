/*
    Copyright (C) 2012 William Hart
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

void flint_mpn_mulmod_preinv1(mp_ptr r,
        mp_srcptr a, mp_srcptr b, mp_size_t n,
        mp_srcptr d, mp_limb_t dinv, ulong norm)
{
    mp_limb_t ts[150];
    mp_ptr t, q;

    FLINT_ASSERT(n > 0);

    if (n <= 30)
        t = ts;
    else
        t = flint_malloc(5 * n * sizeof(mp_limb_t));

    q = t + 2 * n;

    if (a == b)
        flint_mpn_sqr(t, a, n);
    else
        flint_mpn_mul_n(t, a, b, n);

    if (norm)
        mpn_rshift(t, t, 2 * n, norm);

    if (n == 1)
    {
        /* the 3/2 inverse of (d, 0) is the 2/1 inverse of d */
        mp_limb_t qq, rr;
        udiv_qrnnd_preinv(qq, rr, t[1], t[0], d[0], dinv);
        (void) qq;
        r[0] = rr;
    }
    else if (n == 2)
    {
        /* t < d^2, so {t + 2, 2} < d and two 3/2 steps reduce t */
        mp_limb_t qq, r1 = t[3], r0 = t[2];
        FLINT_MPN_UDIV_QR_3BY2(qq, r1, r0, r1, r0, t[1], d[1], d[0], dinv);
        FLINT_MPN_UDIV_QR_3BY2(qq, r1, r0, r1, r0, t[0], d[1], d[0], dinv);
        (void) qq;
        r[0] = r0;
        r[1] = r1;
    }
    else
    {
        /* q needs n limbs and the divide and conquer code n more */
        if (n >= 3 && n < FLINT_MPN_DIV_DC_CUTOFF)
            _flint_mpn_divrem_basecase_preinv1(q, t, 2 * n, d, n, dinv);
        else
            _flint_mpn_divrem_preinv1(q, t, 2 * n, d, n, dinv, q + n);
        flint_mpn_copyi(r, t, n);
    }

    if (n > 30)
        flint_free(t);
}
