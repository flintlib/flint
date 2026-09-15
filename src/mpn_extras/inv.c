/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fixed.h"

#ifndef FLINT_MPN_INV_NEWTON_CUTOFF
#define FLINT_MPN_INV_NEWTON_CUTOFF 256
#endif

/*
    Correctly truncated reciprocal: q = floor(B^n / x) for (x, xn) with
    x[xn-1] != 0 and n >= xn, written to (q, n - xn + 2) (the top limb is
    nonzero only when x is a power of B).

    Viewing x as a fixed-point number alpha in [1/B, 1) with xn fraction
    limbs, B^n / x = (1/alpha) B^(n-xn). fixed_inv_newton with
    p = n - xn + 3 fraction limbs has error at most 4 B^(-p) / alpha <=
    4 B^(-p+1), i.e. 4 B^-2 at the integer scale, so the integer part is
    certified when the first fraction limb lies in [2, B-2]; otherwise the
    quotient is corrected against the explicit numerator B^n.
*/
void
flint_mpn_inv(mp_ptr q, mp_srcptr x, mp_size_t xn, mp_size_t n)
{
    mp_size_t qn = n - xn + 2, p;
    mp_ptr U, qq;
    TMP_INIT;

    FLINT_ASSERT(xn >= 1);
    FLINT_ASSERT(n >= xn);
    FLINT_ASSERT(x[xn - 1] != 0);

    TMP_START;

    if (FLINT_MIN(xn, qn) < FLINT_MPN_INV_NEWTON_CUTOFF)
    {
        mp_ptr N, r;
        N = TMP_ALLOC((n + 1 + xn) * sizeof(mp_limb_t));
        r = N + n + 1;
        flint_mpn_zero(N, n);
        N[n] = 1;
        mpn_tdiv_qr(q, r, 0, N, n + 1, x, xn);
        TMP_END;
        return;
    }

    p = n - xn + 3;
    U = TMP_ALLOC((p + 3) * sizeof(mp_limb_t));
    U[p + 2] = 0;

    fixed_inv_newton(U, x, xn, p);

    /* integer limb 0 sits at U[p - (n - xn)] = U[3] */
    qq = U + 3;

    if (qq[-1] > 1 && qq[-1] < UWORD_MAX - 1)
    {
        flint_mpn_copyi(q, qq, qn);
    }
    else
    {
        /* verify against B^n: (qq, qn + 1) is within a few units of the
           quotient of the (n + 1)-limb numerator by x */
        mp_ptr N, r;
        N = TMP_ALLOC((2 * (n + 2)) * sizeof(mp_limb_t));
        r = N + n + 2;
        flint_mpn_zero(N, n);
        N[n] = 1;

        if (qq[qn] != 0)
            mpn_sub_1(qq, qq, qn + 1, 1);

        if (qn >= xn)
            flint_mpn_mul(r, qq, qn, x, xn);
        else
            flint_mpn_mul(r, x, xn, qq, qn);

        /* r = qq x has qn + xn = n + 2 limbs */
        while (r[n + 1] != 0 || mpn_cmp(r, N, n + 1) > 0)
        {
            mpn_sub(r, r, n + 2, x, xn);
            mpn_sub_1(qq, qq, qn, 1);
        }

        mpn_sub_n(r, N, r, n + 1);
        r[n + 1] = 0;

        while (!flint_mpn_zero_p(r + xn, n + 2 - xn) || mpn_cmp(r, x, xn) >= 0)
        {
            mpn_add_1(qq, qq, qn, 1);
            mpn_sub(r, r, n + 2, x, xn);
        }

        flint_mpn_copyi(q, qq, qn);
    }

    TMP_END;
}
