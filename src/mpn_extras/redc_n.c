/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/*
    Montgomery reduction: t = X * B^(-n) mod m, canonical in [0, m),
    for odd m of n limbs and 0 <= X < m * B^n (2n limbs).
    minv = -m^(-1) mod B^n (negation of _flint_mpn_binvert's output).
    scratch: 2n limbs. t must not alias X, m or minv; X is preserved.

    With q = X_lo * minv mod B^n we have X + q*m = t * B^n exactly and
    t = X_hi + floor(q*m / B^n) + (X_lo != 0) < 2m, since the low halves
    satisfy X_lo + lo(q*m) = 0 or B^n exactly.

    The floor is obtained from flint_mpn_mulhigh_n, whose result errs
    low by at most n + 2 ulp of the returned guard limb C [1]. The exact
    guard limb is known here: it is the top limb T of lo(q*m) = B^n - X_lo.
    Writing T = C + e - B*carry with 0 <= e <= n + 2, the omitted carry
    into the high limbs occurred iff T < C, so a single comparison
    resolves the rounding whenever C is close enough to B for a wrap to
    be possible at all.

    [1] see the documentation of flint_mpn_mulhigh_n; T. Kime,
        P. Zimmermann, "Exact high part of integer multiplication",
        https://hal.science/hal-04861755 gives sharper variant-specific
        bounds (Corollary 4).
*/
void
_flint_mpn_redc_n(nn_ptr t, nn_srcptr X, nn_srcptr m, nn_srcptr minv,
                  mp_size_t n, nn_ptr scratch)
{
    nn_ptr q = scratch, h = scratch + n;
    ulong C, cy;
    mp_size_t i;

    if (n < 16)
    {
        /* small sizes: exact high product, no guard logic */
        nn_ptr full;
        TMP_INIT;

        for (i = 0; i < n; i++)
            if (X[i] != 0)
                break;
        if (i == n)
        {
            flint_mpn_copyi(t, X + n, n);
            goto reduce;
        }

        flint_mpn_mullow_n(q, X, minv, n);

        TMP_START;
        full = TMP_ALLOC(2 * n * sizeof(ulong));
        flint_mpn_mul_n(full, q, m, n);
        cy = mpn_add_n(t, X + n, full + n, n);
        cy += mpn_add_1(t, t, n, 1);
        TMP_END;

        if (cy)
            mpn_sub_n(t, t, m, n);
        goto reduce;
    }

    for (i = 0; i < n; i++)
        if (X[i] != 0)
            break;
    if (i == n)
    {
        /* q = 0 */
        flint_mpn_copyi(t, X + n, n);
        goto reduce;
    }

    flint_mpn_mullow_n(q, X, minv, n);
    C = flint_mpn_mulhigh_n(h, q, m, n);

    /* a wrap of the true guard past B is possible only for C within
       the error bound of B */
    if (C >= UWORD_MAX - (ulong) (n + 1))
    {
        /* T = top limb of B^n - X_lo, with X_lo != 0 */
        ulong T = ~X[n - 1];
        int lower_zero = 1;

        for (i = 0; i < n - 1; i++)
        {
            if (X[i] != 0)
            {
                lower_zero = 0;
                break;
            }
        }
        T += (ulong) lower_zero;

        if (T < C)
            mpn_add_1(h, h, n, 1);
    }

    cy = mpn_add_n(t, X + n, h, n);
    cy += mpn_add_1(t, t, n, 1);   /* X_lo + lo(q*m) = B^n */

    if (cy)
        mpn_sub_n(t, t, m, n);     /* t < 2m: the borrow matches cy */

reduce:
    if (mpn_cmp(t, m, n) >= 0)
        mpn_sub_n(t, t, m, n);
}
