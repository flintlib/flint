/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_hypgeom.h"

static void
factor(acb_poly_t A, acb_poly_t tmp,
    const acb_poly_struct * a, slong p,
    const acb_poly_t z, slong k, slong len, slong prec)
{
    slong i;

    if (p == 0)
    {
        if (z == NULL)
            acb_poly_one(A);
        else
            acb_poly_set(A, z);
    }
    else
    {
        acb_poly_add_si(A, a, k, prec);

        for (i = 1; i < p; i++)
        {
            acb_poly_add_si(tmp, a + i, k, prec);
            acb_poly_mullow(A, A, tmp, len, prec);
        }

        if (z != NULL)
        {
            acb_poly_mullow(A, A, z, len, prec);
        }
    }
}

static void
bsplit(acb_poly_t A1, acb_poly_t B1, acb_poly_t C1,
        const acb_poly_struct * a, slong p,
        const acb_poly_struct * b, slong q,
        const acb_poly_t z,
        slong aa,
        slong bb,
        slong len, slong prec)
{
    if (bb - aa == 1)
    {
        factor(A1, B1, a, p, z, aa, len, prec);
        factor(C1, B1, b, q, NULL, aa, len, prec);
        /* acb_poly_set(B1, C1);   but we skip this */
    }
    else
    {
        slong m;

        acb_poly_t A2, B2, C2, tmp;

        acb_poly_init(A2);
        acb_poly_init(B2);
        acb_poly_init(C2);
        acb_poly_init(tmp);

        m = aa + (bb - aa) / 2;

        bsplit(A1, B1, C1, a, p, b, q, z, aa, m, len, prec);
        bsplit(A2, B2, C2, a, p, b, q, z, m, bb, len, prec);

        if (bb - m == 1)  /* B2 = C2 */
        {
            if (m - aa == 1)
                acb_poly_add(B2, A1, C1, prec);
            else
                acb_poly_add(B2, A1, B1, prec);

            acb_poly_mullow(B1, B2, C2, len, prec);
        }
        else
        {
            if (m - aa == 1)
            {
                acb_poly_mullow(B1, C1, C2, len, prec);
            }
            else
            {
                acb_poly_mullow(tmp, B1, C2, len, prec);
                acb_poly_swap(B1, tmp);
            }

            acb_poly_mullow(tmp, A1, B2, len, prec);
            acb_poly_add(B1, B1, tmp, prec);
        }

        acb_poly_mullow(tmp, A1, A2, len, prec);
        acb_poly_swap(A1, tmp);
        acb_poly_mullow(tmp, C1, C2, len, prec);
        acb_poly_swap(C1, tmp);

        acb_poly_clear(A2);
        acb_poly_clear(B2);
        acb_poly_clear(C2);
        acb_poly_clear(tmp);
    }
}

void
acb_hypgeom_pfq_series_sum_bs(acb_poly_t s, acb_poly_t t,
    const acb_poly_struct * a, slong p,
    const acb_poly_struct * b, slong q,
    const acb_poly_t z, int regularized,
    slong n, slong len, slong prec)
{
    acb_poly_t u, v, w;
    slong i, start;

    if (n == 0)
    {
        acb_hypgeom_pfq_series_sum_forward(s, t, a, p, b, q, z,
            regularized, n, len, prec);
        return;
    }

    start = 0;

    if (regularized)
    {
        /*
        Use basecase algorithm to get past any poles.

        G(x) = 1/Gamma(x)
        n = 1: s = 0; t = G(b); b = 0 requires skipping
        n = 2: s = G(b); t = G(b) G(b+1); b = 0, -1 require skipping...

        TODO: when pole is near the end, do the *start* (or middle segment)
        using fast algorithm; use basecase algorithm for the end.
        */

        for (i = 0; i < q; i++)
        {
            if (acb_poly_is_zero(b + i))
            {
                start = FLINT_MAX(start, 1);
            }
            else
            {
                /* todo: use a fuzzier test? */
                if (acb_contains_int((b + i)->coeffs) &&
                    !arb_is_positive(acb_realref((b + i)->coeffs)) &&
                    arf_cmpabs_2exp_si(arb_midref(acb_realref((b + i)->coeffs)),
                        FLINT_BITS - 2) < 0)
                {
                    slong c = -arf_get_si(arb_midref(acb_realref((b + i)->coeffs)),
                        ARF_RND_NEAR);

                    /* if c >= n, terminates earlier, so no problem */
                    if (c < n)
                    {
                        start = FLINT_MAX(start, c + 1);
                    }
                }
            }
        }

        /* We should now have start <= n. */
        if (start > n)
            flint_throw(FLINT_ERROR, "(%s)\n", __func__);

        acb_hypgeom_pfq_series_sum_forward(s, t, a, p, b, q, z,
            regularized, start, len, prec);
    }
    else
    {
        acb_poly_zero(s);
        acb_poly_one(t);
    }

    if (start == n)
        return;

    acb_poly_init(u);
    acb_poly_init(v);
    acb_poly_init(w);

    bsplit(u, v, w, a, p, b, q, z, start, n, len, prec);

    if (n - start == 1)
        acb_poly_set(v, w);  /* B1 not set */

    acb_poly_mullow(v, v, t, len, prec);
    acb_poly_div_series(v, v, w, len, prec);
    acb_poly_add(s, s, v, prec);

    acb_poly_mullow(t, t, u, len, prec);
    acb_poly_div_series(t, t, w, len, prec);

    acb_poly_clear(u);
    acb_poly_clear(v);
    acb_poly_clear(w);
}

