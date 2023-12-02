/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "acb_poly.h"
#include "acb_hypgeom.h"

#define FAST_T 0

void acb_poly_reciprocal_majorant(arb_poly_t res, const acb_poly_t poly, slong prec);

static void
rsplit(acb_poly_t res, acb_poly_t term,
    const acb_poly_struct * a, slong p,
    const acb_poly_struct * b, slong q,
    const acb_poly_t z, slong offset, slong n, slong len, slong prec)
{
    acb_poly_struct * zpow;
    acb_poly_t s, t, u;
    slong i, j, k, m, tprec;
    int is_real;

#if FAST_T
    arb_poly_t B, C, D;
#else
    acb_poly_t B, C, D;
#endif

    if (n == 0)
    {
        acb_poly_zero(res);
        acb_poly_one(term);
        return;
    }

    if (n < 0)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    m = n_sqrt(n);
    m = FLINT_MIN(m, 150);

    acb_poly_init(s);
    acb_poly_init(t);
    acb_poly_init(u);

#if FAST_T
    tprec = MAG_BITS;

    arb_poly_init(B);
    arb_poly_init(C);
    arb_poly_init(D);

    arb_poly_one(B);
    arb_poly_one(D);

    is_real = 1;

    for (i = 0; i < p; i++)
        for (j = 0; j < FLINT_MIN(acb_poly_length(a + i), len); j++)
            is_real = is_real && acb_is_real((a + i)->coeffs + j);

    for (i = 0; i < q; i++)
        for (j = 0; j < FLINT_MIN(acb_poly_length(b + i), len); j++)
            is_real = is_real && acb_is_real((b + i)->coeffs + j);

    for (j = 0; j < FLINT_MIN(acb_poly_length(z), len); j++)
        is_real = is_real && acb_is_real(z->coeffs + j);
#else
    tprec = 2 * MAG_BITS;
    is_real = 0; (void) is_real;

    acb_poly_init(B);
    acb_poly_init(C);
    acb_poly_init(D);

    acb_poly_one(B);
    acb_poly_one(D);
#endif

    zpow = flint_malloc(sizeof(acb_poly_struct) * (m + 1));

    for (i = 0; i <= m; i++)
        acb_poly_init(zpow + i);

    for (i = 0; i <= m; i++)
    {
        if (i == 0)
            acb_poly_one(zpow + i);
        else if (i == 1)
            acb_poly_set_round(zpow + i, z, prec);
        else if (i % 2 == 0)
            acb_poly_mullow(zpow + i, zpow + i / 2, zpow + i / 2, len, prec);
        else
            acb_poly_mullow(zpow + i, zpow + i - 1, zpow + 1, len, prec);
    }

    for (k = n; k >= 0; k--)
    {
        j = k % m;

        if (k < n)
            acb_poly_add(s, s, zpow + j, prec);

        if (k > 0)
        {
            if (p > 0)
            {
                acb_poly_add_si(u, a, offset + k - 1, prec);

                for (i = 1; i < p; i++)
                {
                    acb_poly_add_si(t, a + i, offset + k - 1, prec);
                    acb_poly_mullow(u, u, t, len, prec);
                }

                if (k < n)
                    acb_poly_mullow(s, s, u, len, prec);

#if FAST_T
                acb_poly_majorant(C, u, tprec);
                arb_poly_mullow(B, B, C, len, tprec);
#else
                acb_poly_set_round(C, u, tprec);
                acb_poly_mullow(B, B, C, len, tprec);

#endif
            }

            if (q > 0)
            {
                acb_poly_add_si(u, b, offset + k - 1, prec);

                for (i = 1; i < q; i++)
                {
                    acb_poly_add_si(t, b + i, offset + k - 1, prec);
                    acb_poly_mullow(u, u, t, len, prec);
                }

                if (k < n)
                    acb_poly_div_series(s, s, u, len, prec);

#if FAST_T
                acb_poly_reciprocal_majorant(C, u, tprec);
                arb_poly_mullow(D, D, C, len, tprec);
#else
                acb_poly_set_round(C, u, tprec);
                acb_poly_mullow(D, D, C, len, tprec);
#endif
            }

            if (j == 0 && k < n)
            {
                acb_poly_mullow(s, s, zpow + m, len, prec);
            }
        }
    }

#if FAST_T
    arb_poly_div_series(B, B, D, len, tprec);
    acb_poly_majorant(C, z, tprec);
    arb_poly_pow_ui_trunc_binexp(C, C, n, len, tprec);
    arb_poly_mullow(C, B, C, len, tprec);

    acb_poly_fit_length(term, arb_poly_length(C));

    for (i = 0; i < arb_poly_length(C); i++)
    {
        arb_get_mag(arb_radref(acb_realref(term->coeffs + i)), C->coeffs + i);
        arf_zero(arb_midref(acb_realref(term->coeffs + i)));

        if (is_real)
            arb_zero(acb_imagref(term->coeffs + i));
        else
            arb_set(acb_imagref(term->coeffs + i), acb_realref(term->coeffs + i));
    }

    _acb_poly_set_length(term, arb_poly_length(C));
    _acb_poly_normalise(term);

#else
    acb_poly_div_series(B, B, D, len, tprec);
    acb_poly_set_round(C, z, tprec);
    acb_poly_pow_ui_trunc_binexp(C, C, n, len, tprec);
    acb_poly_mullow(term, B, C, len, tprec);
#endif

    acb_poly_set(res, s);

#if FAST_T
    arb_poly_clear(B);
    arb_poly_clear(C);
    arb_poly_clear(D);
#else
    acb_poly_clear(B);
    acb_poly_clear(C);
    acb_poly_clear(D);
#endif

    acb_poly_clear(s);
    acb_poly_clear(t);
    acb_poly_clear(u);

    for (i = 0; i <= m; i++)
        acb_poly_clear(zpow + i);
    flint_free(zpow);
}

void
acb_hypgeom_pfq_series_sum_rs(acb_poly_t s, acb_poly_t t,
    const acb_poly_struct * a, slong p,
    const acb_poly_struct * b, slong q,
    const acb_poly_t z, int regularized,
    slong n, slong len, slong prec)
{
    acb_poly_t u, v;
    slong   i, start;

    if (n == 0)
    {
        acb_hypgeom_pfq_series_sum_forward(s, t, a, p, b, q, z,
            regularized, n, len, prec);
        return;
    }

    start = 0;

    if (regularized)
    {
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

    rsplit(u, v, a, p, b, q, z, start, n - start, len, prec);

    acb_poly_mullow(u, u, t, len, prec);
    acb_poly_add(s, s, u, prec);

    acb_poly_mullow(t, t, v, len, prec);

    acb_poly_clear(u);
    acb_poly_clear(v);
}

