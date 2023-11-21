/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_pfq_series_sum_forward(acb_poly_t s, acb_poly_t t,
    const acb_poly_struct * a, slong p,
    const acb_poly_struct * b, slong q,
    const acb_poly_t z, int regularized,
    slong n, slong len, slong prec)
{
    acb_poly_t u, v;
    acb_poly_t tmp;
    slong k, i;

    acb_poly_init(u);
    acb_poly_init(v);
    acb_poly_init(tmp);

    if (!regularized)
    {
        acb_poly_zero(s);
        acb_poly_one(t);

        for (k = 0; k < n && acb_poly_length(t) != 0; k++)
        {
            acb_poly_add(s, s, t, prec);

            if (p > 0)
            {
                acb_poly_add_si(u, a, k, prec);

                for (i = 1; i < p; i++)
                {
                    acb_poly_add_si(v, a + i, k, prec);
                    acb_poly_mullow(u, u, v, len, prec);
                }

                acb_poly_mullow(t, t, u, len, prec);
            }

            if (q > 0)
            {
                acb_poly_add_si(u, b, k, prec);

                for (i = 1; i < q; i++)
                {
                    acb_poly_add_si(v, b + i, k, prec);
                    acb_poly_mullow(u, u, v, len, prec);
                }

                acb_poly_div_series(t, t, u, len, prec);
            }

            acb_poly_mullow(t, t, z, len, prec);
        }
    }
    else
    {
        acb_poly_zero(s);

        if (q == 0)
            acb_poly_one(t);

        for (i = 0; i < q; i++)
        {
            if (i == 0)
            {
                acb_poly_rgamma_series(t, b + i, len, prec);
            }
            else
            {
                acb_poly_rgamma_series(u, b + i, len, prec);
                acb_poly_mullow(tmp, t, u, len, prec);
                acb_poly_swap(tmp, t);
            }
        }

        for (k = 0; k < n; k++)
        {
            acb_poly_add(s, s, t, prec);

            if (p > 0)
            {
                acb_poly_add_si(u, a, k, prec);

                for (i = 1; i < p; i++)
                {
                    acb_poly_add_si(v, a + i, k, prec);
                    acb_poly_mullow(tmp, u, v, len, prec);
                    acb_poly_swap(tmp, u);
                }

                acb_poly_mullow(tmp, t, u, len, prec);
                acb_poly_swap(tmp, t);
            }

            if (q > 0)
            {
                acb_poly_add_si(u, b, k, prec);

                for (i = 1; i < q; i++)
                {
                    acb_poly_add_si(v, b + i, k, prec);
                    acb_poly_mullow(tmp, u, v, len, prec);
                    acb_poly_swap(tmp, u);
                }

                if (acb_poly_length(u) > 0 && !acb_contains_zero(u->coeffs))
                {
                    acb_poly_div_series(tmp, t, u, len, prec);
                    acb_poly_mullow(t, tmp, z, len, prec);
                }
                else
                {
                    /* compute term from scratch */
                    acb_poly_one(t);

                    for (i = 0; i < p; i++)
                    {
                        acb_poly_rising_ui_series(v, a + i, k + 1, len, prec);
                        acb_poly_mullow(t, t, v, len, prec);
                    }

                    for (i = 0; i < q; i++)
                    {
                        acb_poly_add_si(v, b + i, k + 1, prec);
                        acb_poly_rgamma_series(v, v, len, prec);
                        acb_poly_mullow(t, t, v, len, prec);
                    }

                    acb_poly_pow_ui_trunc_binexp(v, z, k + 1, len, prec);
                    acb_poly_mullow(t, t, v, len, prec);
                }
            }
            else
            {
                acb_poly_mullow(tmp, t, z, len, prec);
                acb_poly_swap(tmp, t);
            }
        }
    }

    acb_poly_clear(u);
    acb_poly_clear(v);
    acb_poly_clear(tmp);
}

