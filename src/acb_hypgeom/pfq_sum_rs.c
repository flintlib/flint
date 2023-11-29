/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "acb.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_pfq_sum_rs(acb_t res, acb_t term, acb_srcptr a, slong p,
                                              acb_srcptr b, slong q, const acb_t z, slong n, slong prec)
{
    acb_ptr zpow;
    acb_t s, t, u;
    slong i, j, k, m;
    mag_t B, C;

    if (n == 0)
    {
        acb_zero(res);
        acb_one(term);
        return;
    }

    if (n < 0)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    m = n_sqrt(n);
    m = FLINT_MIN(m, 150);

    mag_init(B);
    mag_init(C);
    acb_init(s);
    acb_init(t);
    acb_init(u);
    zpow = _acb_vec_init(m + 1);

    _acb_vec_set_powers(zpow, z, m + 1, prec);

    mag_one(B);

    for (k = n; k >= 0; k--)
    {
        j = k % m;

        if (k < n)
            acb_add(s, s, zpow + j, prec);

        if (k > 0)
        {
            if (p > 0)
            {
                acb_add_ui(u, a, k - 1, prec);

                for (i = 1; i < p; i++)
                {
                    acb_add_ui(t, a + i, k - 1, prec);
                    acb_mul(u, u, t, prec);
                }

                if (k < n)
                    acb_mul(s, s, u, prec);

                acb_get_mag(C, u);
                mag_mul(B, B, C);
            }

            if (q > 0)
            {
                acb_add_ui(u, b, k - 1, prec);

                for (i = 1; i < q; i++)
                {
                    acb_add_ui(t, b + i, k - 1, prec);
                    acb_mul(u, u, t, prec);
                }

                if (k < n)
                    acb_div(s, s, u, prec);

                acb_get_mag_lower(C, u);
                mag_div(B, B, C);
            }

            if (j == 0 && k < n)
            {
                acb_mul(s, s, zpow + m, prec);
            }
        }
    }

    acb_get_mag(C, z);
    mag_pow_ui(C, C, n);
    mag_mul(B, B, C);

    acb_zero(term);
    if (_acb_vec_is_real(a, p) && _acb_vec_is_real(b, q) && acb_is_real(z))
        arb_add_error_mag(acb_realref(term), B);
    else
        acb_add_error_mag(term, B);

    acb_set(res, s);

    mag_clear(B);
    mag_clear(C);
    acb_clear(s);
    acb_clear(t);
    acb_clear(u);
    _acb_vec_clear(zpow, m + 1);
}

