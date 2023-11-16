/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void _acb_poly_acb_invpow_cpx(acb_ptr res, const acb_t N, const acb_t c, slong trunc, slong prec)
{
    slong i;
    acb_t logN;

    acb_init(logN);
    acb_log(logN, N, prec);
    acb_mul(res + 0, logN, c, prec);
    acb_neg(res + 0, res + 0);
    acb_exp(res + 0, res + 0, prec);

    for (i = 1; i < trunc; i++)
    {
        acb_mul(res + i, res + i - 1, logN, prec);
        acb_div_si(res + i, res + i, -i, prec);
    }

    acb_clear(logN);
}

void
_acb_poly_powsum_series_naive(acb_ptr z,
    const acb_t s, const acb_t a, const acb_t q, slong n, slong len, slong prec)
{
    slong k, i;
    int q_one, s_int;
    acb_t ak, logak, t, qpow, negs;

    acb_init(ak);
    acb_init(logak);
    acb_init(t);
    acb_init(qpow);
    acb_init(negs);

    _acb_vec_zero(z, len);
    acb_one(qpow);
    acb_neg(negs, s);

    q_one = acb_is_one(q);
    s_int = arb_is_int(acb_realref(s)) && arb_is_zero(acb_imagref(s));

    for (k = 0; k < n; k++)
    {
        acb_add_ui(ak, a, k, prec);

        if (len == 1)
        {
            acb_pow(t, ak, negs, prec);
        }
        else
        {
            acb_log(logak, ak, prec);

            if (s_int)
            {
                acb_pow(t, ak, negs, prec);
            }
            else
            {
                acb_mul(t, logak, negs, prec);
                acb_exp(t, t, prec);
            }
        }

        if (!q_one)
        {
            acb_mul(t, t, qpow, prec);
            if (k < n - 1)
                acb_mul(qpow, qpow, q, prec);
        }

        acb_add(z, z, t, prec);

        for (i = 1; i < len; i++)
        {
            acb_mul(t, t, logak, prec);
            acb_div_si(t, t, -i, prec);
            acb_add(z + i, z + i, t, prec);
        }
    }

    acb_clear(ak);
    acb_clear(logak);
    acb_clear(t);
    acb_clear(qpow);
    acb_clear(negs);
}

