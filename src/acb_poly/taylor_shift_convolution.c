/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_taylor_shift_convolution(acb_ptr p, const acb_t c, slong len, slong prec)
{
    slong i, n = len - 1;
    acb_t d;
    arb_t f;
    acb_ptr t, u;

    if (acb_is_zero(c) || len <= 1)
        return;

    t = _acb_vec_init(len);
    u = _acb_vec_init(len);

    arb_init(f);
    acb_init(d);

    arb_one(f);
    for (i = 2; i <= n; i++)
    {
        arb_mul_ui(f, f, i, prec);
        acb_mul_arb(p + i, p + i, f, prec);
    }

    _acb_poly_reverse(p, p, len, len);

    acb_one(t + n);
    for (i = n; i > 0; i--)
        acb_mul_ui(t + i - 1, t + i, i, prec);

    if (acb_equal_si(c, -1))
    {
        for (i = 1; i <= n; i += 2)
            acb_neg(t + i, t + i);
    }
    else if (!acb_is_one(c))
    {
        acb_set(d, c);

        for (i = 1; i <= n; i++)
        {
            acb_mul(t + i, t + i, d, prec);
            acb_mul(d, d, c, prec);
        }
    }

    _acb_poly_mullow(u, p, len, t, len, len, prec);

    arb_mul(f, f, f, prec);

    if (arb_bits(f) > 0.25 * prec)
    {
        arb_inv(f, f, prec);
    }
    else
    {
        for (i = 0; i <= n; i++)
            acb_div_arb(u + i, u + i, f, prec);

        arb_one(f);
    }

    for (i = n; i >= 0; i--)
    {
        acb_mul_arb(p + i, u + n - i, f, prec);
        arb_mul_ui(f, f, (i == 0) ? 1 : i, prec);
    }

    _acb_vec_clear(t, len);
    _acb_vec_clear(u, len);

    arb_clear(f);
    acb_clear(d);
}

void
acb_poly_taylor_shift_convolution(acb_poly_t g, const acb_poly_t f,
        const acb_t c, slong prec)
{
    if (f != g)
        acb_poly_set_round(g, f, prec);

    _acb_poly_taylor_shift_convolution(g->coeffs, c, g->length, prec);
}

