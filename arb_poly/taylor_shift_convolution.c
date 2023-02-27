/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_taylor_shift_convolution(arb_ptr p, const arb_t c, slong len, slong prec)
{
    slong i, n = len - 1;
    arb_t f, d;
    arb_ptr t, u;

    if (arb_is_zero(c) || len <= 1)
        return;

    t = _arb_vec_init(len);
    u = _arb_vec_init(len);

    arb_init(f);
    arb_init(d);

    arb_one(f);
    for (i = 2; i <= n; i++)
    {
        arb_mul_ui(f, f, i, prec);
        arb_mul(p + i, p + i, f, prec);
    }

    _arb_poly_reverse(p, p, len, len);

    arb_one(t + n);
    for (i = n; i > 0; i--)
        arb_mul_ui(t + i - 1, t + i, i, prec);

    if (arb_equal_si(c, -1))
    {
        for (i = 1; i <= n; i += 2)
            arb_neg(t + i, t + i);
    }
    else if (!arb_is_one(c))
    {
        arb_set(d, c);

        for (i = 1; i <= n; i++)
        {
            arb_mul(t + i, t + i, d, prec);
            arb_mul(d, d, c, prec);
        }
    }

    _arb_poly_mullow(u, p, len, t, len, len, prec);

    arb_mul(f, f, f, prec);

    if (arb_bits(f) > 0.25 * prec)
    {
        arb_inv(f, f, prec);
    }
    else
    {
        for (i = 0; i <= n; i++)
            arb_div(u + i, u + i, f, prec);

        arb_one(f);
    }

    for (i = n; i >= 0; i--)
    {
        arb_mul(p + i, u + n - i, f, prec);
        arb_mul_ui(f, f, (i == 0) ? 1 : i, prec);
    }

    _arb_vec_clear(t, len);
    _arb_vec_clear(u, len);

    arb_clear(f);
    arb_clear(d);
}

void
arb_poly_taylor_shift_convolution(arb_poly_t g, const arb_poly_t f,
        const arb_t c, slong prec)
{
    if (f != g)
        arb_poly_set_round(g, f, prec);

    _arb_poly_taylor_shift_convolution(g->coeffs, c, g->length, prec);
}

