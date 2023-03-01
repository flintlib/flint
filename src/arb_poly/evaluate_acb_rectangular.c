/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_evaluate_acb_rectangular(acb_t y, arb_srcptr poly,
    slong len, const acb_t x, slong prec)
{
    slong i, j, m, r;
    acb_ptr xs;
    acb_t s, t, c;

    if (len < 3)
    {
        if (len == 0)
        {
            acb_zero(y);
        }
        else if (len == 1)
        {
            acb_set_round_arb(y, poly + 0, prec);
        }
        else if (len == 2)
        {
            acb_mul_arb(y, x, poly + 1, prec);
            acb_add_arb(y, y, poly + 0, prec);
        }
        return;
    }

    m = n_sqrt(len) + 1;
    r = (len + m - 1) / m;

    xs = _acb_vec_init(m + 1);
    acb_init(s);
    acb_init(t);
    acb_init(c);

    _acb_vec_set_powers(xs, x, m + 1, prec);

    acb_set_arb(y, poly + (r - 1) * m);
    for (j = 1; (r - 1) * m + j < len; j++)
        acb_addmul_arb(y, xs + j, poly + (r - 1) * m + j, prec);

    for (i = r - 2; i >= 0; i--)
    {
        acb_set_arb(s, poly + i * m);
        for (j = 1; j < m; j++)
            acb_addmul_arb(s, xs + j, poly + i * m + j, prec);

        acb_mul(y, y, xs + m, prec);
        acb_add(y, y, s, prec);
    }

    _acb_vec_clear(xs, m + 1);
    acb_clear(s);
    acb_clear(t);
    acb_clear(c);
}

void
arb_poly_evaluate_acb_rectangular(acb_t res, const arb_poly_t f, const acb_t a, slong prec)
{
    _arb_poly_evaluate_acb_rectangular(res, f->coeffs, f->length, a, prec);
}

