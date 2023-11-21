/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "arb_poly.h"

void
_arb_poly_evaluate2_rectangular(arb_t y, arb_t z, arb_srcptr poly,
    slong len, const arb_t x, slong prec)
{
    slong i, j, m, r;
    arb_ptr xs;
    arb_t s, t, c;

    if (len < 3)
    {
        if (len == 0)
        {
            arb_zero(y);
            arb_zero(z);
        }
        else if (len == 1)
        {
            arb_set_round(y, poly + 0, prec);
            arb_zero(z);
        }
        else if (len == 2)
        {
            arb_mul(y, x, poly + 1, prec);
            arb_add(y, y, poly + 0, prec);
            arb_set_round(z, poly + 1, prec);
        }
        return;
    }

    m = n_sqrt(len) + 1;
    m *= 1;

    r = (len + m - 1) / m;

    xs = _arb_vec_init(m + 1);
    arb_init(s);
    arb_init(t);
    arb_init(c);

    _arb_vec_set_powers(xs, x, m + 1, prec);

    arb_dot(y, poly + (r - 1) * m, 0, xs + 1, 1,
        poly + (r - 1) * m + 1, 1, len - (r - 1) * m - 1, prec);

    for (i = r - 2; i >= 0; i--)
    {
        arb_dot(s, poly + i * m, 0, xs + 1, 1,
            poly + i * m + 1, 1, m - 1, prec);
        arb_mul(y, y, xs + m, prec);
        arb_add(y, y, s, prec);
    }

    /* todo: rewrite using arb_dot */
    len -= 1;
    r = (len + m - 1) / m;
    arb_mul_ui(z, poly + (r - 1) * m + 1, (r - 1) * m + 1, ARF_PREC_EXACT);
    for (j = 1; (r - 1) * m + j < len; j++)
    {
        arb_mul_ui(c, poly + (r - 1) * m + j + 1, (r - 1) * m + j + 1, ARF_PREC_EXACT);
        arb_addmul(z, xs + j, c, prec);
    }

    for (i = r - 2; i >= 0; i--)
    {
        arb_mul_ui(s, poly + i * m + 1, i * m + 1, ARF_PREC_EXACT);

        for (j = 1; j < m; j++)
        {
            arb_mul_ui(c, poly + i * m + j + 1, i * m + j + 1, ARF_PREC_EXACT);
            arb_addmul(s, xs + j, c, prec);
        }

        arb_mul(z, z, xs + m, prec);
        arb_add(z, z, s, prec);
    }

    _arb_vec_clear(xs, m + 1);
    arb_clear(s);
    arb_clear(t);
    arb_clear(c);
}

void
arb_poly_evaluate2_rectangular(arb_t r, arb_t s, const arb_poly_t f, const arb_t a, slong prec)
{
    _arb_poly_evaluate2_rectangular(r, s, f->coeffs, f->length, a, prec);
}

