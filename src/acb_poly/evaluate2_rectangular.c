/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "acb_poly.h"

void
_acb_poly_evaluate2_rectangular(acb_t y, acb_t z, acb_srcptr poly,
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
            acb_zero(z);
        }
        else if (len == 1)
        {
            acb_set_round(y, poly + 0, prec);
            acb_zero(z);
        }
        else if (len == 2)
        {
            acb_mul(y, x, poly + 1, prec);
            acb_add(y, y, poly + 0, prec);
            acb_set_round(z, poly + 1, prec);
        }
        return;
    }

    m = n_sqrt(len) + 1;
    m *= 1;

    r = (len + m - 1) / m;

    xs = _acb_vec_init(m + 1);
    acb_init(s);
    acb_init(t);
    acb_init(c);

    _acb_vec_set_powers(xs, x, m + 1, prec);
    acb_dot(y, poly + (r - 1) * m, 0, xs + 1, 1,
        poly + (r - 1) * m + 1, 1, len - (r - 1) * m - 1, prec);

    for (i = r - 2; i >= 0; i--)
    {
        acb_dot(s, poly + i * m, 0, xs + 1, 1,
            poly + i * m + 1, 1, m - 1, prec);
        acb_mul(y, y, xs + m, prec);
        acb_add(y, y, s, prec);
    }

    /* todo: rewrite using acb_dot */
    len -= 1;
    r = (len + m - 1) / m;
    acb_mul_ui(z, poly + (r - 1) * m + 1, (r - 1) * m + 1, ARF_PREC_EXACT);
    for (j = 1; (r - 1) * m + j < len; j++)
    {
        acb_mul_ui(c, poly + (r - 1) * m + j + 1, (r - 1) * m + j + 1, ARF_PREC_EXACT);
        acb_addmul(z, xs + j, c, prec);
    }

    for (i = r - 2; i >= 0; i--)
    {
        acb_mul_ui(s, poly + i * m + 1, i * m + 1, ARF_PREC_EXACT);

        for (j = 1; j < m; j++)
        {
            acb_mul_ui(c, poly + i * m + j + 1, i * m + j + 1, ARF_PREC_EXACT);
            acb_addmul(s, xs + j, c, prec);
        }

        acb_mul(z, z, xs + m, prec);
        acb_add(z, z, s, prec);
    }

    _acb_vec_clear(xs, m + 1);
    acb_clear(s);
    acb_clear(t);
    acb_clear(c);
}

void
acb_poly_evaluate2_rectangular(acb_t r, acb_t s, const acb_poly_t f, const acb_t a, slong prec)
{
    _acb_poly_evaluate2_rectangular(r, s, f->coeffs, f->length, a, prec);
}

