/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_evaluate2_horner(acb_t y, acb_t z, acb_srcptr poly,
        slong len, const acb_t x, slong prec)
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
    else if (acb_is_zero(x))
    {
        acb_set_round(y, poly + 0, prec);
        acb_set_round(z, poly + 1, prec);
    }
    else if (len == 2)
    {
        acb_mul(y, x, poly + 1, prec);
        acb_add(y, y, poly + 0, prec);
        acb_set_round(z, poly + 1, prec);
    }
    else
    {
        acb_t t, u, v;
        slong i;

        acb_init(t);
        acb_init(u);
        acb_init(v);

        acb_set_round(u, poly + len - 1, prec);
        acb_zero(v);

        for (i = len - 2; i >= 0; i--)
        {
            acb_mul(t, v, x, prec);
            acb_add(v, u, t, prec);
            acb_mul(t, u, x, prec);
            acb_add(u, t, poly + i, prec);
        }

        acb_swap(y, u);
        acb_swap(z, v);

        acb_clear(t);
        acb_clear(u);
        acb_clear(v);
    }
}

void
acb_poly_evaluate2_horner(acb_t r, acb_t s, const acb_poly_t f, const acb_t a, slong prec)
{
    _acb_poly_evaluate2_horner(r, s, f->coeffs, f->length, a, prec);
}

