/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_evaluate2_horner(arb_t y, arb_t z, arb_srcptr poly,
        slong len, const arb_t x, slong prec)
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
    else if (arb_is_zero(x))
    {
        arb_set_round(y, poly + 0, prec);
        arb_set_round(z, poly + 1, prec);
    }
    else if (len == 2)
    {
        arb_mul(y, x, poly + 1, prec);
        arb_add(y, y, poly + 0, prec);
        arb_set_round(z, poly + 1, prec);
    }
    else
    {
        arb_t t, u, v;
        slong i;

        arb_init(t);
        arb_init(u);
        arb_init(v);

        arb_set_round(u, poly + len - 1, prec);
        arb_zero(v);

        for (i = len - 2; i >= 0; i--)
        {
            arb_mul(t, v, x, prec);
            arb_add(v, u, t, prec);
            arb_mul(t, u, x, prec);
            arb_add(u, t, poly + i, prec);
        }

        arb_swap(y, u);
        arb_swap(z, v);

        arb_clear(t);
        arb_clear(u);
        arb_clear(v);
    }
}

void
arb_poly_evaluate2_horner(arb_t r, arb_t s, const arb_poly_t f, const arb_t a, slong prec)
{
    _arb_poly_evaluate2_horner(r, s, f->coeffs, f->length, a, prec);
}

