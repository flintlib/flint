/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_evaluate_horner(arb_t y, arb_srcptr f, slong len,
                           const arb_t x, slong prec)
{
    if (len == 0)
    {
        arb_zero(y);
    }
    else if (len == 1 || arb_is_zero(x))
    {
        arb_set_round(y, f, prec);
    }
    else if (len == 2)
    {
        arb_mul(y, x, f + 1, prec);
        arb_add(y, y, f + 0, prec);
    }
    else
    {
        slong i = len - 1;
        arb_t t, u;

        arb_init(t);
        arb_init(u);
        arb_set(u, f + i);

        for (i = len - 2; i >= 0; i--)
        {
            arb_mul(t, u, x, prec);
            arb_add(u, f + i, t, prec);
        }

        arb_swap(y, u);

        arb_clear(t);
        arb_clear(u);
    }
}

void
arb_poly_evaluate_horner(arb_t res, const arb_poly_t f, const arb_t a, slong prec)
{
    _arb_poly_evaluate_horner(res, f->coeffs, f->length, a, prec);
}

