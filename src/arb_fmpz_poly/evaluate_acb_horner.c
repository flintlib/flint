/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_fmpz_poly.h"

void
_arb_fmpz_poly_evaluate_acb_horner(acb_t y, const fmpz * f, slong len,
                           const acb_t x, slong prec)
{
    if (len == 0)
    {
        acb_zero(y);
    }
    else if (len == 1 || acb_is_zero(x))
    {
        acb_set_round_fmpz(y, f, prec);
    }
    else if (len == 2)
    {
        acb_mul_fmpz(y, x, f + 1, prec);
        acb_add_fmpz(y, y, f + 0, prec);
    }
    else
    {
        slong i = len - 1;
        acb_t t, u;

        acb_init(t);
        acb_init(u);
        acb_set_fmpz(u, f + i);

        for (i = len - 2; i >= 0; i--)
        {
            acb_mul(t, u, x, prec);
            acb_add_fmpz(u, t, f + i, prec);
        }

        acb_swap(y, u);

        acb_clear(t);
        acb_clear(u);
    }
}

void
arb_fmpz_poly_evaluate_acb_horner(acb_t res, const fmpz_poly_t f, const acb_t a, slong prec)
{
    _arb_fmpz_poly_evaluate_acb_horner(res, f->coeffs, f->length, a, prec);
}

