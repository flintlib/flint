/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "arb.h"
#include "arb_fmpz_poly.h"

void
_arb_fmpz_poly_evaluate_arb_rectangular(arb_t y, const fmpz * poly,
    slong len, const arb_t x, slong prec)
{
    slong i, m, r;
    arb_ptr xs;
    arb_t s, t, c;

    if (len < 3)
    {
        _arb_fmpz_poly_evaluate_arb_horner(y, poly, len, x, prec);
        return;
    }

    m = n_sqrt(len) + 1;
    r = (len + m - 1) / m;

    xs = _arb_vec_init(m + 1);
    arb_init(s);
    arb_init(t);
    arb_init(c);

    _arb_vec_set_powers(xs, x, m + 1, prec);

    arb_set_fmpz(y, poly + (r - 1) * m);
    arb_dot_fmpz(y, y, 0, xs + 1, 1,
        poly + (r - 1) * m + 1, 1, len - (r - 1) * m - 1, prec);

    for (i = r - 2; i >= 0; i--)
    {
        arb_set_fmpz(s, poly + i * m);
        arb_dot_fmpz(s, s, 0, xs + 1, 1,
            poly + i * m + 1, 1, m - 1, prec);
        arb_mul(y, y, xs + m, prec);
        arb_add(y, y, s, prec);
    }

    _arb_vec_clear(xs, m + 1);
    arb_clear(s);
    arb_clear(t);
    arb_clear(c);
}

void
arb_fmpz_poly_evaluate_arb_rectangular(arb_t res, const fmpz_poly_t f, const arb_t a, slong prec)
{
    _arb_fmpz_poly_evaluate_arb_rectangular(res, f->coeffs, f->length, a, prec);
}

