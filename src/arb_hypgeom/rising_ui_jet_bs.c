/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "arb_hypgeom.h"

static void
bsplit(arb_ptr res, const arb_t x, ulong a, ulong b, slong trunc, slong prec)
{
    trunc = FLINT_MIN(trunc, b - a + 1);

    if (b - a <= 12)
    {
        if (a == 0)
        {
            arb_hypgeom_rising_ui_jet_powsum(res, x, b - a, FLINT_MIN(trunc, b - a + 1), prec);
        }
        else
        {
            arb_t t;
            arb_init(t);
            arb_add_ui(t, x, a, prec);
            arb_hypgeom_rising_ui_jet_powsum(res, t, b - a, FLINT_MIN(trunc, b - a + 1), prec);
            arb_clear(t);
        }
    }
    else
    {
        arb_ptr L, R;
        slong len1, len2;

        slong m = a + (b - a) / 2;

        len1 = poly_pow_length(2, m - a, trunc);
        len2 = poly_pow_length(2, b - m, trunc);

        L = _arb_vec_init(len1 + len2);
        R = L + len1;

        bsplit(L, x, a, m, trunc, prec);
        bsplit(R, x, m, b, trunc, prec);

        _arb_poly_mullow(res, L, len1, R, len2,
            FLINT_MIN(trunc, len1 + len2 - 1), prec);

        _arb_vec_clear(L, len1 + len2);
    }
}

void
arb_hypgeom_rising_ui_jet_bs(arb_ptr res, const arb_t x, ulong n, slong len, slong prec)
{
    if (len == 0)
        return;

    if (len > n + 1)
    {
        _arb_vec_zero(res + n + 1, len - n - 1);
        len = n + 1;
    }

    if (len == n + 1)
    {
        arb_one(res + n);
        len = n;
    }

    if (n <= 1)
    {
        if (n == 1)
            arb_set_round(res, x, prec);
        return;
    }

    bsplit(res, x, 0, n, len, prec);
}

