/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "arb.h"
#include "hypgeom.h"

static void
atanh_bsplit(arb_t s, ulong c, slong a, slong prec)
{
    arb_t t;
    hypgeom_t series;
    hypgeom_init(series);
    arb_init(t);

    fmpz_poly_set_ui(series->A, 1);
    fmpz_poly_set_coeff_ui(series->B, 0, 1);
    fmpz_poly_set_coeff_ui(series->B, 1, 2);
    fmpz_poly_set_ui(series->P, 1);
    fmpz_poly_set_ui(series->Q, c * c);

    arb_hypgeom_infsum(s, t, series, prec, prec);
    arb_mul_si(s, s, a, prec);
    arb_mul_ui(t, t, c, prec);
    arb_div(s, s, t, prec);

    arb_clear(t);
    hypgeom_clear(series);
}

void
arb_const_log10_eval(arb_t s, slong prec)
{
    arb_t t;
    arb_init(t);

    prec += FLINT_CLOG2(prec);

    atanh_bsplit(s, 31, 46, prec);
    atanh_bsplit(t, 49, 34, prec);
    arb_add(s, s, t, prec);
    atanh_bsplit(t, 161, 20, prec);
    arb_add(s, s, t, prec);

    arb_clear(t);
}

ARB_DEF_CACHED_CONSTANT(arb_const_log10, arb_const_log10_eval)

