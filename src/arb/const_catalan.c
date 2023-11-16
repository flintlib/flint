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

void
arb_const_catalan_eval(arb_t s, slong prec)
{
    hypgeom_t series;
    arb_t t;

    arb_init(t);
    hypgeom_init(series);

    fmpz_poly_set_str(series->A, "7  1999553 21620948 94165776 211938912 260619984 166411584 43203456");
    fmpz_poly_set_str(series->B, "1  1");
    fmpz_poly_set_str(series->P, "9  0 0 0 1280 -17536 86400 -195840 207360 -82944");
    fmpz_poly_set_str(series->Q, "9  363825 12034680 150240200 918651040 3101725520 6073920000 6863040000 4147200000 1036800000");

    prec += 4 + FLINT_CLOG2(prec);
    arb_hypgeom_infsum(s, t, series, prec, prec);
    arb_mul_ui(t, t, 2182950, prec);
    arb_div(s, s, t, prec);

    hypgeom_clear(series);
    arb_clear(t);
}

ARB_DEF_CACHED_CONSTANT(arb_const_catalan, arb_const_catalan_eval)

