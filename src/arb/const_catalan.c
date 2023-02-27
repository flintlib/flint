/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "hypgeom.h"

void
arb_const_catalan_eval(arb_t s, slong prec)
{
    hypgeom_t series;
    arb_t t;

    arb_init(t);
    hypgeom_init(series);

    fmpz_poly_set_str(series->P, "5  32 160 288 224 64");
    fmpz_poly_set_str(series->Q, "5  225 3240 14904 23328 11664");
    fmpz_poly_set_str(series->A, "3  411 976 580");
    fmpz_poly_set_str(series->B, "5  1 5 9 7 2");

    prec += FLINT_CLOG2(prec);
    arb_hypgeom_infsum(s, t, series, prec, prec);
    arb_mul_ui(t, t, 450, prec);
    arb_div(s, s, t, prec);

    hypgeom_clear(series);
    arb_clear(t);
}

ARB_DEF_CACHED_CONSTANT(arb_const_catalan, arb_const_catalan_eval)

