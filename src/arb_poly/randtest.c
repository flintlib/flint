/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
arb_poly_randtest(arb_poly_t poly, flint_rand_t state, slong len, slong prec, slong mag_bits)
{
    slong i;

    arb_poly_fit_length(poly, len);

    if (n_randint(state, 2))
        for (i = 0; i < len; i++)
            arb_randtest(poly->coeffs + i, state, prec, mag_bits);
    else
        for (i = 0; i < len; i++)
            arb_randtest_precise(poly->coeffs + i, state, prec, mag_bits);

    _arb_poly_set_length(poly, len);
    _arb_poly_normalise(poly);
}

