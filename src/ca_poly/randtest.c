/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

void
ca_poly_randtest(ca_poly_t poly, flint_rand_t state, slong len, slong depth, slong bits, ca_ctx_t ctx)
{
    slong i;

    ca_poly_fit_length(poly, len, ctx);

    for (i = 0; i < len; i++)
        ca_randtest(poly->coeffs + i, state, depth, bits, ctx);

    _ca_poly_set_length(poly, len, ctx);
    _ca_poly_normalise(poly, ctx);
}

void
ca_poly_randtest_rational(ca_poly_t poly, flint_rand_t state, slong len, slong bits, ca_ctx_t ctx)
{
    slong i;

    ca_poly_fit_length(poly, len, ctx);

    for (i = 0; i < len; i++)
        ca_randtest_rational(poly->coeffs + i, state, bits, ctx);

    _ca_poly_set_length(poly, len, ctx);
    _ca_poly_normalise(poly, ctx);
}
