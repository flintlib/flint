/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "gr_poly.h"

#ifdef __GNUC__
# define sqrt __builtin_sqrt
#else
# include <math.h>
#endif

int
_gr_arb_poly_taylor_shift(arb_ptr res, arb_srcptr poly, slong n, const arb_t c, gr_ctx_t ctx)
{
    slong prec;

    if (n <= 30)
        return _gr_poly_taylor_shift_horner(res, poly, n, c, ctx);

    prec = _gr_ctx_get_real_prec(ctx);

    if (n <= 30 || (n <= 500 && arb_bits(c) == 1 && n < 30 + 3 * sqrt(prec))
                || (n <= 100 && arb_bits(c) < 0.01 * prec))
    {
        return _gr_poly_taylor_shift_horner(res, poly, n, c, ctx);
    }
    else if (prec > 2 * n)
    {
        return _gr_poly_taylor_shift_convolution(res, poly, n, c, ctx);
    }
    else
    {
        return _gr_poly_taylor_shift_divconquer(res, poly, n, c, ctx);
    }
}

void
_arb_poly_taylor_shift(arb_ptr poly, const arb_t c, slong n, slong prec)
{
    gr_ctx_t ctx;
    gr_ctx_init_real_arb(ctx, prec);
    if (_gr_poly_taylor_shift(poly, poly, n, c, ctx) != GR_SUCCESS)
        _arb_vec_indeterminate(poly, n);
}

void
arb_poly_taylor_shift(arb_poly_t g, const arb_poly_t f,
    const arb_t c, slong prec)
{
    if (f != g)
        arb_poly_set_round(g, f, prec);

    _arb_poly_taylor_shift(g->coeffs, c, g->length, prec);
}

