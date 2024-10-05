/*
    Copyright (C) 2023 Joel Dahne

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr_poly.h"

int
_gr_poly_nth_derivative(gr_ptr res, gr_srcptr poly, ulong n, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op_fmpz mul_fmpz = GR_BINARY_OP_FMPZ(ctx, MUL_FMPZ);
    fmpz_t c;
    slong i;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    /* todo: optimize for small n so that c fits in ulong */
    /* todo: optimize for large n and finite characteristic or precision */

    fmpz_init(c);

    fmpz_fac_ui(c, n);
    status |= mul_fmpz(GR_ENTRY(res, 0, sz), GR_ENTRY(poly, n, sz), c, ctx);
    for (i = n + 1; i < len; i ++)
    {
        fmpz_divexact_ui(c, c, i - n);
        fmpz_mul_ui(c, c, i);
        status |= mul_fmpz(GR_ENTRY(res, i - n, sz), GR_ENTRY(poly, i, sz), c, ctx);
    }

    fmpz_clear(c);

    return status;
}

int
gr_poly_nth_derivative(gr_poly_t res, const gr_poly_t poly, ulong n, gr_ctx_t ctx)
{
    int status;
    const slong len = poly->length;

    if ((ulong) len <= n)
    {
        status = gr_poly_zero(res, ctx);
    }
    else if (n == 0)
    {
        status = gr_poly_set(res, poly, ctx);
    }
    else if (n == 1)
    {
        gr_poly_fit_length(res, len - 1, ctx);
        status = _gr_poly_derivative(res->coeffs, poly->coeffs, len, ctx);
        _gr_poly_set_length(res, len - 1, ctx);
        /* todo: only call in nonzero characteristic */
        _gr_poly_normalise(res, ctx);
    }
    else
    {
        gr_poly_fit_length(res, len - n, ctx);
        status = _gr_poly_nth_derivative(res->coeffs, poly->coeffs, n, len, ctx);
        _gr_poly_set_length(res, len - n, ctx);
        /* todo: only call in nonzero characteristic */
        _gr_poly_normalise(res, ctx);
    }

    return status;
}
