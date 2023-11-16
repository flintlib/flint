/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_pow_ui_binexp(gr_ptr res,
    gr_srcptr f, slong flen, ulong exp, gr_ctx_t ctx)
{
    return _gr_poly_pow_series_ui_binexp(res, f, flen, exp, exp * (flen - 1) + 1, ctx);
}

int
gr_poly_pow_ui_binexp(gr_poly_t res,
    const gr_poly_t poly, ulong exp, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong flen, rlen;

    flen = poly->length;

    if (exp == 0)
    {
        return gr_poly_one(res, ctx);
    }
    else if (flen == 0)
    {
        return gr_poly_zero(res, ctx);
    }
    else
    {
        ulong hi, lo;

        umul_ppmm(hi, lo, exp, flen - 1);

        if (hi != 0 || lo >= WORD_MAX)
            return GR_UNABLE;

        rlen = lo + 1;

        if (res != poly)
        {
            gr_poly_fit_length(res, rlen, ctx);
            status |= _gr_poly_pow_ui_binexp(res->coeffs, poly->coeffs, flen, exp, ctx);
            _gr_poly_set_length(res, rlen, ctx);
            _gr_poly_normalise(res, ctx);
        }
        else
        {
            gr_poly_t t;
            gr_poly_init2(t, rlen, ctx);
            status |= _gr_poly_pow_ui_binexp(t->coeffs, poly->coeffs, flen, exp, ctx);
            _gr_poly_set_length(t, rlen, ctx);
            _gr_poly_normalise(t, ctx);
            gr_poly_swap(res, t, ctx);
            gr_poly_clear(t, ctx);
        }

        return status;
    }
}
