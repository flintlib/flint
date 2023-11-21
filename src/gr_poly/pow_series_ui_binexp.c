/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "gr_vec.h"
#include "gr_poly.h"


#define MUL(z, zlen, x, xlen, y, ylen, trunc, ctx) \
    do { \
        slong slen = FLINT_MIN(xlen + ylen - 1, trunc); \
        status |= _gr_poly_mullow(z, x, xlen, y, ylen, slen, ctx); \
        zlen = slen; \
    } while (0)

int
_gr_poly_pow_series_ui_binexp(gr_ptr res,
    gr_srcptr f, slong flen, ulong exp, slong len, gr_ctx_t ctx)
{
    gr_ptr v, R, S, T;
    slong rlen;
    ulong bit;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (exp <= 1)
    {
        if (exp == 0)
            return gr_one(res, ctx);
        else
            return _gr_vec_set(res, f, len, ctx);
    }

    /* (f * x^r)^m = x^(rm) * f^m */
    while (flen > 1 && gr_is_zero(f, ctx) == T_TRUE)
    {
        if (((ulong) len) > exp)
        {
            status |= _gr_vec_zero(res, exp, ctx);
            len -= exp;
            res = GR_ENTRY(res, exp, sz);
        }
        else
        {
            status |= _gr_vec_zero(res, len, ctx);
            return status;
        }

        f = GR_ENTRY(f, 1, sz);
        flen--;
    }

    if (exp == 2)
    {
        status |= _gr_poly_mullow(res, f, flen, f, flen, len, ctx);
        return status;
    }

    if (flen == 1)
    {
        status |= gr_pow_ui(res, f, exp, ctx);
        return status;
    }

    GR_TMP_INIT_VEC(v, len, ctx);

    bit = UWORD(1) << (FLINT_BIT_COUNT(exp) - 2);

    if (n_zerobits(exp) % 2)
    {
        R = res;
        S = v;
    }
    else
    {
        R = v;
        S = res;
    }

    MUL(R, rlen, f, flen, f, flen, len, ctx);

    if (bit & exp)
    {
        MUL(S, rlen, R, rlen, f, flen, len, ctx);
        T = R;
        R = S;
        S = T;
    }

    while (bit >>= 1)
    {
        if (bit & exp)
        {
            MUL(S, rlen, R, rlen, R, rlen, len, ctx);
            MUL(R, rlen, S, rlen, f, flen, len, ctx);
        }
        else
        {
            MUL(S, rlen, R, rlen, R, rlen, len, ctx);
            T = R;
            R = S;
            S = T;
        }
    }

    GR_TMP_CLEAR_VEC(v, len, ctx);

    return status;
}

int
gr_poly_pow_series_ui_binexp(gr_poly_t res,
    const gr_poly_t poly, ulong exp, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong flen, rlen;

    flen = poly->length;
    len = FLINT_MAX(len, 0);

    if (exp == 0 && len != 0)
    {
        return gr_poly_one(res, ctx);
    }
    else if (flen == 0 || len == 0)
    {
        return gr_poly_zero(res, ctx);
    }
    else
    {
        rlen = poly_pow_length(flen, exp, len);

        if (res != poly)
        {
            gr_poly_fit_length(res, rlen, ctx);
            status |= _gr_poly_pow_series_ui_binexp(res->coeffs,
                poly->coeffs, flen, exp, rlen, ctx);
            _gr_poly_set_length(res, rlen, ctx);
            _gr_poly_normalise(res, ctx);
        }
        else
        {
            gr_poly_t t;
            gr_poly_init2(t, rlen, ctx);
            status |= _gr_poly_pow_series_ui_binexp(t->coeffs,
                poly->coeffs, flen, exp, rlen, ctx);
            _gr_poly_set_length(t, rlen, ctx);
            _gr_poly_normalise(t, ctx);
            gr_poly_swap(res, t, ctx);
            gr_poly_clear(t, ctx);
        }

        return status;
    }
}
