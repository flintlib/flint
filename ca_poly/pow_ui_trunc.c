/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

#define MUL(z, zlen, x, xlen, y, ylen, trunc, ctx) \
    do { \
        slong slen = FLINT_MIN(xlen + ylen - 1, trunc); \
        _ca_poly_mullow(z, x, xlen, y, ylen, slen, ctx); \
        zlen = slen; \
    } while (0)

void
_ca_poly_pow_ui_trunc(ca_ptr res,
    ca_srcptr f, slong flen, ulong exp, slong len, ca_ctx_t ctx)
{
    ca_ptr v, R, S, T;
    slong rlen;
    ulong bit;

    if (exp <= 1)
    {
        if (exp == 0)
            ca_one(res, ctx);
        else if (exp == 1)
            _ca_vec_set(res, f, len, ctx);
        return;
    }

    /* (f * x^r)^m = x^(rm) * f^m */
    while (flen > 1 && (ca_check_is_zero(f, ctx) == T_TRUE))
    {
        if (((ulong) len) > exp)
        {
            _ca_vec_zero(res, exp, ctx);
            len -= exp;
            res += exp;
        }
        else
        {
            _ca_vec_zero(res, len, ctx);
            return;
        }

        f++;
        flen--;
    }

    if (exp == 2)
    {
        _ca_poly_mullow(res, f, flen, f, flen, len, ctx);
        return;
    }

    if (flen == 1)
    {
        ca_pow_ui(res, f, exp, ctx);
        return;
    }

    v = _ca_vec_init(len, ctx);
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
    
    _ca_vec_clear(v, len, ctx);
}

void
ca_poly_pow_ui_trunc(ca_poly_t res,
    const ca_poly_t poly, ulong exp, slong len, ca_ctx_t ctx)
{
    slong flen, rlen;

    flen = poly->length;

    if (exp == 0 && len != 0)
    {
        ca_poly_one(res, ctx);
    }
    else if (flen == 0 || len == 0)
    {
        ca_poly_zero(res, ctx);
    }
    else
    {
        rlen = poly_pow_length(flen, exp, len);

        if (res != poly)
        {
            ca_poly_fit_length(res, rlen, ctx);
            _ca_poly_pow_ui_trunc(res->coeffs,
                poly->coeffs, flen, exp, rlen, ctx);
            _ca_poly_set_length(res, rlen, ctx);
            _ca_poly_normalise(res, ctx);
        }
        else
        {
            ca_poly_t t;
            ca_poly_init2(t, rlen, ctx);
            _ca_poly_pow_ui_trunc(t->coeffs,
                poly->coeffs, flen, exp, rlen, ctx);
            _ca_poly_set_length(t, rlen, ctx);
            _ca_poly_normalise(t, ctx);
            ca_poly_swap(res, t, ctx);
            ca_poly_clear(t, ctx);
        }
    }
}

