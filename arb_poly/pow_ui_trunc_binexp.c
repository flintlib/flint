/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

#define MUL(z, zlen, x, xlen, y, ylen, trunc, prec) \
    do { \
        slong slen = FLINT_MIN(xlen + ylen - 1, trunc); \
        _arb_poly_mullow(z, x, xlen, y, ylen, slen, prec); \
        zlen = slen; \
    } while (0)

void
_arb_poly_pow_ui_trunc_binexp(arb_ptr res,
    arb_srcptr f, slong flen, ulong exp, slong len, slong prec)
{
    arb_ptr v, R, S, T;
    slong rlen;
    ulong bit;

    if (exp <= 1)
    {
        if (exp == 0)
            arb_one(res);
        else if (exp == 1)
            _arb_vec_set_round(res, f, len, prec);
        return;
    }

    /* (f * x^r)^m = x^(rm) * f^m */
    while (flen > 1 && arb_is_zero(f))
    {
        if (((ulong) len) > exp)
        {
            _arb_vec_zero(res, exp);
            len -= exp;
            res += exp;
        }
        else
        {
            _arb_vec_zero(res, len);
            return;
        }

        f++;
        flen--;
    }

    if (exp == 2)
    {
        _arb_poly_mullow(res, f, flen, f, flen, len, prec);
        return;
    }

    if (flen == 1)
    {
        arb_pow_ui(res, f, exp, prec);
        return;
    }

    v = _arb_vec_init(len);
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

    MUL(R, rlen, f, flen, f, flen, len, prec);

    if (bit & exp)
    {
        MUL(S, rlen, R, rlen, f, flen, len, prec);
        T = R;
        R = S;
        S = T;
    }
    
    while (bit >>= 1)
    {
        if (bit & exp)
        {
            MUL(S, rlen, R, rlen, R, rlen, len, prec);
            MUL(R, rlen, S, rlen, f, flen, len, prec);
        }
        else
        {
            MUL(S, rlen, R, rlen, R, rlen, len, prec);
            T = R;
            R = S;
            S = T;
        }
    }
    
    _arb_vec_clear(v, len);
}

void
arb_poly_pow_ui_trunc_binexp(arb_poly_t res,
    const arb_poly_t poly, ulong exp, slong len, slong prec)
{
    slong flen, rlen;

    flen = poly->length;

    if (exp == 0 && len != 0)
    {
        arb_poly_one(res);
    }
    else if (flen == 0 || len == 0)
    {
        arb_poly_zero(res);
    }
    else
    {
        rlen = poly_pow_length(flen, exp, len);

        if (res != poly)
        {
            arb_poly_fit_length(res, rlen);
            _arb_poly_pow_ui_trunc_binexp(res->coeffs,
                poly->coeffs, flen, exp, rlen, prec);
            _arb_poly_set_length(res, rlen);
            _arb_poly_normalise(res);
        }
        else
        {
            arb_poly_t t;
            arb_poly_init2(t, rlen);
            _arb_poly_pow_ui_trunc_binexp(t->coeffs,
                poly->coeffs, flen, exp, rlen, prec);
            _arb_poly_set_length(t, rlen);
            _arb_poly_normalise(t);
            arb_poly_swap(res, t);
            arb_poly_clear(t);
        }
    }
}

