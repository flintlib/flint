/*
    Copyright (C) 2019 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_poly.h"

/* counts zero bits in the binary representation of e */
static int
n_zerobits(mp_limb_t e)
{
    int zeros = 0;

    while (e > 1)
    {
        zeros += !(e & 1);
        e >>= 1;
    }

    return zeros;
}

static slong
poly_pow_length(slong poly_len, ulong exp, slong trunc)
{
    mp_limb_t hi, lo;
    umul_ppmm(hi, lo, poly_len - 1, exp);
    add_ssaaaa(hi, lo, hi, lo, 0, 1);
    if (hi != 0 || lo > (mp_limb_t) WORD_MAX)
        return trunc;
    return FLINT_MIN((slong) lo, trunc);
}

#define MUL(z, zden, zlen, x, xden, xlen, y, yden, ylen, trunc, prec) \
    do { \
        slong slen = FLINT_MIN(xlen + ylen - 1, trunc); \
        if (xlen >= ylen) \
            _fmpz_poly_mullow(z, x, xlen, y, ylen, slen); \
        else \
            _fmpz_poly_mullow(z, y, ylen, x, xlen, slen); \
        zlen = slen; \
        fmpz_mul(zden, xden, yden); \
        _fmpq_poly_canonicalise(z, zden, zlen); \
    } while (0)

void
_fmpq_poly_pow_trunc(fmpz * res, fmpz_t resden,
    const fmpz * f, const fmpz_t fden, slong flen, ulong exp, slong len)
{
    fmpz * v, * R, * S, * T, * Rden, * Sden;
    fmpz_t vden;
    slong rlen;
    ulong bit;

    if (exp <= 1)
    {
        if (exp == 0)
        {
            fmpz_one(res);
            fmpz_one(resden);
        }
        else if (exp == 1)
        {
            _fmpz_vec_set(res, f, len);
            fmpz_set(resden, fden);
            _fmpq_poly_canonicalise(res, resden, len);
        }
        return;
    }

    /* (f * x^r)^m = x^(rm) * f^m */
    while (flen > 1 && fmpz_is_zero(f))
    {
        if (((ulong) len) > exp)
        {
            _fmpz_vec_zero(res, exp);
            len -= exp;
            res += exp;
        }
        else
        {
            _fmpz_vec_zero(res, len);
            fmpz_one(resden);
            return;
        }

        f++;
        flen--;
    }

    if (exp == 2)
    {
        _fmpq_poly_mullow(res, resden, f, fden, flen, f, fden, flen, len);
        _fmpq_poly_canonicalise(res, resden, len);
        return;
    }

    if (flen == 1)
    {
        fmpz_set(res, f);
        fmpz_set(resden, fden);
        _fmpq_canonicalise(res, resden);
        fmpz_pow_ui(res, res, exp);
        fmpz_pow_ui(resden, resden, exp);
        return;
    }

    v = _fmpz_vec_init(len);
    fmpz_init(vden);
    bit = UWORD(1) << (FLINT_BIT_COUNT(exp) - 2);
    
    if (n_zerobits(exp) % 2)
    {
        R = res;
        Rden = resden;
        S = v;
        Sden = vden;
    }
    else
    {
        R = v;
        Rden = vden;
        S = res;
        Sden = resden;
    }

    MUL(R, Rden, rlen, f, fden, flen, f, fden, flen, len, prec);

    if (bit & exp)
    {
        MUL(S, Sden, rlen, R, Rden, rlen, f, fden, flen, len, prec);
        T = R;
        R = S;
        S = T;
        T = Rden;
        Rden = Sden;
        Sden = T;
    }
    
    while (bit >>= 1)
    {
        if (bit & exp)
        {
            MUL(S, Sden, rlen, R, Rden, rlen, R, Rden, rlen, len, prec);
            MUL(R, Rden, rlen, S, Sden, rlen, f, fden, flen, len, prec);
        }
        else
        {
            MUL(S, Sden, rlen, R, Rden, rlen, R, Rden, rlen, len, prec);
            T = R;
            R = S;
            S = T;
            T = Rden;
            Rden = Sden;
            Sden = T;
        }
    }

    _fmpz_vec_clear(v, len);
    fmpz_clear(vden);
}

void
fmpq_poly_pow_trunc(fmpq_poly_t res,
    const fmpq_poly_t poly, ulong exp, slong len)
{
    slong flen, rlen;

    flen = poly->length;

    if (exp == 0 && len != 0)
    {
        fmpq_poly_one(res);
    }
    else if (flen == 0 || len == 0)
    {
        fmpq_poly_zero(res);
    }
    else
    {
        rlen = poly_pow_length(flen, exp, len);

        if (res != poly)
        {
            fmpq_poly_fit_length(res, rlen);
            _fmpq_poly_pow_trunc(res->coeffs, res->den,
                poly->coeffs, poly->den, flen, exp, rlen);
            _fmpq_poly_set_length(res, rlen);
            _fmpq_poly_normalise(res);
        }
        else
        {
            fmpq_poly_t t;
            fmpq_poly_init2(t, rlen);
            _fmpq_poly_pow_trunc(t->coeffs, t->den,
                poly->coeffs, poly->den, flen, exp, rlen);
            _fmpq_poly_set_length(t, rlen);
            _fmpq_poly_normalise(t);
            fmpq_poly_swap(res, t);
            fmpq_poly_clear(t);
        }
    }
}

