/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_fmpz_poly.h"
#include "ca_qqbar.h"

static ulong _deflation(const fmpz * poly, slong len)
{
    fmpz_poly_t t;
    t->alloc = t->length = len;
    t->coeffs = (fmpz *) poly;
    return arb_fmpz_poly_deflation(t);
}

/* todo: use linear algebra instead of naive horner */
void
_ca_qqbar_evaluate_fmpq_poly(ca_qqbar_t res, const fmpz * poly, const fmpz_t den, slong len, const ca_qqbar_t x)
{
    ulong deflation;

    if (len == 0)
    {
        ca_qqbar_zero(res);
    }
    else if (len == 1)
    {
        if (fmpz_is_one(den))
        {
            ca_qqbar_set_fmpz(res, poly);
        }
        else
        {
            fmpq_t t;
            fmpq_init(t);
            fmpq_set_fmpz_frac(t, poly, den);
            ca_qqbar_set_fmpq(res, t);
            fmpq_clear(t);
        }
    }
    else if (ca_qqbar_is_rational(x))
    {
        fmpq_t t, u;
        fmpq_init(t);
        fmpq_init(u);
        fmpz_neg(fmpq_numref(u), CA_QQBAR_COEFFS(x));
        fmpz_set(fmpq_denref(u), CA_QQBAR_COEFFS(x) + 1);
        _fmpq_poly_evaluate_fmpq(fmpq_numref(t), fmpq_denref(t), poly, den, len, fmpq_numref(u), fmpq_denref(u));
        ca_qqbar_set_fmpq(res, t);
        fmpq_clear(t);
        fmpq_clear(u);
    }
    else if (len == 2)
    {
        ca_qqbar_scalar_op(res, x, poly + 1, poly, den);
    }
    else if (fmpz_is_zero(poly))
    {
        slong n;
        ca_qqbar_t t;

        n = 1;
        while (n < len && fmpz_is_zero(poly + n))
            n++;

        ca_qqbar_init(t);
        ca_qqbar_pow_ui(t, x, n);
        _ca_qqbar_evaluate_fmpq_poly(res, poly + n, den, len - n, x);
        ca_qqbar_mul(res, res, t);
        ca_qqbar_clear(t);
    }
    else if ((deflation =_deflation(poly, len)) > 1)
    {
        slong i, len2;
        fmpz * tmp;
        ca_qqbar_t t;

        len2 = (len - 1) / deflation + 1;
        tmp = flint_malloc(sizeof(fmpz) * len2);
        for (i = 0; i < len2; i++)
            tmp[i] = poly[i * deflation];

        ca_qqbar_init(t);
        ca_qqbar_pow_ui(t, x, deflation);
        _ca_qqbar_evaluate_fmpq_poly(res, tmp, den, len2, t);
        ca_qqbar_clear(t);
        flint_free(tmp);
    }
    else if (len > ca_qqbar_degree(x))
    {
        fmpz * tmp;
        fmpz_t r, one;
        slong d, len2;

        d = ca_qqbar_degree(x);
        tmp = _fmpz_vec_init(len);
        fmpz_init(r);
        fmpz_init(one);
        fmpz_one(one);

        _fmpq_poly_rem(tmp, r, poly, den, len,
            CA_QQBAR_COEFFS(x), one, d + 1, NULL);

        len2 = d;
        while (len2 > 1 && fmpz_is_zero(tmp + len2 - 1))
            len2--;

        _ca_qqbar_evaluate_fmpq_poly(res, tmp, r, len2, x);

        fmpz_clear(r);
        fmpz_clear(one);
        _fmpz_vec_clear(tmp, d);
    }
    else
    {
        ca_qqbar_t t;
        slong i;

        ca_qqbar_init(t);

        ca_qqbar_mul_fmpz(t, x, poly + len - 1);
        ca_qqbar_add_fmpz(t, t, poly + len - 2);

        for (i = len - 3; i >= 0; i--)
        {
            ca_qqbar_mul(t, t, x);
            ca_qqbar_add_fmpz(t, t, poly + i);
        }

        ca_qqbar_div_fmpz(res, t, den);

        ca_qqbar_clear(t);
    }
}

void
ca_qqbar_evaluate_fmpq_poly(ca_qqbar_t res, const fmpq_poly_t poly, const ca_qqbar_t x)
{
    _ca_qqbar_evaluate_fmpq_poly(res, poly->coeffs, poly->den, poly->length, x);
}

