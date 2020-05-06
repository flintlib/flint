/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_fmpz_poly.h"
#include "qqbar.h"

static ulong _deflation(const fmpz * poly, slong len)
{
    fmpz_poly_t t;
    t->alloc = t->length = len;
    t->coeffs = (fmpz *) poly;
    return arb_fmpz_poly_deflation(t);
}

/* todo: use linear algebra instead of naive horner */
void
_qqbar_evaluate_fmpq_poly(qqbar_t res, const fmpz * poly, const fmpz_t den, slong len, const qqbar_t x)
{
    ulong deflation;

    if (len == 0)
    {
        qqbar_zero(res);
    }
    else if (len == 1)
    {
        if (fmpz_is_one(den))
        {
            qqbar_set_fmpz(res, poly);
        }
        else
        {
            fmpq_t t;
            fmpq_init(t);
            fmpq_set_fmpz_frac(t, poly, den);
            qqbar_set_fmpq(res, t);
            fmpq_clear(t);
        }
    }
    else if (qqbar_is_rational(x))
    {
        fmpq_t t, u;
        fmpq_init(t);
        fmpq_init(u);
        fmpz_neg(fmpq_numref(u), QQBAR_COEFFS(x));
        fmpz_set(fmpq_denref(u), QQBAR_COEFFS(x) + 1);
        _fmpq_poly_evaluate_fmpq(fmpq_numref(t), fmpq_denref(t), poly, den, len, fmpq_numref(u), fmpq_denref(u));
        qqbar_set_fmpq(res, t);
        fmpq_clear(t);
        fmpq_clear(u);
    }
    else if (len == 2)
    {
        qqbar_scalar_op(res, x, poly + 1, poly, den);
    }
    else if (fmpz_is_zero(poly))
    {
        slong n;
        qqbar_t t;

        n = 1;
        while (n < len && fmpz_is_zero(poly + n))
            n++;

        qqbar_init(t);
        qqbar_pow_ui(t, x, n);
        _qqbar_evaluate_fmpq_poly(res, poly + n, den, len - n, x);
        qqbar_mul(res, res, t);
        qqbar_clear(t);
    }
    else if ((deflation =_deflation(poly, len)) > 1)
    {
        slong i, len2;
        fmpz * tmp;
        qqbar_t t;

        len2 = (len - 1) / deflation + 1;
        tmp = flint_malloc(sizeof(fmpz) * len2);
        for (i = 0; i < len2; i++)
            tmp[i] = poly[i * deflation];

        qqbar_init(t);
        qqbar_pow_ui(t, x, deflation);
        _qqbar_evaluate_fmpq_poly(res, tmp, den, len2, t);
        qqbar_clear(t);
        flint_free(tmp);
    }
    else if (len > qqbar_degree(x))
    {
        fmpz * tmp;
        fmpz_t r, one;
        slong d, len2;

        d = qqbar_degree(x);
        tmp = _fmpz_vec_init(len);
        fmpz_init(r);
        fmpz_init(one);
        fmpz_one(one);

        _fmpq_poly_rem(tmp, r, poly, den, len,
            QQBAR_COEFFS(x), one, d + 1, NULL);

        len2 = d;
        while (len2 > 1 && fmpz_is_zero(tmp + len2 - 1))
            len2--;

        _qqbar_evaluate_fmpq_poly(res, tmp, r, len2, x);

        fmpz_clear(r);
        fmpz_clear(one);
        _fmpz_vec_clear(tmp, d);
    }
    else
    {
        qqbar_t t;
        slong i;

        qqbar_init(t);

        qqbar_mul_fmpz(t, x, poly + len - 1);
        qqbar_add_fmpz(t, t, poly + len - 2);

        for (i = len - 3; i >= 0; i--)
        {
            qqbar_mul(t, t, x);
            qqbar_add_fmpz(t, t, poly + i);
        }

        qqbar_div_fmpz(res, t, den);

        qqbar_clear(t);
    }
}

void
qqbar_evaluate_fmpq_poly(qqbar_t res, const fmpq_poly_t poly, const qqbar_t x)
{
    _qqbar_evaluate_fmpq_poly(res, poly->coeffs, poly->den, poly->length, x);
}

