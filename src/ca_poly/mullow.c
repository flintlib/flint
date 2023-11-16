/*
    Copyright (C) 2020, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

ca_field_ptr
_ca_vec_same_field2(ca_srcptr A, slong Alen, ca_srcptr B, slong Blen, ca_ctx_t ctx)
{
    ca_field_ptr K;
    ca_field_ptr QQ;
    slong i;

    QQ = ctx->field_qq;
    K = QQ;

    for (i = 0; i < Alen; i++)
    {
        if (CA_IS_QQ(A + i, ctx))
            continue;

        if (CA_IS_SPECIAL(A + i))
            return NULL;

        if (K == QQ)
            K = CA_FIELD(A + i, ctx);
        else if (K != CA_FIELD(A + i, ctx))
            return NULL;
    }

    if (B != NULL)
    {
        for (i = 0; i < Blen; i++)
        {
            if (CA_IS_QQ(B + i, ctx))
                continue;

            if (CA_IS_SPECIAL(B + i))
                return NULL;

            if (K == QQ)
                K = CA_FIELD(B + i, ctx);
            else if (K != CA_FIELD(B + i, ctx))
                return NULL;
        }
    }

    return K;
}

static void
_ca_addmul(ca_t x, ca_t tmp, const ca_t a, const ca_t b, ca_ctx_t ctx)
{
    ca_mul(tmp, a, b, ctx);
    ca_add(x, x, tmp, ctx);
}

static void
_ca_poly_sqrlow_classical(ca_ptr res, ca_srcptr poly1, slong len1,
    slong n, ca_ctx_t ctx)
{
    slong i, start, stop;
    ca_t t;

    /* Basecase squaring */
    ca_init(t, ctx);

    ca_sqr(res, poly1, ctx);
    ca_mul(res + 1, poly1, poly1 + 1, ctx);
    ca_mul_ui(res + 1, res + 1, 2, ctx);

    for (i = 2; i < FLINT_MIN(n, 2 * len1 - 3); i++)
    {
        start = FLINT_MAX(0, i - len1 + 1);
        stop = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

        ca_dot(res + i, NULL, 0, poly1 + start, 1,
            poly1 + i - start, -1, stop - start + 1, ctx);
        ca_mul_ui(res + i, res + i, 2, ctx);
        if (i % 2 == 0 && i / 2 < len1)
            _ca_addmul(res + i, t, poly1 + i / 2, poly1 + i / 2, ctx);
    }

    if (len1 > 2 && n >= 2 * len1 - 2)
    {
        ca_mul(res + 2 * len1 - 3, poly1 + len1 - 1, poly1 + len1 - 2, ctx);
        ca_mul_ui(res + 2 * len1 - 3, res + 2 * len1 - 3, 2, ctx);
    }

    if (n >= 2 * len1 - 1)
        ca_sqr(res + 2 * len1 - 2, poly1 + len1 - 1, ctx);

    ca_clear(t, ctx);
}

static void
_ca_poly_sqrlow_fmpqs(ca_ptr res, ca_srcptr poly1, slong len1,
    slong n, ca_ctx_t ctx)
{
    fmpz *z1, *z3;
    fmpz_t den1;

    if (_ca_vec_fmpq_vec_is_fmpz_vec(poly1, len1, ctx))
    {
        slong i;

        z1 = _fmpz_vec_init(len1 + n);
        z3 = z1 + len1;

        for (i = 0; i < len1; i++)
            z1[i] = *CA_FMPQ_NUMREF(poly1 + i);

        _fmpz_poly_sqrlow(z3, z1, len1, n);

        for (i = 0; i < n; i++)
        {
            _ca_make_fmpq(res + i, ctx);
            fmpz_one(CA_FMPQ_DENREF(res + i));
            fmpz_clear(CA_FMPQ_NUMREF(res + i));
            *CA_FMPQ_NUMREF(res + i) = z3[i];
        }

        flint_free(z1);
    }
    else
    {
        fmpz_init(den1);
        z1 = _fmpz_vec_init(len1 + n);
        z3 = z1 + len1;

        _ca_vec_fmpq_vec_get_fmpz_vec_den(z1, den1, poly1, len1, ctx);

        fmpz_mul(den1, den1, den1);
        _fmpz_poly_sqrlow(z3, z1, len1, n);
        _ca_vec_set_fmpz_vec_div_fmpz(res, z3, den1, n, ctx);

        _fmpz_vec_clear(z1, len1 + n);
        fmpz_clear(den1);
    }
}

static void
_ca_poly_mullow_fmpqs(ca_ptr res,
    ca_srcptr poly1, slong len1,
    ca_srcptr poly2, slong len2, slong n, ca_ctx_t ctx)
{
    fmpz *z1, *z2, *z3;
    fmpz_t den1, den2;

    /* Todo: handle mixed cases */
    if (_ca_vec_fmpq_vec_is_fmpz_vec(poly1, len1, ctx) &&
        _ca_vec_fmpq_vec_is_fmpz_vec(poly2, len2, ctx))
    {
        slong i;

        z1 = _fmpz_vec_init(len1 + len2 + n);
        z2 = z1 + len1;
        z3 = z2 + len2;

        for (i = 0; i < len1; i++)
            z1[i] = *CA_FMPQ_NUMREF(poly1 + i);
        for (i = 0; i < len2; i++)
            z2[i] = *CA_FMPQ_NUMREF(poly2 + i);

        if (len1 >= len2)
            _fmpz_poly_mullow(z3, z1, len1, z2, len2, n);
        else
            _fmpz_poly_mullow(z3, z2, len2, z1, len1, n);

        for (i = 0; i < n; i++)
        {
            _ca_make_fmpq(res + i, ctx);
            fmpz_one(CA_FMPQ_DENREF(res + i));
            fmpz_clear(CA_FMPQ_NUMREF(res + i));
            *CA_FMPQ_NUMREF(res + i) = z3[i];
        }

        flint_free(z1);
    }
    else
    {
        fmpz_init(den1);
        fmpz_init(den2);
        z1 = _fmpz_vec_init(len1 + len2 + n);
        z2 = z1 + len1;
        z3 = z2 + len2;

        _ca_vec_fmpq_vec_get_fmpz_vec_den(z1, den1, poly1, len1, ctx);
        _ca_vec_fmpq_vec_get_fmpz_vec_den(z2, den2, poly2, len2, ctx);

        fmpz_mul(den1, den1, den2);
        if (len1 >= len2)
            _fmpz_poly_mullow(z3, z1, len1, z2, len2, n);
        else
            _fmpz_poly_mullow(z3, z2, len2, z1, len1, n);

        _ca_vec_set_fmpz_vec_div_fmpz(res, z3, den1, n, ctx);

        _fmpz_vec_clear(z1, len1 + len2 + n);
        fmpz_clear(den1);
        fmpz_clear(den2);
    }
}


void
_ca_poly_mullow(ca_ptr res,
    ca_srcptr poly1, slong len1,
    ca_srcptr poly2, slong len2, slong n, ca_ctx_t ctx)
{
    ca_field_ptr K;
    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    if (n == 1)
    {
        ca_mul(res, poly1, poly2, ctx);
        return;
    }

    if (len1 == 1)
    {
        _ca_vec_scalar_mul_ca(res, poly2, n, poly1, ctx);
        return;
    }

    if (len2 == 1)
    {
        _ca_vec_scalar_mul_ca(res, poly1, n, poly2, ctx);
        return;
    }

    /* Squaring */
    if (poly1 == poly2 && len1 == len2)
    {
        /* Integers and rationals. */
        if (len1 >= 4 && _ca_vec_is_fmpq_vec(poly1, len1, ctx))
        {
            _ca_poly_sqrlow_fmpqs(res, poly1, len1, n, ctx);
            return;
        }

        /* Square over the same number field. */
        if (len1 >= 4)
        {
            K = _ca_vec_same_field2(poly1, len1, NULL, 0, ctx);

            if (K != NULL && CA_FIELD_IS_NF(K) && (FLINT_MIN(len1, len2) >= CA_FIELD_NF(K)->pol->length || FLINT_MIN(len1, len2) >= 10))
            {
                _ca_poly_mullow_same_nf(res, poly1, len1, poly2, len2, n, K, ctx);
                return;
            }
        }

        _ca_poly_sqrlow_classical(res, poly1, len1, n, ctx);
        return;
    }

    if (len1 >= 4 && len2 >= 4 && _ca_vec_is_fmpq_vec(poly1, len1, ctx) && _ca_vec_is_fmpq_vec(poly2, len2, ctx))
    {
        _ca_poly_mullow_fmpqs(res, poly1, len1, poly2, len2, n, ctx);
        return;
    }

    /* Multiply over the same number field. */
    if (len1 >= 4)
    {
        K = _ca_vec_same_field2(poly1, len1, poly2, len2, ctx);

        if (K != NULL && CA_FIELD_IS_NF(K) && (FLINT_MIN(len1, len2) >= CA_FIELD_NF(K)->pol->length || FLINT_MIN(len1, len2) >= 10))
        {
            _ca_poly_mullow_same_nf(res, poly1, len1, poly2, len2, n, K, ctx);
            return;
        }
    }

    /* General case */
    {
        slong i, top1, top2;

        ca_mul(res, poly1, poly2, ctx);

        for (i = 1; i < n; i++)
        {
            top1 = FLINT_MIN(len1 - 1, i);
            top2 = FLINT_MIN(len2 - 1, i);

            ca_dot(res + i, NULL, 0, poly1 + i - top2, 1,
                poly2 + top2, -1, top1 + top2 - i + 1, ctx);
        }
    }
}

void
ca_poly_mullow(ca_poly_t res, const ca_poly_t poly1,
                                            const ca_poly_t poly2,
                                                slong n, ca_ctx_t ctx)
{
    slong len_out;

    if (poly1->length == 0 || poly2->length == 0 || n == 0)
    {
        ca_poly_zero(res, ctx);
        return;
    }

    len_out = poly1->length + poly2->length - 1;
    n = FLINT_MIN(n, len_out);

    if (res == poly1 || res == poly2)
    {
        ca_poly_t t;
        ca_poly_init2(t, n, ctx);
        _ca_poly_mullow(t->coeffs, poly1->coeffs, poly1->length,
                                    poly2->coeffs, poly2->length, n, ctx);
        ca_poly_swap(res, t, ctx);
        ca_poly_clear(t, ctx);
    }
    else
    {
        ca_poly_fit_length(res, n, ctx);
        _ca_poly_mullow(res->coeffs, poly1->coeffs, poly1->length,
                                    poly2->coeffs, poly2->length, n, ctx);
    }

    _ca_poly_set_length(res, n, ctx);
    _ca_poly_normalise(res, ctx);
}
