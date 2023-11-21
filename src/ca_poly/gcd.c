/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"
#include "ca_poly.h"

int
_ca_poly_check_coprime_numerical(ca_srcptr A, slong lenA, ca_srcptr B, slong lenB, ca_ctx_t ctx)
{
    acb_t D;
    slong degA, degB, i, j;
    slong prec;
    int result;
    acb_ptr TA, TB;

    degA = lenA - 1;
    degB = lenB - 1;

    TA = _acb_vec_init(lenA);
    TB = _acb_vec_init(lenA);
    acb_init(D);

    prec = ctx->options[CA_OPT_LOW_PREC];

    for (i = 0; i <= degA; i++)
        ca_get_acb(TA + i, A + i, prec, ctx);
    for (i = 0; i <= degB; i++)
        ca_get_acb(TB + i, B + i, prec, ctx);

    if (_acb_vec_is_real(TA, lenA) && _acb_vec_is_real(TB, lenB))
    {
        arb_mat_t R;
        arb_mat_init(R, degA + degB, degA + degB);

        for (i = 0; i < degB; i++)
        {
            for (j = 0; j <= degA; j++)
            {
                if (i == 0)
                    arb_swap(acb_mat_entry(R, 0, j), acb_realref(TA + j));
                else
                    arb_set(arb_mat_entry(R, i, i + j), arb_mat_entry(R, 0, j));
            }
        }

        for (i = 0; i < degA; i++)
        {
            for (j = 0; j <= degB; j++)
            {
                if (i == 0)
                    arb_swap(arb_mat_entry(R, degB, j), acb_realref(TB + j));
                else
                    arb_set(acb_mat_entry(R, degB + i, i + j), arb_mat_entry(R, degB, j));
            }
        }

        arb_mat_det(acb_realref(D), R, prec);
        arb_mat_clear(R);
    }
    else
    {
        acb_mat_t C;
        acb_mat_init(C, degA + degB, degA + degB);

        for (i = 0; i < degB; i++)
        {
            for (j = 0; j <= degA; j++)
            {
                if (i == 0)
                    acb_swap(acb_mat_entry(C, 0, j), TA + j);
                else
                    acb_set(acb_mat_entry(C, i, i + j), acb_mat_entry(C, 0, j));
            }
        }

        for (i = 0; i < degA; i++)
        {
            for (j = 0; j <= degB; j++)
            {
                if (i == 0)
                    acb_swap(acb_mat_entry(C, degB, j), TB + j);
                else
                    acb_set(acb_mat_entry(C, degB + i, i + j), acb_mat_entry(C, degB, j));
            }
        }

        acb_mat_det(D, C, prec);
        acb_mat_clear(C);
    }

    result = !acb_contains_zero(D);

    _acb_vec_clear(TA, lenA);
    _acb_vec_clear(TB, lenB);
    acb_clear(D);

    return result;
}

/* assumes lenA >= lenB >= 1, and both A and B have nonzero leading
   coefficient */

slong
_ca_poly_gcd(ca_ptr G, ca_srcptr A, slong lenA,
                                ca_srcptr B, slong lenB, ca_ctx_t ctx)
{
    if (_ca_vec_is_fmpq_vec(A, lenA, ctx) && _ca_vec_is_fmpq_vec(B, lenB, ctx))
    {
        fmpz * zA, * zB, *zG;
        fmpz_t den;
        slong i, lenA2;

        fmpz_init(den);
        zA = _fmpz_vec_init(lenA);
        zB = _fmpz_vec_init(lenB);
        zG = _fmpz_vec_init(lenA);

        _ca_vec_fmpq_vec_get_fmpz_vec_den(zA, den, A, lenA, ctx);
        _ca_vec_fmpq_vec_get_fmpz_vec_den(zB, den, B, lenB, ctx);

        _fmpz_poly_gcd(zG, zA, lenA, zB, lenB);

        lenA2 = lenA;
        while (lenA2 > 1 && fmpz_is_zero(zG + lenA2 - 1))
            lenA2--;

        for (i = 0; i < lenA2; i++)
            ca_set_fmpz(G + i, zG + i, ctx);

        _fmpz_vec_clear(zA, lenA);
        _fmpz_vec_clear(zB, lenB);
        _fmpz_vec_clear(zG, lenA);
        fmpz_clear(den);

        return lenA2;
    }

    if (_ca_poly_check_coprime_numerical(A, lenA, B, lenB, ctx))
    {
        ca_one(G, ctx);
        return 1;
    }

    return _ca_poly_gcd_euclidean(G, A, lenA, B, lenB, ctx);
}

int
ca_poly_gcd(ca_poly_t G, const ca_poly_t A,
                        const ca_poly_t B, ca_ctx_t ctx)
{
    slong lenA = A->length, lenB = B->length, lenG;
    ca_ptr g;

    if (A->length == 0 && B->length == 0)
    {
        ca_poly_zero(G, ctx);
        return 1;
    }

    if (A->length == 0)
        return ca_poly_make_monic(G, B, ctx);

    if (B->length == 0)
        return ca_poly_make_monic(G, A, ctx);

    if (A->length < B->length)
        return ca_poly_gcd(G, B, A, ctx);

    if (ca_check_is_zero(A->coeffs + A->length - 1, ctx) != T_FALSE ||
        ca_check_is_zero(B->coeffs + B->length - 1, ctx) != T_FALSE)
    {
        return 0;
    }

    /* lenA >= lenB >= 1 */
    if (G == A || G == B)
    {
        g = _ca_vec_init(FLINT_MIN(lenA, lenB), ctx);
    }
    else
    {
        ca_poly_fit_length(G, FLINT_MIN(lenA, lenB), ctx);
        g = G->coeffs;
    }

    lenG = _ca_poly_gcd(g, A->coeffs, lenA, B->coeffs, lenB, ctx);

    if (G == A || G == B)
    {
        _ca_vec_clear(G->coeffs, G->alloc, ctx);
        G->coeffs = g;
        G->alloc = FLINT_MIN(lenA, lenB);
        G->length = FLINT_MIN(lenA, lenB);
    }
    _ca_poly_set_length(G, lenG, ctx);

    if (lenG == 0)
    {
        return 0;
    }
    else
    {
        if (G->length == 1)
            ca_one(G->coeffs, ctx);
        else
            ca_poly_make_monic(G, G, ctx);
        return 1;
    }
}
