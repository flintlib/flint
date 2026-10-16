/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

static int
_gr_poly_pseudo_rem_cohen(gr_ptr R, gr_srcptr A, slong lenA,
                                gr_srcptr B, slong lenB, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    gr_srcptr leadB = GR_ENTRY(B, lenB - 1, sz);
    slong e;

    if (lenB == 1)
        return _gr_vec_zero(R, lenA, ctx);

    if (R != A)
        status |= _gr_vec_set(R, A, lenA, ctx);

    e = lenA - lenB + 1;

    while (lenA >= lenB)
    {
        status |= _gr_vec_mul_scalar(R, R, lenA - 1, leadB, ctx);
        status |= _gr_vec_submul_scalar(GR_ENTRY(R, lenA - lenB, sz),
                                        B, lenB - 1,
                                        GR_ENTRY(R, lenA - 1, sz), ctx);
        status |= gr_zero(GR_ENTRY(R, lenA - 1, sz), ctx);
        lenA--;
        status |= _gr_vec_normalise(&lenA, R, lenA, ctx);

        if (status != GR_SUCCESS)
            return status;

        e--;
    }

    if (e == 1)
    {
        status |= _gr_vec_mul_scalar(R, R, lenA, leadB, ctx);
    }
    else if (e > 1)
    {
        gr_ptr pow;
        GR_TMP_INIT(pow, ctx);
        status |= gr_pow_ui(pow, leadB, (ulong) e, ctx);
        status |= _gr_vec_mul_scalar(R, R, lenA, pow, ctx);
        GR_TMP_CLEAR(pow, ctx);
    }

    return status;
}

static int
_gr_vec_content(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i, sz = ctx->sizeof_elem;
    gr_method_binary_op op_gcd = GR_BINARY_OP(ctx, GCD);

    if (len == 0)
        return gr_zero(res, ctx);
    if (len == 1)
        return gr_set(res, vec, ctx);

    status |= op_gcd(res, GR_ENTRY(vec, 0, sz), GR_ENTRY(vec, 1, sz), ctx);
    for (i = 2; i < len && status == GR_SUCCESS; i++)
        status |= op_gcd(res, res, GR_ENTRY(vec, i, sz), ctx);

    return status;
}

/* Adapted from _fmpz_poly_resultant_euclidean */
int
_gr_poly_resultant_subresultant(gr_ptr res, gr_srcptr poly1, slong len1,
                                gr_srcptr poly2, slong len2, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (len2 == 1)
        return _gr_poly_resultant_small(res, poly1, len1, poly2, len2, ctx);

    {
        gr_ptr W, A, B, a, b, g_sc, h, t;
        slong lenA, lenB;
        int sgn;   /* tracks accumulated sign: 0 = positive, 1 = negative */

        GR_TMP_INIT_VEC(W, len1 + len2 + 5, ctx);
        A    = W;
        B    = GR_ENTRY(W, len1,         sz);
        a    = GR_ENTRY(W, len1 + len2,     sz);
        b    = GR_ENTRY(W, len1 + len2 + 1, sz);
        g_sc = GR_ENTRY(W, len1 + len2 + 2, sz);
        h    = GR_ENTRY(W, len1 + len2 + 3, sz);
        t    = GR_ENTRY(W, len1 + len2 + 4, sz);

        lenA = len1;
        lenB = len2;
        sgn  = 0;

        status |= _gr_vec_content(a, poly1, lenA, ctx);
        status |= _gr_vec_content(b, poly2, lenB, ctx);
        status |= _gr_vec_divexact_scalar(A, poly1, lenA, a, ctx);
        status |= _gr_vec_divexact_scalar(B, poly2, lenB, b, ctx);

        status |= gr_pow_ui(a, a, (ulong)(lenB - 1), ctx);
        status |= gr_pow_ui(b, b, (ulong)(lenA - 1), ctx);
        status |= gr_mul(t, a, b, ctx);

        if (status != GR_SUCCESS)
            goto cleanup;

        status |= gr_one(g_sc, ctx);
        status |= gr_one(h,    ctx);

        do
        {
            const slong d = lenA - lenB;

            if (!(lenA & WORD(1)) && !(lenB & WORD(1)))
                sgn ^= 1;

            status |= _gr_poly_pseudo_rem_cohen(A, A, lenA, B, lenB, ctx);
            status |= _gr_vec_normalise(&lenA, A, lenA, ctx);

            if (status != GR_SUCCESS)
                goto cleanup;

            if (lenA == 0)
            {
                status |= gr_zero(res, ctx);
                goto cleanup;
            }

            FLINT_SWAP(gr_ptr, A, B);
            FLINT_SWAP(slong, lenA, lenB);

            status |= gr_pow_ui(a, h, (ulong) d, ctx);   /* a = h^d */
            status |= gr_mul(b, g_sc, a, ctx);             /* b = g_sc * h^d */
            status |= _gr_vec_divexact_scalar(B, B, lenB, b, ctx);

            status |= gr_pow_ui(b, GR_ENTRY(A, lenA - 1, sz), (ulong) d, ctx);
            status |= gr_mul(h, h, b, ctx);
            status |= gr_divexact(h, h, a, ctx);
            status |= gr_set(g_sc, GR_ENTRY(A, lenA - 1, sz), ctx);

            if (status != GR_SUCCESS)
                goto cleanup;
        }
        while (lenB > 1);

        {
            const slong e = lenA - 1;
            status |= gr_pow_ui(a, h, (ulong) e, ctx);          /* a = h^e */
            status |= gr_pow_ui(b, B, (ulong) e, ctx);          /* b = B^e */
            status |= gr_mul(h, h, b, ctx);                       /* h = h * b */
            status |= gr_divexact(h, h, a, ctx);                  /* h = h / h^e */
        }

        status |= gr_mul(res, t, h, ctx);
        if (sgn)
            status |= gr_neg(res, res, ctx);

cleanup:
        GR_TMP_CLEAR_VEC(W, len1 + len2 + 5, ctx);
    }

    return status;
}

int
gr_poly_resultant_subresultant(gr_ptr r, const gr_poly_t f,
                               const gr_poly_t g, gr_ctx_t ctx)
{
    slong len1 = f->length;
    slong len2 = g->length;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (len1 == 0 || len2 == 0)
        return gr_zero(r, ctx);

    if (gr_is_zero(GR_ENTRY(f->coeffs, len1 - 1, sz), ctx) != T_FALSE ||
        gr_is_zero(GR_ENTRY(g->coeffs, len2 - 1, sz), ctx) != T_FALSE)
        return GR_UNABLE;

    if (len1 >= len2)
    {
        status |= _gr_poly_resultant_subresultant(r, f->coeffs, len1,
                                                     g->coeffs, len2, ctx);
    }
    else
    {
        status |= _gr_poly_resultant_subresultant(r, g->coeffs, len2,
                                                     f->coeffs, len1, ctx);
        if (!(len1 & WORD(1)) && !(len2 & WORD(1)))
            status |= gr_neg(r, r, ctx);
    }

    return status;
}
