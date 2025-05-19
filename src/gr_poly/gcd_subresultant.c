/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

/* todo: make this public and tested */
int
_gr_poly_pseudo_rem_cohen(gr_ptr R, gr_srcptr A, slong lenA,
                            gr_srcptr B, slong lenB, gr_ctx_t ctx)
{
    slong e;
    gr_ptr pow;
    slong sz = ctx->sizeof_elem; 
    int status = GR_SUCCESS;
    gr_srcptr leadB = GR_ENTRY(B, lenB - 1, sz);

    if (lenB == 1)
        return _gr_vec_zero(R, lenA, ctx);

    if (R != A)
        status |= _gr_vec_set(R, A, lenA, ctx);

    e = lenA - lenB + 1;

    while (lenA >= lenB)
    {
        status |= _gr_vec_mul_scalar(R, R, lenA - 1, leadB, ctx);
        status |= _gr_vec_submul_scalar(GR_ENTRY(R, lenA - lenB, sz), B, lenB - 1, GR_ENTRY(R, lenA - 1, sz), ctx);
        status |= gr_zero(GR_ENTRY(R, lenA - 1, sz), ctx);
        lenA--;
        status |= _gr_vec_normalise(&lenA, R, lenA, ctx);

        if (status != GR_SUCCESS)
        {
            /* todo: do we want to zero R? */
            goto cleanup;
        }

        e--;
    }
    
    if (e == 1)
        status |= _gr_vec_mul_scalar(R, R, lenA, leadB, ctx);
    else if (e != 0)
    {
        GR_TMP_INIT(pow, ctx);
        status |= gr_pow_ui(pow, leadB, e, ctx);
        status |= _gr_vec_mul_scalar(R, R, lenA, pow, ctx);
        GR_TMP_CLEAR(pow, ctx);
    }

cleanup:
    return status;
}

/* todo: check leading coeff? */
int
gr_poly_pseudo_rem_cohen(gr_poly_t R, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
{
    slong lenr;
    int status;

    if (B->length == 0)
        return GR_DOMAIN;

    if (A->length < B->length)
        return gr_poly_set(R, A, ctx);

    lenr = A->length;

    if (R == B)
    {
        gr_poly_t T;
        gr_poly_init2(T, lenr, ctx);
        status = _gr_poly_pseudo_rem_cohen(T->coeffs, A->coeffs, A->length, B->coeffs, B->length, ctx);
        gr_poly_swap(T, R, ctx);
        gr_poly_clear(T, ctx);
    }
    else
    {
        status = _gr_poly_pseudo_rem_cohen(R->coeffs, A->coeffs, A->length, B->coeffs, B->length, ctx);
    }

    _gr_poly_set_length(R, B->length - 1, ctx);
    _gr_poly_normalise(R, ctx);
    return status;
}

/* todo: before making this public, when len == 1, call canonical_associate */
int
_gr_vec_content(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op gcd = GR_BINARY_OP(ctx, GCD);
    int status = GR_SUCCESS;
    slong i, sz = ctx->sizeof_elem;

    if (len <= 2)
    {
        if (len == 2)
            return gcd(res, vec, GR_ENTRY(vec, 1, ctx->sizeof_elem), ctx);
        else if (len == 1)
            return gr_set(res, vec, ctx);
        else
            return gr_zero(res, ctx);
    }

    /* todo: parallel computation (as in _gr_vec_product_serial) */
    /* todo: bidirectional or randomized sampling to try to
             reach a trivial gcd faster */
    status |= gcd(res, GR_ENTRY(vec, 0, sz), GR_ENTRY(vec, 1, sz), ctx);
    for (i = 2; i < len; i++)
        status |= gcd(res, res, GR_ENTRY(vec, i, sz), ctx);

    return status;
}

int
_gr_poly_gcd_subresultant(gr_ptr res, slong * len_res, gr_srcptr poly1, slong len1,
                                        gr_srcptr poly2, slong len2, gr_ctx_t ctx)
{
    int status;

    if (len2 == 1)
    {
        if (len1 == 1)
        {
            status = gr_gcd(res, poly1, poly2, ctx);
        }
        else
        {
            gr_ptr c;
            GR_TMP_INIT(c, ctx);
            status = _gr_vec_content(c, poly1, len1, ctx);
            status |= gr_gcd(res, c, poly2, ctx);
            GR_TMP_CLEAR(c, ctx);
        }

        *len_res = 1;
    }
    else
    {
        gr_ptr a, b, d, g, h;
        gr_ptr A, B, W;
        slong lenA, lenB;
        status = GR_SUCCESS;
        slong sz = ctx->sizeof_elem;

        GR_TMP_INIT_VEC(W, len1 + len2 + 5, ctx);
        A = W;
        B = GR_ENTRY(W, len1, sz);
        a = GR_ENTRY(B, len2, sz);
        b = GR_ENTRY(a, 1, sz);
        d = GR_ENTRY(b, 1, sz);
        g = GR_ENTRY(d, 1, sz);
        h = GR_ENTRY(g, 1, sz);

        lenA = len1;
        lenB = len2;

        status |= _gr_vec_content(a, poly1, lenA, ctx);
        status |= _gr_vec_content(b, poly2, lenB, ctx);
        /* todo: check for 1? */
        status |= _gr_vec_divexact_scalar(A, poly1, lenA, a, ctx);
        status |= _gr_vec_divexact_scalar(B, poly2, lenB, b, ctx);
        status |= gr_gcd(d, a, b, ctx);
        status |= gr_one(g, ctx);
        status |= gr_one(h, ctx);

        while (1)
        {
            slong delta = lenA - lenB;

            status |= _gr_poly_pseudo_rem_cohen(A, A, lenA, B, lenB, ctx);

            /* GR_VEC_NORM(A, lenA); */
            while (lenA >= 1)
            {
                truth_t is_zero = gr_is_zero(GR_ENTRY(A, lenA - 1, sz), ctx);

                if (is_zero == T_UNKNOWN)
                {
                    status = GR_UNABLE;
                    goto cleanup;
                }

                if (is_zero == T_FALSE)
                    break;

                lenA--;
            }

            if (lenA <= 1)
                break;

            /* Swap A and B */
            FLINT_SWAP(gr_ptr, A, B);
            FLINT_SWAP(slong, lenA, lenB);

            if (delta == 1)
            {
                status |= gr_mul(b, g, h, ctx);
                status |= _gr_vec_divexact_scalar(B, B, lenB, b, ctx);
                status |= gr_set(g, GR_ENTRY(A, lenA - 1, sz), ctx);
                status |= gr_set(h, g, ctx);
            }
            else
            {
                status |= gr_pow_ui(a, h, delta, ctx);
                status |= gr_mul(b, g, a, ctx);
                status |= _gr_vec_divexact_scalar(B, B, lenB, b, ctx);
                status |= gr_pow_ui(b, GR_ENTRY(A, lenA - 1, sz), delta, ctx);
                status |= gr_mul(g, h, b, ctx);
                status |= gr_divexact(h, g, a, ctx);
                status |= gr_set(g, GR_ENTRY(A, lenA - 1, sz), ctx);
            }
        }

        if (lenA == 1)
        {
            status |= gr_set(res, d, ctx);
            *len_res = 1;
        }
        else
        {
            status |= _gr_vec_content(b, B, lenB, ctx);
            status |= _gr_vec_divexact_scalar(B, B, lenB, b, ctx);
            status |= _gr_vec_mul_scalar(res, B, lenB, d, ctx);
            *len_res = lenB;
        }

cleanup:
        GR_TMP_CLEAR_VEC(W, len1 + len2 + 5, ctx);
    }

    if (status != GR_SUCCESS)
        *len_res = 0;

    return status;
}

int
gr_poly_gcd_subresultant(gr_poly_t G, const gr_poly_t A,
                        const gr_poly_t B, gr_ctx_t ctx)
{
    return gr_poly_gcd_wrapper((gr_method_poly_gcd_op) _gr_poly_gcd_subresultant, 1, G, A, B, ctx);
}

