/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_divrem_divconquer_recursive(gr_ptr Q, gr_ptr BQ, gr_ptr W,
    gr_srcptr A, gr_srcptr B, slong lenB, gr_srcptr invB, slong cutoff, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (lenB < FLINT_MAX(2, cutoff))
    {
        status |= _gr_vec_zero(BQ, lenB - 1, ctx);
        status |= _gr_vec_set(GR_ENTRY(BQ, lenB - 1, sz), GR_ENTRY(A, lenB - 1, sz), lenB, ctx);

        if (invB == NULL)
            status |= _gr_poly_divrem_basecase_noinv(Q, BQ, BQ, 2 * lenB - 1, B, lenB, ctx);
        else
            status |= _gr_poly_divrem_basecase_preinv1(Q, BQ, BQ, 2 * lenB - 1, B, lenB, invB, ctx);

        status |= _gr_vec_neg(BQ, BQ, lenB - 1, ctx);
        status |= _gr_vec_set(GR_ENTRY(BQ, lenB - 1, sz), GR_ENTRY(A, lenB - 1, sz), lenB, ctx);
    }
    else
    {
        const slong n2 = lenB / 2;
        const slong n1 = lenB - n2;

        gr_ptr W1 = W;
        gr_ptr W2 = GR_ENTRY(W, lenB, sz);

        gr_srcptr p1 = GR_ENTRY(A, 2 * n2, sz);
        gr_srcptr p2;
        gr_srcptr d1 = GR_ENTRY(B, n2, sz);
        gr_srcptr d2 = B;
        gr_srcptr d3 = GR_ENTRY(B, n1, sz);
        gr_srcptr d4 = B;

        gr_ptr q1 = GR_ENTRY(Q, n2, sz);
        gr_ptr q2 = Q;
        gr_ptr dq1 = GR_ENTRY(BQ, n2, sz);
        gr_ptr d1q1 = GR_ENTRY(BQ, 2 * n2, sz);

        gr_ptr d2q1, d3q2, d4q2, t;

        /*
           Set q1 to p1 div d1, a 2 n1 - 1 by n1 division so q1 ends up
           being of length n1;  d1q1 = d1 q1 is of length 2 n1 - 1
         */

        status |= _gr_poly_divrem_divconquer_recursive(q1, d1q1, W1, p1, d1, n1, invB, cutoff, ctx);

        /*
           Compute d2q1 = d2 q1, of length lenB - 1
         */

        d2q1 = W1;
        status |= _gr_poly_mul(d2q1, q1, n1, d2, n2, ctx);

        /*
           Compute dq1 = d1 q1 x^n2 + d2 q1, of length 2 n1 + n2 - 1
         */

        _gr_vec_swap(dq1, d2q1, n2, ctx);
        status |= _gr_poly_add(GR_ENTRY(dq1, n2, sz), GR_ENTRY(dq1, n2, sz), n1 - 1, GR_ENTRY(d2q1, n2, sz), n1 - 1, ctx);

        /*
           Compute t = A/x^n2 - dq1, which has length 2 n1 + n2 - 1, but we
           are not interested in the top n1 coeffs as they will be zero, so
           this has effective length n1 + n2 - 1
           For the following division, we want to set {p2, 2 n2 - 1} to the
           top 2 n2 - 1 coeffs of this
           Since the bottom n2 - 1 coeffs of p2 are irrelevant for the
           division, we in fact set {t, n2} to the relevant coeffs
         */

        t = BQ;
        status |= _gr_poly_sub(t, GR_ENTRY(A, n2 + (n1 - 1), sz), n2, GR_ENTRY(dq1, (n1 - 1), sz), n2, ctx);
        p2 = GR_ENTRY(t, - (n2 - 1), sz);

        /*
           Compute q2 = t div d3, a 2 n2 - 1 by n2 division, so q2 will have
           length n2; let d3q2 = d3 q2, of length 2 n2 - 1
         */

        d3q2 = W1;
        status |= _gr_poly_divrem_divconquer_recursive(q2, d3q2, W2, p2, d3, n2, invB, cutoff, ctx);

        /*
           Compute d4q2 = d4 q2, of length n1 + n2 - 1 = lenB - 1
         */

        d4q2 = W2;
        status |= _gr_poly_mul(d4q2, d4, n1, q2, n2, ctx);

        /*
           Compute dq2 = d3q2 x^n1 + d4q2, of length n1 + 2 n2 - 1
         */

        _gr_vec_swap(BQ, d4q2, n2, ctx);
        status |= _gr_poly_add(GR_ENTRY(BQ, n2, sz), GR_ENTRY(BQ, n2, sz), n1 - 1, GR_ENTRY(d4q2, n2, sz), n1 - 1, ctx);
        status |= _gr_poly_add(GR_ENTRY(BQ, n1, sz), GR_ENTRY(BQ, n1, sz), 2 * n2 - 1, d3q2, 2 * n2 - 1, ctx);

        /*
           Note Q = q1 x^n2 + q2, and BQ = dq1 x^n2 + dq2
         */
    }

    return status;
}

static int
__gr_poly_divrem_divconquer(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr invB, slong cutoff, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (lenA < 2 * lenB - 1)
    {
        /*
           Convert unbalanced division into a 2 n1 - 1 by n1 division
         */

        const slong n1 = lenA - lenB + 1;
        const slong n2 = lenB - n1;

        gr_srcptr p1 = GR_ENTRY(A, n2, sz);
        gr_srcptr d1 = GR_ENTRY(B, n2, sz);
        gr_srcptr d2 = B;
        gr_ptr W, d1q1, d2q1;

        GR_TMP_INIT_VEC(W, (2 * n1 - 1) + lenB - 1, ctx);

        d1q1 = GR_ENTRY(R, n2, sz);
        d2q1 = GR_ENTRY(W, (2 * n1 - 1), sz);

        status |= _gr_poly_divrem_divconquer_recursive(Q, d1q1, W, p1, d1, n1, invB, cutoff, ctx);

        /*
           Compute d2q1 = Q d2, of length lenB - 1
         */

        if (n1 >= n2)
            status |= _gr_poly_mul(d2q1, Q, n1, d2, n2, ctx);
        else
            status |= _gr_poly_mul(d2q1, d2, n2, Q, n1, ctx);

        /*
           Compute BQ = d1q1 * x^n1 + d2q1, of length lenB - 1;
           then compute R = A - BQ
         */

        _gr_vec_swap(R, d2q1, n2, ctx);
        status |= _gr_poly_add(GR_ENTRY(R, n2, sz), GR_ENTRY(R, n2, sz), n1 - 1, GR_ENTRY(d2q1, n2, sz), n1 - 1, ctx);
        status |= _gr_poly_sub(R, A, lenA, R, lenA, ctx);

        GR_TMP_CLEAR_VEC(W, (2 * n1 - 1) + lenB - 1, ctx);
    }
    else                        /* lenA = 2 * lenB - 1 */
    {
        gr_ptr W;

        GR_TMP_INIT_VEC(W, lenA, ctx);

        status |= _gr_poly_divrem_divconquer_recursive(Q, R, W, A, B, lenB, invB, cutoff, ctx);
        status |= _gr_poly_sub(R, A, lenB - 1, R, lenB - 1, ctx);

        GR_TMP_CLEAR_VEC(W, lenA, ctx);
    }

    return status;
}

int
_gr_poly_divrem_divconquer_preinv1(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr invB, slong cutoff, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (lenA <= 2 * lenB - 1)
    {
        gr_ptr W;
        /* todo: the fq version avoids this temporary, writing directly to R */
        GR_TMP_INIT_VEC(W, lenA, ctx);
        status |= __gr_poly_divrem_divconquer(Q, W, A, lenA, B, lenB, invB, cutoff, ctx);
        _gr_vec_swap(R, W, lenB - 1, ctx);
        GR_TMP_CLEAR_VEC(W, lenA, ctx);
    }
    else                        /* lenA > 2 * lenB - 1 */
    {
        slong shift, n = 2 * lenB - 1, len1;
        gr_ptr QB, W, S;

        len1 = 2 * n + lenA;
        GR_TMP_INIT_VEC(W, len1, ctx);
        S = GR_ENTRY(W, 2 * n, sz);
        status |= _gr_vec_set(S, A, lenA, ctx);
        QB = GR_ENTRY(W, n, sz);

        while (lenA >= n)
        {
            shift = lenA - n;
            status |= _gr_poly_divrem_divconquer_recursive(GR_ENTRY(Q, shift, sz), QB, W, GR_ENTRY(S, shift, sz), B, lenB, invB, cutoff, ctx);
            status |= _gr_poly_sub(GR_ENTRY(S, shift, sz), GR_ENTRY(S, shift, sz), n, QB, n, ctx);
            lenA -= lenB;
        }

        if (lenA >= lenB)
        {
            status |= __gr_poly_divrem_divconquer(Q, W, S, lenA, B, lenB, invB, cutoff, ctx);
            _gr_vec_swap(W, S, lenA, ctx);
        }

        _gr_vec_swap(R, S, lenB - 1, ctx);
        GR_TMP_CLEAR_VEC(W, len1, ctx);
    }

    return status;
}

int
_gr_poly_divrem_divconquer_noinv(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong cutoff, gr_ctx_t ctx)
{
    return _gr_poly_divrem_divconquer_preinv1(Q, R, A, lenA, B, lenB, NULL, cutoff, ctx);
}

int
_gr_poly_divrem_divconquer(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong cutoff, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    gr_ptr invB;

    GR_TMP_INIT(invB, ctx);

    status = gr_inv(invB, GR_ENTRY(B, lenB - 1, sz), ctx);

    if (status == GR_SUCCESS)
        status = _gr_poly_divrem_divconquer_preinv1(Q, R, A, lenA, B, lenB, invB, cutoff, ctx);
    else
        status = _gr_poly_divrem_divconquer_noinv(Q, R, A, lenA, B, lenB, cutoff, ctx);

    GR_TMP_CLEAR(invB, ctx);

    return status;
}

int
gr_poly_divrem_divconquer(gr_poly_t Q, gr_poly_t R,
    const gr_poly_t A, const gr_poly_t B, slong cutoff, gr_ctx_t ctx)
{
    slong lenA = A->length, lenB = B->length, lenQ = lenA - lenB + 1;
    slong sz = ctx->sizeof_elem;
    gr_poly_t tQ, tR;
    gr_ptr q, r;
    int status = GR_SUCCESS;

    if (lenB == 0)
        return GR_DOMAIN;

    if (gr_is_zero(GR_ENTRY(B->coeffs, lenB - 1, sz), ctx) != T_FALSE)
        return GR_UNABLE;

    if (lenA < lenB)
    {
        status |= gr_poly_set(R, A, ctx);
        status |= gr_poly_zero(Q, ctx);
        return status;
    }

    if (Q == A || Q == B)
    {
        gr_poly_init2(tQ, lenQ, ctx);
        q = tQ->coeffs;
    }
    else
    {
        gr_poly_fit_length(Q, lenQ, ctx);
        q = Q->coeffs;
    }

    if (R == B)
    {
        gr_poly_init2(tR, lenB - 1, ctx);
        r = tR->coeffs;
    }
    else
    {
        gr_poly_fit_length(R, lenB - 1, ctx);
        r = R->coeffs;
    }

    status |= _gr_poly_divrem_divconquer(q, r, A->coeffs, lenA, B->coeffs, lenB, cutoff, ctx);

    if (Q == A || Q == B)
    {
        gr_poly_swap(tQ, Q, ctx);
        gr_poly_clear(tQ, ctx);
    }

    if (R == B)
    {
        gr_poly_swap(tR, R, ctx);
        gr_poly_clear(tR, ctx);
    }

    _gr_poly_set_length(Q, lenQ, ctx);
    _gr_poly_set_length(R, lenB - 1, ctx);
    _gr_poly_normalise(R, ctx);

    return status;
}
