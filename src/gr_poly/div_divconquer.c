/*
    Copyright (C) 2008, 2009, 2011 William Hart
    Copyright (C) 2012 Fredrik Johansson

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
    gr_srcptr A, gr_srcptr B, slong lenB, gr_srcptr invB, slong cutoff, gr_ctx_t ctx);

int
_gr_poly_div_divconquer_recursive(gr_ptr Q, gr_ptr W,
                                gr_srcptr A, gr_srcptr B, slong lenB, gr_srcptr invB, slong cutoff, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (lenB < FLINT_MAX(2, cutoff))
    {
        if (invB == NULL)
            return _gr_poly_div_basecase_noinv(Q, A, 2 * lenB - 1, B, lenB, ctx);
        else
            return _gr_poly_div_basecase_preinv1(Q, A, 2 * lenB - 1, B, lenB, invB, ctx);
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

        gr_ptr q1   = GR_ENTRY(Q, n2, sz);
        gr_ptr q2   = Q;
        gr_ptr d1q1 = W2;

        gr_ptr d2q1, t;

        /*
           Set q1 to p1 div d1, a 2 n1 - 1 by n1 division so q1 ends up
           being of length n1;  low(d1q1) = d1 q1 is of length n1 - 1
         */

        status |= _gr_poly_divrem_divconquer_recursive(q1, d1q1, W1, p1, d1, n1, invB, cutoff, ctx);
        /*
           Compute bottom n1 + n2 - 1 coeffs of d2q1 = d2 q1
         */

        d2q1 = W1;
        status |= _gr_poly_mullow(d2q1, q1, n1, d2, n2, n1 + n2 - 1, ctx);

        /*
           Compute dq1 = d1 q1 x^n2 + d2 q1, of length n1 + n2 - 1
           Split it into a segment of length n1 - 1 at which is ignored
           and a piece of length n2 at BQ.
         */

        if (n2 > n1 - 1)
            status |= gr_set(W1, GR_ENTRY(d2q1, n1 - 1, sz), ctx);

        status |= _gr_poly_add(GR_ENTRY(W1, n2 - (n1 - 1), sz), d1q1, n1 - 1, GR_ENTRY(d2q1, n2, sz), n1 - 1, ctx);

        /*
           Compute t = A/x^n2 - dq1, which has length 2 n1 + n2 - 1, but we
           are not interested in the top n1 coeffs as they will be zero, so
           this has effective length n1 + n2 - 1
           For the following division, we want to set {p2, 2 n2 - 1} to the
           top 2 n2 - 1 coeffs of this
           Since the bottom n2 - 1 coeffs of p2 are irrelevant for the
           division, we in fact set {t, n2} to the relevant coeffs
         */

        t = W1;
        status |= _gr_poly_sub(t, GR_ENTRY(A, n2 + (n1 - 1), sz), n2, t, n2, ctx);
        p2 = GR_ENTRY(t, - (n2 - 1), sz);

        /*
           Compute q2 = t div d3, a 2 n2 - 1 by n2 division, so q2 will have
           length n2;
         */
        status |= _gr_poly_div_divconquer_recursive(q2, W2, p2, d3, n2, invB, cutoff, ctx);

        /*
           Note Q = q1 x^n2 + q2
         */

        return status;
    }
}


static int
__gr_poly_div_divconquer(gr_ptr Q, gr_srcptr A, slong lenA,
                  gr_srcptr B, slong lenB, gr_srcptr invB, slong cutoff, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (lenA < 2 * lenB - 1)
    {
        /*
           Convert unbalanced division into a 2 n1 - 1 by n1 division
         */

        const slong n1 = lenA - lenB + 1;
        const slong n2 = lenB - n1;

        gr_srcptr p1 = GR_ENTRY(A, n2, sz);
        gr_srcptr d1 = GR_ENTRY(B, n2, sz);
        gr_ptr W;

        GR_TMP_INIT_VEC(W, 2 * lenB, ctx);
        status |= _gr_poly_div_divconquer_recursive(Q, W, p1, d1, n1, invB, cutoff, ctx);
        GR_TMP_CLEAR_VEC(W, 2 * lenB, ctx);
    }
    else  /* lenA = 2 * lenB - 1 */
    {
        gr_ptr W;
        GR_TMP_INIT_VEC(W, 2 * lenB, ctx);
        status |= _gr_poly_div_divconquer_recursive(Q, W, A, B, lenB, invB, cutoff, ctx);
        GR_TMP_CLEAR_VEC(W, 2 * lenB, ctx);
    }

    return status;
}

/* needed due to partial overlap */
static int
_gr_vec_sub_dec(gr_ptr a, gr_srcptr b, gr_srcptr c, slong n, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    slong i;

    for (i = n - 1; i >= 0; i--)
        status |= gr_sub(GR_ENTRY(a, i, sz), GR_ENTRY(b, i, sz), GR_ENTRY(c, i, sz), ctx);

    return status;
}

int
_gr_poly_div_divconquer_preinv1(gr_ptr Q, gr_srcptr A, slong lenA,
                  gr_srcptr B, slong lenB, gr_srcptr invB, slong cutoff, gr_ctx_t ctx)
{
    if (lenA <= 2 * lenB - 1)
    {
        return __gr_poly_div_divconquer(Q, A, lenA, B, lenB, invB, cutoff, ctx);
    }
    else  /* lenA > 2 * lenB - 1 */
    {
        gr_ptr S, T, R;
        slong shift, next, n = 2 * lenB - 1;
        slong sz = ctx->sizeof_elem;
        int status = GR_SUCCESS;

        GR_TMP_INIT_VEC(S, 3 * n, ctx);
        T = GR_ENTRY(S, n, sz);
        R = GR_ENTRY(T, n, sz);

        shift = lenA - n;
        status |= _gr_vec_set(S, GR_ENTRY(A, shift, sz), n, ctx);

        while (lenA >= n)
        {
            shift = lenA - n;
            status |= _gr_poly_divrem_divconquer_recursive(GR_ENTRY(Q, shift, sz), T, R, S, B, lenB, invB, cutoff, ctx);
            next = FLINT_MIN(lenB, shift);
            status |= _gr_vec_sub_dec(GR_ENTRY(S, next, sz), S, T, lenB - 1, ctx);
            status |= _gr_vec_set(S, GR_ENTRY(A, shift - next, sz), next, ctx);
            lenA -= lenB;
        }

        if (lenA >= lenB)
            status |= __gr_poly_div_divconquer(Q, S, lenA, B, lenB, invB, cutoff, ctx);

        GR_TMP_CLEAR_VEC(S, 3 * n, ctx);

        return status;
    }
}

int
_gr_poly_div_divconquer_noinv(gr_ptr Q, gr_srcptr A, slong lenA,
                  gr_srcptr B, slong lenB, slong cutoff, gr_ctx_t ctx)
{
    return _gr_poly_div_divconquer_preinv1(Q, A, lenA, B, lenB, NULL, cutoff, ctx);
}


int
_gr_poly_div_divconquer(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong cutoff, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    gr_ptr invB;

    GR_TMP_INIT(invB, ctx);

    status = gr_inv(invB, GR_ENTRY(B, lenB - 1, sz), ctx);

    if (status == GR_SUCCESS)
        status = _gr_poly_div_divconquer_preinv1(Q, A, lenA, B, lenB, invB, cutoff, ctx);
    else
        status = _gr_poly_div_divconquer_noinv(Q, A, lenA, B, lenB, cutoff, ctx);

    GR_TMP_CLEAR(invB, ctx);

    return status;
}

int
gr_poly_div_divconquer(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, slong cutoff, gr_ctx_t ctx)
{
    slong Alen, Blen, Qlen;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    Alen = A->length;
    Blen = B->length;

    if (Blen == 0)
        return GR_DOMAIN;

    if (gr_is_zero(GR_ENTRY(B->coeffs, Blen - 1, sz), ctx) != T_FALSE)
        return GR_UNABLE;

    if (Alen < Blen)
        return gr_poly_zero(Q, ctx);

    Qlen = Alen - Blen + 1;

    if (Q == A || Q == B)
    {
        gr_poly_t t;
        gr_poly_init2(t, Qlen, ctx);
        status = _gr_poly_div_divconquer(t->coeffs, A->coeffs, A->length, B->coeffs, B->length, cutoff, ctx);
        gr_poly_swap(Q, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(Q, Qlen, ctx);
        status = _gr_poly_div_divconquer(Q->coeffs, A->coeffs, A->length, B->coeffs, B->length, cutoff, ctx);
    }

    _gr_poly_set_length(Q, Qlen, ctx);
    _gr_poly_normalise(Q, ctx);
    return status;
}
