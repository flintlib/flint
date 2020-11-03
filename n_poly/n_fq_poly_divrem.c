/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"


#define N_FQ_POLY_DIVREM_BASECASE_ITCH \
    FLINT_MAX(FLINT_MAX(4, N_FQ_MUL_ITCH), 2 + (N_FQ_REDUCE_ITCH))

void _n_fq_poly_rem_basecase_(
    mp_limb_t * Q,
    mp_limb_t * A,
    const mp_limb_t * AA, slong Alen,
    const mp_limb_t * B, slong Blen,
    const mp_limb_t * invB,
    const fq_nmod_ctx_t ctx,
    n_poly_stack_t St)
{
    slong i;
    slong d = fq_nmod_ctx_degree(ctx);
    nmod_t mod = fq_nmod_ctx_mod(ctx);
    mp_limb_t * tmp = n_poly_stack_vec_init(St, d*(3 + N_FQ_POLY_DIVREM_BASECASE_ITCH));
    mp_limb_t * u = tmp + d*N_FQ_POLY_DIVREM_BASECASE_ITCH;
    mp_limb_t * q0 = u + d;
    mp_limb_t * q1 = q0 + d;

    if (A != AA)
        _nmod_vec_set(A, AA, d*Alen);

    while (Alen - Blen > 3 && Blen > 1)
    {
        _n_fq_mul(q1, A + d*(Alen - 1), invB, ctx, tmp);
        _n_fq_mul(q0, q1, B + d*(Blen - 2), ctx, tmp);
        _n_fq_sub(q0, q0, A + d*(Alen - 2), d, mod);
        _n_fq_mul(q0, q0, invB, ctx, tmp);
        _nmod_vec_neg(q1, q1, d, ctx->mod);

        i = -1;
        _n_fq_mul(u, q0, B + d*(i + 1), ctx, tmp);
        _n_fq_add(A + d*(i + Alen - Blen), A + d*(i + Alen - Blen), u, d, mod);
        for (i = 0; i + 2 < Blen; i++)
        {
            _n_fq_mul2(tmp, q1, B + d*i, ctx);
            _n_fq_madd2(tmp, q0, B + d*(i + 1), ctx, tmp + 2*d);
            _n_fq_reduce2(u, tmp, ctx, tmp + 2*d);
            _n_fq_add(A + d*(i + Alen - Blen), A + d*(i + Alen - Blen), u, d, mod);
        }

        Alen -= 2;
        _nmod_vec_zero(A + d*Alen, 2*d);
    }

    while (Alen - Blen >= 0)
    {
        _n_fq_mul(q0, A + d*(Alen - 1), invB, ctx, tmp);

        for (i = 0; i + 1 < Blen; i++)
        {
            _n_fq_mul(u, q0, B + d*i, ctx, tmp);
            _n_fq_sub(A + d*(i + Alen - Blen), A + d*(i + Alen - Blen), u, d, mod);
        }

        Alen -= 1;
        _nmod_vec_zero(A + d*Alen, 1*d);
    }

    n_poly_stack_vec_clear(St);
}


void _n_fq_poly_divrem_basecase_(
    mp_limb_t * Q,
    mp_limb_t * A,
    const mp_limb_t * AA, slong Alen,
    const mp_limb_t * B, slong Blen,
    const mp_limb_t * invB,
    const fq_nmod_ctx_t ctx,
    n_poly_stack_t St)
{
    slong i;
    slong d = fq_nmod_ctx_degree(ctx);
    nmod_t mod = fq_nmod_ctx_mod(ctx);
    mp_limb_t * tmp = n_poly_stack_vec_init(St, d*(1 + N_FQ_POLY_DIVREM_BASECASE_ITCH));
    mp_limb_t * u = tmp + d*N_FQ_POLY_DIVREM_BASECASE_ITCH;

    if (A != AA)
        _nmod_vec_set(A, AA, d*Alen);

    while (Alen - Blen > 3 && Blen > 1)
    {
        mp_limb_t * q1 = Q + d*(Alen - Blen);
        mp_limb_t * q0 = Q + d*(Alen - Blen - 1);

        _n_fq_mul(q1, A + d*(Alen - 1), invB, ctx, tmp);
        _n_fq_mul(q0, q1, B + d*(Blen - 2), ctx, tmp);
        _n_fq_sub(q0, q0, A + d*(Alen - 2), d, mod);
        _n_fq_mul(q0, q0, invB, ctx, tmp);
        _nmod_vec_neg(q1, q1, d, ctx->mod);

        i = -1;
        _n_fq_mul(u, q0, B + d*(i + 1), ctx, tmp);
        _n_fq_add(A + d*(i + Alen - Blen), A + d*(i + Alen - Blen), u, d, mod);
        for (i = 0; i + 2 < Blen; i++)
        {
            _n_fq_mul2(tmp, q1, B + d*i, ctx);
            _n_fq_madd2(tmp, q0, B + d*(i + 1), ctx, tmp + 2*d);
            _n_fq_reduce2(u, tmp, ctx, tmp + 2*d);
            _n_fq_add(A + d*(i + Alen - Blen), A + d*(i + Alen - Blen), u, d, mod);
        }

        _nmod_vec_neg(q0, q0, 2*d, mod); /* q0 and q1 */

        Alen -= 2;
        _nmod_vec_zero(A + d*Alen, 2*d);
    }

    while (Alen - Blen >= 0)
    {
        mp_limb_t * q0 = Q + d*(Alen - Blen);

        _n_fq_mul(q0, A + d*(Alen - 1), invB, ctx, tmp);

        for (i = 0; i + 1 < Blen; i++)
        {
            _n_fq_mul(u, q0, B + d*i, ctx, tmp);
            _n_fq_sub(A + d*(i + Alen - Blen), A + d*(i + Alen - Blen), u, d, mod);
        }

        Alen -= 1;
        _nmod_vec_zero(A + d*Alen, 1*d);
    }

    n_poly_stack_vec_clear(St);
}

void _n_fq_poly_divrem_divconquer_recursive_(
    mp_limb_t * Q,
    mp_limb_t * BQ,
    mp_limb_t * W,
    const mp_limb_t * A,
    const mp_limb_t * B, slong lenB,
    const mp_limb_t * invB,
    const fq_nmod_ctx_t ctx,
    n_poly_stack_t St)
{
    slong d = fq_nmod_ctx_degree(ctx);

    if (lenB <= N_FQ_POLY_DIVREM_DIVCONQUER_CUTOFF)
    {
        _nmod_vec_zero(BQ, d*(lenB - 1));
        _nmod_vec_set(BQ + d*(lenB - 1), A + d*(lenB - 1), d*lenB);

        _n_fq_poly_divrem_basecase_(Q, BQ, BQ, 2*lenB - 1, B, lenB, invB, ctx, St);

        _nmod_vec_neg(BQ, BQ, d*(lenB - 1), ctx->mod);
        _nmod_vec_set(BQ + d*(lenB - 1), A + d*(lenB - 1), d*lenB);
    }
    else
    {
        const slong n2 = lenB / 2;
        const slong n1 = lenB - n2;
        mp_limb_t * W1 = W;
        mp_limb_t * W2 = W + d*lenB;
        const mp_limb_t * p1 = A + d*2*n2;
        const mp_limb_t * p2;
        const mp_limb_t * d1 = B + d*n2;
        const mp_limb_t * d2 = B;
        const mp_limb_t * d3 = B + d*n1;
        const mp_limb_t * d4 = B;
        mp_limb_t * q1 = Q + d*n2;
        mp_limb_t * q2 = Q;
        mp_limb_t * dq1 = BQ + d*n2;
        mp_limb_t * d1q1 = BQ + d*2*n2;
        mp_limb_t * d2q1, * d3q2, * d4q2, * t;

        _n_fq_poly_divrem_divconquer_recursive_(q1, d1q1, W1, p1, d1, n1, invB, ctx, St);

        d2q1 = W1;
        _n_fq_poly_mul_(d2q1, q1, n1, d2, n2, ctx, St);

        _nmod_vec_swap(dq1, d2q1, d*n2);
        _nmod_vec_add(dq1 + d*n2, dq1 + d*n2, d2q1 + d*n2, d*(n1 - 1), ctx->mod);

        t = BQ;
        _nmod_vec_sub(t, A + d*n2 + d*(n1 - 1), dq1 + d*(n1 - 1), d*n2, ctx->mod);
        p2 = t - d*(n2 - 1);

        d3q2 = W1;
        _n_fq_poly_divrem_divconquer_recursive_(q2, d3q2, W2, p2, d3, n2, invB, ctx, St);

        d4q2 = W2;
        _n_fq_poly_mul_(d4q2, d4, n1, q2, n2, ctx, St);

        _nmod_vec_swap(BQ, d4q2, d*n2);
        _nmod_vec_add(BQ + d*n2, BQ + d*n2, d4q2 + d*n2, d*(n1 - 1), ctx->mod);
        _nmod_vec_add(BQ + d*n1, BQ + d*n1, d3q2, d*(2*n2 - 1), ctx->mod);
    }
}

static void __n_fq_poly_divrem_divconquer_(
    mp_limb_t * Q,
    mp_limb_t * R,
    mp_limb_t * A, slong lenA,
    mp_limb_t * B, slong lenB,
    mp_limb_t * invB,
    const fq_nmod_ctx_t ctx,
    n_poly_stack_t St)
{
    slong d = fq_nmod_ctx_degree(ctx);

    FLINT_ASSERT(Q != A && Q != B);
    FLINT_ASSERT(R != A && R != B);
    FLINT_ASSERT(lenB > 0);
    FLINT_ASSERT(lenA <= 2*lenB - 1);

    if (lenA < 2*lenB - 1)
    {
        const slong n1 = lenA - lenB + 1;
        const slong n2 = lenB - n1;
        const mp_limb_t * p1 = A + d*n2;
        const mp_limb_t * d1 = B + d*n2;
        const mp_limb_t * d2 = B;
        mp_limb_t * W = n_poly_stack_vec_init(St, d*((2*n1 - 1) + lenB - 1));
        mp_limb_t * d1q1 = R + d*n2;
        mp_limb_t * d2q1 = W + d*(2*n1 - 1);

        _n_fq_poly_divrem_divconquer_recursive_(Q, d1q1, W, p1, d1, n1, invB, ctx, St);

        _n_fq_poly_mul_(d2q1, Q, n1, d2, n2, ctx, St);

        _nmod_vec_swap(R, d2q1, d*n2);
        _nmod_vec_add(R + d*n2, R + d*n2, d2q1 + d*n2, d*(n1 - 1), ctx->mod);
        _nmod_vec_sub(R, A, R, d*lenA, ctx->mod);

        n_poly_stack_vec_clear(St);
    }
    else
    {
        mp_limb_t * W = n_poly_stack_vec_init(St, d*lenA);

        _n_fq_poly_divrem_divconquer_recursive_(Q, R, W, A, B, lenB, invB, ctx, St);

        _nmod_vec_sub(R, A, R, d*(lenB - 1), ctx->mod);

        n_poly_stack_vec_clear(St);
    }
}


void _n_fq_poly_divrem_divconquer_(
    mp_limb_t * Q,
    mp_limb_t * R,
    mp_limb_t * A, slong lenA,
    mp_limb_t * B, slong lenB,
    mp_limb_t * invB,
    const fq_nmod_ctx_t ctx,
    n_poly_stack_t St)
{
    slong d = fq_nmod_ctx_degree(ctx);

    if (lenA <= 2*lenB - 1)
    {
        __n_fq_poly_divrem_divconquer_(Q, R, A, lenA, B, lenB, invB, ctx, St);
    }
    else
    {
        slong shift, n = 2*lenB - 1;
        mp_limb_t * QB, * W;

        _nmod_vec_set(R, A, d*lenA);
        W = n_poly_stack_vec_init(St, d*2*n);
        QB = W + d*n;

        while (lenA >= n)
        {
            shift = lenA - n;
            _n_fq_poly_divrem_divconquer_recursive_(Q + d*shift, QB,
                                       W, R + d*shift, B, lenB, invB, ctx, St);
            _nmod_vec_sub(R + d*shift, R + d*shift, QB, d*n, ctx->mod);
            lenA -= lenB;
        }

        if (lenA >= lenB)
        {
            __n_fq_poly_divrem_divconquer_(Q, W, R, lenA, B, lenB, invB, ctx, St);
            _nmod_vec_swap(W, R, d*lenA);
        }

        n_poly_stack_vec_clear(St);
    }
}

void n_fq_poly_divrem_divconquer_(
    n_fq_poly_t Q,
    n_fq_poly_t R,
    const n_fq_poly_t A,
    const n_fq_poly_t B,
    const fq_nmod_ctx_t ctx,
    n_poly_stack_t St)
{
    slong d = fq_nmod_ctx_degree(ctx);
    const slong lenA = A->length;
    const slong lenB = B->length;
    const slong lenQ = lenA - lenB + 1;
    mp_limb_t * tmp, * invB;
    n_poly_t Qt, Rt;
    mp_limb_t * q, * r;
#if FLINT_WANT_ASSERT
    fq_nmod_poly_t QQ, RR, AA, BB;
#endif

    if (lenQ < 1)
    {
        n_fq_poly_set(R, A, ctx);
        n_poly_zero(Q);
        return;
    }

#if FLINT_WANT_ASSERT
    fq_nmod_poly_init(QQ, ctx);
    fq_nmod_poly_init(RR, ctx);
    fq_nmod_poly_init(AA, ctx);
    fq_nmod_poly_init(BB, ctx);
    n_fq_poly_get_fq_nmod_poly(AA, A, ctx);
    n_fq_poly_get_fq_nmod_poly(BB, B, ctx);
    fq_nmod_poly_divrem(QQ, RR, AA, BB, ctx);
#endif

    tmp = n_poly_stack_vec_init(St, d*N_FQ_INV_ITCH + d);
    invB = tmp + d*N_FQ_INV_ITCH;

    _n_fq_inv(invB, B->coeffs + d*(lenB - 1), ctx, tmp);

    if (Q == A || Q == B)
    {
        n_fq_poly_init(Qt);
        n_poly_fit_length(Qt, d*lenQ);
        q = Qt->coeffs;
    }
    else
    {
        n_poly_fit_length(Q, d*lenQ);
        q = Q->coeffs;
    }

    /* TODO why lenA here and not lenB ? */
    if (R == A || R == B)
    {
        n_fq_poly_init(Rt);
        n_poly_fit_length(Rt, d*lenA);
        r = Rt->coeffs;
    }
    else
    {
        n_poly_fit_length(R, d*lenA);
        r = R->coeffs;
    }

    _n_fq_poly_divrem_divconquer_(q, r, A->coeffs, lenA, B->coeffs, lenB, invB, ctx, St);

    if (Q == A || Q == B)
    {
        n_fq_poly_swap(Q, Qt);
        n_fq_poly_clear(Qt);
    }

    Q->length = lenQ;

    if (R == A || R == B)
    {
        n_fq_poly_swap(R, Rt);
        n_fq_poly_clear(Rt);
    }

    R->length = lenB - 1;
    _n_fq_poly_normalise(R, d);

    n_poly_stack_vec_clear(St);

#if FLINT_WANT_ASSERT
    n_fq_poly_get_fq_nmod_poly(AA, Q, ctx);
    n_fq_poly_get_fq_nmod_poly(BB, R, ctx);
    FLINT_ASSERT(fq_nmod_poly_equal(AA, QQ, ctx));
    FLINT_ASSERT(fq_nmod_poly_equal(BB, RR, ctx));
    fq_nmod_poly_clear(QQ, ctx);
    fq_nmod_poly_clear(RR, ctx);
    fq_nmod_poly_clear(AA, ctx);
    fq_nmod_poly_clear(BB, ctx);
#endif
}

void n_fq_poly_divrem(
    n_fq_poly_t Q,
    n_fq_poly_t R,
    const n_fq_poly_t A,
    const n_fq_poly_t B,
    const fq_nmod_ctx_t ctx)
{
    n_poly_stack_t St;
    n_poly_stack_init(St);
    n_fq_poly_divrem_divconquer_(Q, R, A, B, ctx, St);
    n_poly_stack_clear(St);
}
