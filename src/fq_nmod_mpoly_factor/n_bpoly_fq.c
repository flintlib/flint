/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"

void n_fq_poly_mullow_(
    n_fq_poly_t A,
    const n_fq_poly_t B,
    const n_fq_poly_t C,
    slong order,
    const fq_nmod_ctx_t ctx,
    n_poly_stack_t St);

void n_fq_bpoly_add(
    n_bpoly_t A,
    const n_bpoly_t B,
    const n_bpoly_t C,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    slong Alen = FLINT_MAX(B->length, C->length);

    n_bpoly_fit_length(A, Alen);

    for (i = 0; i < Alen; i++)
    {
        if (i < B->length)
        {
            if (i < C->length)
            {
                n_fq_poly_add(A->coeffs + i, B->coeffs + i, C->coeffs + i, ctx);
            }
            else
            {
                n_fq_poly_set(A->coeffs + i, B->coeffs + i, ctx);
            }
        }
        else
        {
            FLINT_ASSERT(i < C->length);
            n_fq_poly_set(A->coeffs + i, C->coeffs + i, ctx);
        }
    }

    A->length = Alen;
    n_bpoly_normalise(A);
}


void n_fq_bpoly_sub(
    n_bpoly_t A,
    const n_bpoly_t B,
    const n_bpoly_t C,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    slong Alen = FLINT_MAX(B->length, C->length);

    n_bpoly_fit_length(A, Alen);

    for (i = 0; i < Alen; i++)
    {
        if (i < B->length)
        {
            if (i < C->length)
            {
                n_fq_poly_sub(A->coeffs + i, B->coeffs + i, C->coeffs + i, ctx);
            }
            else
            {
                n_fq_poly_set(A->coeffs + i, B->coeffs + i, ctx);
            }
        }
        else
        {
            FLINT_ASSERT(i < C->length);
            n_fq_poly_neg(A->coeffs + i, C->coeffs + i, ctx);
        }
    }

    A->length = Alen;
    n_bpoly_normalise(A);
}


void n_fq_bpoly_mul(
    n_bpoly_t A,
    const n_bpoly_t B,
    const n_bpoly_t C,
    const fq_nmod_ctx_t ctx)
{
    slong i, j;
    n_poly_struct * t;
    n_poly_stack_t St;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    if (B->length < 1 || C->length < 1)
    {
        A->length = 0;
        return;
    }

    n_poly_stack_init(St);
    n_poly_stack_fit_request(St, 1);
    t = n_poly_stack_take_top(St);

    n_bpoly_fit_length(A, B->length + C->length - 1);
    for (i = 0; i < B->length + C->length - 1; i++)
        n_poly_zero(A->coeffs + i);

    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < C->length; j++)
        {
            n_fq_poly_mul_(t, B->coeffs + i, C->coeffs + j, ctx, St);
            n_fq_poly_add(A->coeffs + i + j, A->coeffs + i + j, t, ctx);
        }
    }

    A->length = B->length + C->length - 1;
    n_bpoly_normalise(A);

    n_poly_stack_give_back(St, 1);
    n_poly_stack_clear(St);
}


void n_fq_bpoly_mul_series(
    n_fq_bpoly_t A,
    const n_fq_bpoly_t B,
    const n_fq_bpoly_t C,
    slong order,
    const fq_nmod_ctx_t ctx)
{
    slong i, j;
    n_poly_struct * t;
    n_poly_stack_t St;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    n_poly_stack_init(St);
    n_poly_stack_fit_request(St, 1);
    t = n_poly_stack_take_top(St);

    n_fq_bpoly_fit_length(A, B->length + C->length - 1);
    for (i = 0; i < B->length + C->length - 1; i++)
        n_fq_poly_zero(A->coeffs + i);

    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < C->length; j++)
        {
            n_fq_poly_mullow_(t, B->coeffs + i, C->coeffs + j, order, ctx, St);
            n_fq_poly_add(A->coeffs + i + j, A->coeffs + i + j, t, ctx);
        }
    }

    A->length = B->length + C->length - 1;
    n_bpoly_normalise(A);

    n_poly_stack_give_back(St, 1);
    n_poly_stack_clear(St);
}

/*
    division in (Fq[y]/y^order)[x]
    inputs need not be reduced mod y^order
*/
void n_fq_bpoly_divrem_series(
    n_bpoly_t Q,
    n_bpoly_t R,
    const n_bpoly_t A,
    const n_bpoly_t B,
    slong order,
    const fq_nmod_ctx_t ctx)
{
    slong i, qoff;
    n_poly_struct * q, * t, * binv;
    n_poly_stack_t St;

    FLINT_ASSERT(R != A);
    FLINT_ASSERT(R != B);
    FLINT_ASSERT(Q != A);
    FLINT_ASSERT(Q != B);
    FLINT_ASSERT(B->length > 0);

    n_poly_stack_init(St);
    n_poly_stack_fit_request(St, 3);
    q = n_poly_stack_take_top(St);
    t = n_poly_stack_take_top(St);
    binv = n_poly_stack_take_top(St);

    n_fq_bpoly_set(R, A, ctx);
    for (i = 0; i < R->length; i++)
        n_fq_poly_truncate(R->coeffs + i, order, ctx);
    n_bpoly_normalise(R);

    n_fq_poly_inv_series(binv, B->coeffs + B->length - 1, order, ctx);

    Q->length = 0;

    while (R->length >= B->length)
    {
        n_fq_poly_mullow_(q, R->coeffs + R->length - 1, binv, order, ctx, St);

        for (i = 0; i < B->length; i++)
        {
            n_fq_poly_mullow_(t, B->coeffs + i, q, order, ctx, St);
            n_fq_poly_sub(R->coeffs + i + R->length - B->length,
                                R->coeffs + i + R->length - B->length, t, ctx);
        }

        qoff = R->length - B->length;

        FLINT_ASSERT(qoff >= 0);
        if (qoff >= Q->length)
        {
            n_bpoly_fit_length(Q, qoff + 1);
            for (i = Q->length; i <= qoff; i++)
                n_poly_zero(Q->coeffs + i);
            Q->length = qoff + 1;
        }

        n_fq_poly_set(Q->coeffs + qoff, q, ctx);

        FLINT_ASSERT(n_poly_is_zero(R->coeffs + R->length - 1));
        n_bpoly_normalise(R);
    }

    n_poly_stack_give_back(St, 3);
    n_poly_stack_clear(St);
}


int n_fq_bpoly_divides(
    n_bpoly_t Q,
    const n_bpoly_t A,
    const n_bpoly_t B,
    const fq_nmod_ctx_t ctx)
{
    slong i, qoff;
    int divides;
    n_fq_poly_struct * q, * t;
    n_fq_bpoly_t R;
    n_poly_stack_t St;

    FLINT_ASSERT(Q != A);
    FLINT_ASSERT(Q != B);
    FLINT_ASSERT(B->length > 0);

    n_poly_stack_init(St);
    n_poly_stack_fit_request(St, 2);
    q = n_poly_stack_take_top(St);
    t = n_poly_stack_take_top(St);

    n_bpoly_init(R);
    n_fq_bpoly_set(R, A, ctx);

    Q->length = 0;

    while (R->length >= B->length)
    {
        n_fq_poly_divrem_(q, t, R->coeffs + R->length - 1,
                                           B->coeffs + B->length - 1, ctx, St);
        if (!n_poly_is_zero(t))
        {
            divides = 0;
            goto cleanup;
        }

        for (i = 0; i < B->length; i++)
        {
            n_fq_poly_mul_(t, B->coeffs + i, q, ctx, St);
            n_fq_poly_sub(R->coeffs + i + R->length - B->length,
                                R->coeffs + i + R->length - B->length, t, ctx);
        }

        qoff = R->length - B->length;

        FLINT_ASSERT(qoff >= 0);
        if (qoff >= Q->length)
        {
            n_bpoly_fit_length(Q, qoff + 1);
            for (i = Q->length; i <= qoff; i++)
                n_poly_zero(Q->coeffs + i);
            Q->length = qoff + 1;
        }

        n_fq_poly_set(Q->coeffs + qoff, q, ctx);

        FLINT_ASSERT(n_fq_poly_is_zero(R->coeffs + R->length - 1));
        n_bpoly_normalise(R);
    }

    divides = (R->length == 0);

cleanup:

    n_poly_stack_give_back(St, 2);
    n_poly_stack_clear(St);
    n_bpoly_clear(R);

    return divides;
}


void n_fq_bpoly_make_primitive(
    n_poly_t g,
    n_bpoly_t A,
    const fq_nmod_ctx_t ctx)
{
    slong Alen = A->length;
    slong i;
    n_poly_t q, r;

    n_poly_init(q);
    n_poly_init(r);

    n_poly_zero(g);
    for (i = 0; i < Alen; i++)
	{
        n_fq_poly_gcd(q, g, A->coeffs + i, ctx);
		n_poly_swap(g, q);
	}

    for (i = 0; i < Alen; i++)
    {
        n_fq_poly_divrem(q, r, A->coeffs + i, g, ctx);
		n_fq_poly_set(A->coeffs + i, q, ctx);
    }

    /* make lc_xy(A) one */
    if (Alen > 0)
    {
        slong d = fq_nmod_ctx_degree(ctx);
        n_poly_struct * Alead = A->coeffs + Alen - 1;
        mp_limb_t * c, * c_ = Alead->coeffs + d*(Alead->length - 1);
        c = FLINT_ARRAY_ALLOC(d, mp_limb_t);
        if (!_n_fq_is_one(c_, d))
        {
            n_fq_poly_scalar_mul_n_fq(g, g, c_, ctx);
            n_fq_inv(c, c_, ctx);
            for (i = 0; i < Alen; i++)
                n_fq_poly_scalar_mul_n_fq(A->coeffs + i, A->coeffs + i, c, ctx);
        }
        flint_free(c);
    }

    n_poly_clear(q);
    n_poly_clear(r);
}


void fq_nmod_mpoly_get_n_fq_bpoly(
    n_bpoly_t A,
    const fq_nmod_mpoly_t B,
    slong varx,
    slong vary,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong j;
    slong NB;
    ulong Bexpx, Bexpy;
    slong Boffx, Bshiftx, Boffy, Bshifty;
    ulong mask;

    FLINT_ASSERT(B->bits <= FLINT_BITS);
    NB = mpoly_words_per_exp_sp(B->bits, ctx->minfo);

    mpoly_gen_offset_shift_sp(&Boffx, &Bshiftx, varx, B->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&Boffy, &Bshifty, vary, B->bits, ctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - B->bits);

    n_bpoly_zero(A);
    for (j = 0; j < B->length; j++)
    {
        Bexpx = ((B->exps + NB*j)[Boffx] >> Bshiftx) & mask;
        Bexpy = ((B->exps + NB*j)[Boffy] >> Bshifty) & mask;
        n_fq_bpoly_set_coeff_n_fq(A, Bexpx, Bexpy, B->coeffs + d*j, ctx->fqctx);
    }
}


void fq_nmod_mpoly_set_n_fq_bpoly(
    fq_nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const n_bpoly_t B,
    slong varx,
    slong vary,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong n = ctx->minfo->nvars;
    slong i, j;
    slong NA;
    ulong * Aexps;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    TMP_START;

    Aexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));
    for (i = 0; i < n; i++)
        Aexps[i] = 0;

    fq_nmod_mpoly_fit_length_reset_bits(A, 4, Abits, ctx);
    NA = mpoly_words_per_exp(Abits, ctx->minfo);

    A->length = 0;
    for (i = 0; i < B->length; i++)
    {
        n_poly_struct * Bc = B->coeffs + i;
        for (j = 0; j < Bc->length; j++)
        {
            Aexps[varx] = i;
            Aexps[vary] = j;

            if (_n_fq_is_zero(Bc->coeffs + d*j, d))
                continue;

            fq_nmod_mpoly_fit_length(A, A->length + 1, ctx);

            _n_fq_set(A->coeffs + d*A->length, Bc->coeffs + d*j, d);
            mpoly_set_monomial_ui(A->exps + NA*A->length, Aexps, Abits, ctx->minfo);
            A->length++;
        }
    }

    TMP_END;

    fq_nmod_mpoly_sort_terms(A, ctx);
}

