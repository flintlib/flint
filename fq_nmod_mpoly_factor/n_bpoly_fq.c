/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"


void n_bpoly_fq_get_coeff_fq_nmod(
    fq_nmod_t c,
    const n_bpoly_t A,
    slong e0,
    slong e1,
    const fq_nmod_ctx_t ctx)
{
    if (e0 >= A->length)
        fq_nmod_zero(c, ctx);
    else
        n_poly_fq_get_coeff_fq_nmod(c, A->coeffs + e0, e1, ctx);
}


void n_bpoly_fq_print_pretty(
    const n_bpoly_t A,
    const char * var0,
    const char * var1,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    int first;

    first = 1;
    for (i = A->length - 1; i >= 0; i--)
    {
        if (n_poly_is_zero(A->coeffs + i))
            continue;

        if (!first)
            flint_printf(" + ");

        first = 0;

        flint_printf("(");
        n_poly_fq_print_pretty(A->coeffs + i, var1, ctx);
        flint_printf(")*%s^%wd", var0, i);
    }

    if (first)
        flint_printf("0");
}

int n_bpoly_fq_is_canonical(const n_bpoly_t A, const fq_nmod_ctx_t ctx)
{
    slong i;

    if (A->length < 0)
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (!n_poly_fq_is_canonical(A->coeffs + i, ctx))
            return 0;

        if (i + 1 == A->length && n_poly_is_zero(A->coeffs + i))
            return 0;
    }

    return 1;
}

int n_bpoly_fq_equal(
    const n_bpoly_t A,
    const n_bpoly_t B,
    const fq_nmod_ctx_t ctx)
{
    slong i;

    if (A->length != B->length)
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (!n_poly_fq_equal(A->coeffs + i, B->coeffs + i, ctx))
            return 0;
    }

    return 1;
}

void n_bpoly_fq_set_coeff_n_fq(
    n_bpoly_t A,
    slong xi,
    slong yi,
    const mp_limb_t *c,
    const fq_nmod_ctx_t ctx)
{
    slong i;

    if (xi >= A->length)
    {
        n_bpoly_fit_length(A, xi + 1);
        for (i = A->length; i <= xi; i++)
            n_poly_zero(A->coeffs + i);
        A->length = xi + 1;
    }

    n_poly_fq_set_coeff_n_fq(A->coeffs + xi, yi, c, ctx);
    n_bpoly_normalise(A);
}

void n_bpoly_fq_set_coeff_fq_nmod(
    n_bpoly_t A,
    slong xi,
    slong yi,
    const fq_nmod_t c,
    const fq_nmod_ctx_t ctx)
{
    slong i;

    if (xi >= A->length)
    {
        n_bpoly_fit_length(A, xi + 1);
        for (i = A->length; i <= xi; i++)
            n_poly_zero(A->coeffs + i);
        A->length = xi + 1;
    }

    n_poly_fq_set_coeff_fq_nmod(A->coeffs + xi, yi, c, ctx);
    n_bpoly_normalise(A);
}


void n_bpoly_fq_set_fq_nmod_poly_var0(
    n_bpoly_t A,
    const fq_nmod_poly_t B,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    n_bpoly_fit_length(A, B->length);
    A->length = 0;
    for (i = 0; i < B->length; i++)
        n_poly_fq_set_fq_nmod(A->coeffs + i, B->coeffs + i, ctx);
    A->length = B->length;
    n_bpoly_normalise(A);
}

void n_bpoly_fq_set_n_poly_fq_var0(
    n_bpoly_t A,
    const n_poly_t B,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    n_bpoly_fit_length(A, B->length);
    for (i = 0; i < B->length; i++)
        n_poly_fq_set_n_fq(A->coeffs + i, B->coeffs + d*i, ctx);
    A->length = B->length;
    n_bpoly_normalise(A);
}

void n_bpoly_fq_set_fq_nmod_poly_var1(
    n_bpoly_t A,
    const fq_nmod_poly_t B,
    const fq_nmod_ctx_t ctx)
{
    n_bpoly_fit_length(A, 1);
	n_poly_fq_set_fq_nmod_poly(A->coeffs + 0, B, ctx);
	A->length = !n_poly_is_zero(A->coeffs + 0);
}

void n_bpoly_fq_set_n_poly_fq_var1(
    n_bpoly_t A,
    const n_poly_t B,
    const fq_nmod_ctx_t ctx)
{
    n_bpoly_fit_length(A, 1);
	n_poly_fq_set(A->coeffs + 0, B, ctx);
	A->length = !n_poly_is_zero(A->coeffs + 0);
}


void n_bpoly_fq_mul(
    n_bpoly_t A,
    const n_bpoly_t B,
    const n_bpoly_t C,
    const fq_nmod_ctx_t ctx)
{
    slong i, j;
    n_poly_struct * t;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    if (B->length < 1 || C->length < 1)
    {
        A->length = 0;
        return;
    }

    n_bpoly_fit_length(A, B->length + C->length);
    for (i = 0; i < B->length + C->length - 1; i++)
        n_poly_zero(A->coeffs + i);

    t = A->coeffs + B->length + C->length - 1;

    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < C->length; j++)
        {
            n_poly_fq_mul(t, B->coeffs + i, C->coeffs + j, ctx);
            n_poly_fq_add(A->coeffs + i + j, A->coeffs + i + j, t, ctx);
        }
    }

    A->length = B->length + C->length - 1;
    n_bpoly_normalise(A);
}


void n_bpoly_fq_mul_series(
    n_bpoly_t A,
    const n_bpoly_t B,
    const n_bpoly_t C,
    slong order,
    const fq_nmod_ctx_t ctx)
{
    slong i, j;
    n_poly_t t;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    n_poly_init(t);

    n_bpoly_fit_length(A, B->length + C->length - 1);
    for (i = 0; i < B->length + C->length - 1; i++)
        n_poly_zero(A->coeffs + i);

    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < C->length; j++)
        {
            n_poly_fq_mullow(t, B->coeffs + i, C->coeffs + j, order, ctx);
            n_poly_fq_add(A->coeffs + i + j, A->coeffs + i + j, t, ctx);
        }
    }

    n_poly_clear(t);

    A->length = B->length + C->length - 1;
    n_bpoly_normalise(A);
}


void n_bpoly_fq_add(
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
                n_poly_fq_add(A->coeffs + i, B->coeffs + i, C->coeffs + i, ctx);
            }
            else
            {
                n_poly_fq_set(A->coeffs + i, B->coeffs + i, ctx);
            }
        }
        else
        {
            FLINT_ASSERT(i < C->length);
            n_poly_fq_set(A->coeffs + i, C->coeffs + i, ctx);
        }
    }

    A->length = Alen;
    n_bpoly_normalise(A);
}

void n_bpoly_fq_one(n_bpoly_t A, const fq_nmod_ctx_t ctx)
{
    n_bpoly_fit_length(A, 1);
    A->length = 1;
    n_poly_fq_one(A->coeffs + 0, ctx);
}

void n_bpoly_fq_sub(
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
                n_poly_fq_sub(A->coeffs + i, B->coeffs + i, C->coeffs + i, ctx);
            }
            else
            {
                n_poly_fq_set(A->coeffs + i, B->coeffs + i, ctx);
            }
        }
        else
        {
            FLINT_ASSERT(i < C->length);
            n_poly_fq_neg(A->coeffs + i, C->coeffs + i, ctx);
        }
    }

    A->length = Alen;
    n_bpoly_normalise(A);
}

void n_bpoly_fq_derivative(
    n_bpoly_t A,
    const n_bpoly_t B,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    slong Blen = B->length;

    if (Blen < 2)
    {
        n_bpoly_zero(A);
        return;
    }

    n_bpoly_fit_length(A, Blen - 1);

    for (i = 1; i < Blen; i++)
        n_poly_fq_scalar_mul_ui(A->coeffs + i - 1, B->coeffs + i, i, ctx);

    A->length = Blen - 1;
    n_bpoly_normalise(A);
}

/*
    division in (Fq[y]/y^order)[x]
    inputs need not be reduced mod y^order
*/
void n_bpoly_fq_divrem_series(
    n_bpoly_t Q,
    n_bpoly_t R,
    const n_bpoly_t A,
    const n_bpoly_t B,
    slong order,
    const fq_nmod_ctx_t ctx)
{
    slong i, qoff;
    n_poly_t q, t, binv;

    FLINT_ASSERT(R != A);
    FLINT_ASSERT(R != B);
    FLINT_ASSERT(Q != A);
    FLINT_ASSERT(Q != B);
    FLINT_ASSERT(B->length > 0);

    n_poly_init(q);
    n_poly_init(t);
    n_poly_init(binv);

    n_bpoly_fq_set(R, A, ctx);
    for (i = 0; i < R->length; i++)
        n_poly_fq_truncate(R->coeffs + i, order, ctx);
    n_bpoly_normalise(R);

    n_poly_fq_inv_series(binv, B->coeffs + B->length - 1, order, ctx);

    Q->length = 0;

    while (R->length >= B->length)
    {
        n_poly_fq_mullow(q, R->coeffs + R->length - 1, binv, order, ctx);

        for (i = 0; i < B->length; i++)
        {
            n_poly_fq_mullow(t, B->coeffs + i, q, order, ctx);
            n_poly_fq_sub(R->coeffs + i + R->length - B->length,
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

        n_poly_fq_set(Q->coeffs + qoff, q, ctx);

        FLINT_ASSERT(n_poly_is_zero(R->coeffs + R->length - 1));
        n_bpoly_normalise(R);
    }

    n_poly_clear(q);
    n_poly_clear(t);
    n_poly_clear(binv);
}


int n_bpoly_fq_divides(
    n_bpoly_t Q,
    const n_bpoly_t A,
    const n_bpoly_t B,
    const fq_nmod_ctx_t ctx)
{
    slong i, qoff;
    int divides;
    n_poly_t q, t;
    n_bpoly_t R;

    FLINT_ASSERT(Q != A);
    FLINT_ASSERT(Q != B);
    FLINT_ASSERT(B->length > 0);

    n_poly_init(q);
    n_poly_init(t);
    n_bpoly_init(R);
    n_bpoly_fq_set(R, A, ctx);

    Q->length = 0;

    while (R->length >= B->length)
    {
        n_poly_fq_divrem(q, t, R->coeffs + R->length - 1,
                                               B->coeffs + B->length - 1, ctx);
        if (!n_poly_is_zero(t))
        {
            divides = 0;
            goto cleanup;
        }

        for (i = 0; i < B->length; i++)
        {
            n_poly_fq_mul(t, B->coeffs + i, q, ctx);
            n_poly_fq_sub(R->coeffs + i + R->length - B->length,
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

        n_poly_fq_set(Q->coeffs + qoff, q, ctx);

        FLINT_ASSERT(n_poly_is_zero(R->coeffs + R->length - 1));
        n_bpoly_normalise(R);
    }

    divides = (R->length == 0);

cleanup:

    n_poly_clear(q);
    n_poly_clear(t);
    n_bpoly_clear(R);

    return divides;
}

void n_bpoly_fq_set(
    n_bpoly_t A,
    const n_bpoly_t B,
    const fq_nmod_ctx_t ctx)
{
    slong i;

    if (A == B)
        return;

    n_bpoly_fit_length(A, B->length);
    A->length = B->length;

    for (i = 0; i < B->length; i++)
        n_poly_fq_set(A->coeffs + i, B->coeffs + i, ctx);
}

void n_bpoly_fq_make_primitive(
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
        n_poly_fq_gcd(q, g, A->coeffs + i, ctx);
		n_poly_swap(g, q);
	}

    for (i = 0; i < Alen; i++)
    {
        n_poly_fq_divrem(q, r, A->coeffs + i, g, ctx);
		n_poly_fq_set(A->coeffs + i, q, ctx);
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
            n_poly_fq_scalar_mul_n_fq(g, g, c_, ctx);
            n_fq_inv(c, c_, ctx);
            for (i = 0; i < Alen; i++)
                n_poly_fq_scalar_mul_n_fq(A->coeffs + i, A->coeffs + i, c, ctx);
        }
        flint_free(c);
    }

    n_poly_clear(q);
    n_poly_clear(r);
}


void _n_poly_fq_taylor_shift_horner_n_fq(
    mp_limb_t * poly,
    const mp_limb_t * c,
    slong n,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i, j;
    mp_limb_t * p = FLINT_ARRAY_ALLOC(d, mp_limb_t);

    for (i = n - 2; i >= 0; i--)
    {
        for (j = i; j < n - 1; j++)
        {
            n_fq_mul(p, poly + d*(j + 1), c, ctx);
            n_fq_add(poly + d*j, poly + d*j, p, ctx);
        }
    }

    flint_free(p);
}

void n_poly_fq_taylor_shift_horner_n_fq(
    n_poly_t g,
    const n_poly_t f,
    const mp_limb_t * c,
    const fq_nmod_ctx_t ctx)
{
    if (f != g)
        n_poly_fq_set(g, f, ctx);

    _n_poly_fq_taylor_shift_horner_n_fq(g->coeffs, c, g->length, ctx);
}


void n_bpoly_fq_taylor_shift_var1_fq_nmod(
    n_bpoly_t A,
    const n_bpoly_t B,
    const fq_nmod_t c_,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    mp_limb_t * c = FLINT_ARRAY_ALLOC(d, mp_limb_t);

    n_fq_set_fq_nmod(c, c_, ctx);
    n_bpoly_fq_set(A, B, ctx);
    for (i = A->length - 1; i >= 0; i--)
        _n_poly_fq_taylor_shift_horner_n_fq(A->coeffs[i].coeffs, c, A->coeffs[i].length, ctx);

    flint_free(c);  
}

void n_bpoly_fq_taylor_shift_var0_fq_nmod(
    n_bpoly_t A,
    const fq_nmod_t alpha,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong n, i, j;
    mp_limb_t * c;
    n_poly_t t;

    if (fq_nmod_is_zero(alpha, ctx))
        return;

    c = FLINT_ARRAY_ALLOC(d, mp_limb_t);
    n_fq_set_fq_nmod(c, alpha, ctx);

    n_poly_init(t);
    n = A->length;

    for (i = n - 2; i >= 0; i--)
    {
        for (j = i; j < n - 1; j++)
        {
            n_poly_fq_scalar_mul_n_fq(t, A->coeffs + j + 1, c, ctx);
            n_poly_fq_add(A->coeffs + j, A->coeffs + j, t, ctx);
        }
    }

    n_poly_clear(t);

    flint_free(c);
}


void fq_nmod_mpoly_get_n_bpoly_fq(
    n_bpoly_t A,
    const fq_nmod_mpoly_t B,
    slong varx,
    slong vary,
    const fq_nmod_mpoly_ctx_t ctx)
{
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
        n_bpoly_fq_set_coeff_fq_nmod(A, Bexpx, Bexpy, B->coeffs + j, ctx->fqctx);
    }
}


void fq_nmod_mpoly_set_n_bpoly_fq(
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
    slong Alen;
    fq_nmod_struct * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    ulong * Aexps;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    TMP_START;

    Aexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));
    for (i = 0; i < n; i++)
        Aexps[i] = 0;

    NA = mpoly_words_per_exp(Abits, ctx->minfo);
    fq_nmod_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        n_poly_struct * Bc = B->coeffs + i;
        _fq_nmod_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Bc->length, NA, ctx->fqctx);

        for (j = 0; j < Bc->length; j++)
        {
            if (_n_fq_is_zero(Bc->coeffs + d*j, d))
                continue;

            Aexps[varx] = i;
            Aexps[vary] = j;
            n_fq_get_fq_nmod(Acoeff + Alen, Bc->coeffs + d*j, ctx->fqctx);
            mpoly_set_monomial_ui(Aexp + NA*Alen, Aexps, Abits, ctx->minfo);
            Alen++;
        }
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    A->length = Alen;

    TMP_END;

    fq_nmod_mpoly_sort_terms(A, ctx);
}

