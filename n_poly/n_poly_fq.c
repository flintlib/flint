/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"
#include "nmod_vec.h"


void n_poly_fq_print_pretty(
    const n_poly_t A,
    const char * x,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    int first;

    first = 1;
    for (i = A->length - 1; i >= 0; i--)
    {
        if (_n_fq_is_zero(A->coeffs + d*i, d))
            continue;

        if (!first)
            flint_printf(" + ");

        first = 0;

        flint_printf("(");
        n_fq_print_pretty(A->coeffs + d*i, ctx);
        flint_printf(")*%s^%wd", x, i);
    }

    if (first)
        flint_printf("0");
}

int n_poly_fq_is_canonical(
    const n_poly_t A,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;

    if (A->length < 0)
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (!n_fq_is_canonical(A->coeffs + d*i, ctx))
            return 0;
        if (_n_fq_is_zero(A->coeffs + d*i, d) && i + 1 == A->length)
            return 0;
    }
    return 1;
}

void n_poly_fq_randtest(
    n_poly_t A,
    flint_rand_t state,
    slong len,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;

    if (len < 1)
    {
        A->length = 0;
        return;
    }

    n_poly_fit_length(A, d*len);

    for (i = 0; i < d*len; i++)
        A->coeffs[i] = n_randint(state, ctx->mod.n);

    A->length = len;
    _n_poly_fq_normalise(A, d);
}


void _n_poly_fq_one(n_poly_t A, slong d)
{
    n_poly_fit_length(A, d);
    A->length = 1;
    _n_fq_one(A->coeffs + d*0, d);
}


int n_poly_fq_is_one(n_poly_t A, const fq_nmod_ctx_t ctx)
{
    return A->length == 1 && _n_fq_is_one(A->coeffs, fq_nmod_ctx_degree(ctx));
}


void n_poly_fq_get_coeff_fq_nmod(
    fq_nmod_t c,
    const n_poly_t A,
    slong e,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);

    if (e >= A->length)
        fq_nmod_zero(c, ctx);
    else
        n_fq_get_fq_nmod(c, A->coeffs + d*e, ctx);
}

void n_poly_fq_set_coeff_n_fq(
    n_poly_t A,
    slong j,
    const mp_limb_t * c,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);

    n_poly_fit_length(A, d*(j + 1));

    if (j + 1 <= A->length)
    {
        _n_fq_set(A->coeffs + d*j, c, d);
        if (j + 1 == A->length)
            _n_poly_fq_normalise(A, d);
    }
    else if (!_n_fq_is_zero(c, d)) /* extend polynomial */
    {
        flint_mpn_zero(A->coeffs + A->length, d*(j - A->length));
        _n_fq_set(A->coeffs + d*j, c, d);
        A->length = j + 1;
    }
}

void n_poly_fq_set_coeff_fq_nmod(
    n_poly_t A,
    slong j,
    const fq_nmod_t c,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);

    n_poly_fit_length(A, d*(j + 1));

    if (j + 1 <= A->length)
    {
        n_fq_set_fq_nmod(A->coeffs + d*j, c, ctx);
        if (j + 1 == A->length)
            _n_poly_fq_normalise(A, d);
    }
    else if (!fq_nmod_is_zero(c, ctx))
    {
        flint_mpn_zero(A->coeffs + A->length, d*(j - A->length));
        n_fq_set_fq_nmod(A->coeffs + d*j, c, ctx);
        A->length = j + 1;
    }
}


void n_poly_fq_scalar_mul_ui(
    n_poly_t A,
    const n_poly_t B,
    ulong c,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);

    if (c >= ctx->mod.n)
    {
        NMOD_RED(c, c, ctx->mod);
    }

    if (B->length < 1 || c == 0)
    {
        n_poly_zero(A);
        return;
    }

    n_poly_fit_length(A, d*B->length);
    _nmod_vec_scalar_mul_nmod(A->coeffs, B->coeffs, d*B->length, c, ctx->mod);
    A->length = B->length;
    _n_poly_fq_normalise(A, d);
}


int n_poly_fq_equal(
    const n_poly_t A,
    const n_poly_t B,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;

    if (A->length != B->length)
        return 0;

    for (i = 0; i < d*B->length; i++)
        if (A->coeffs[i] != B->coeffs[i])
            return 0;

    return 1;
}

void n_poly_fq_set(
    n_poly_t A,
    const n_poly_t B,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);

    if (A == B)
        return;

    n_poly_fit_length(A, d*B->length);
    _nmod_vec_set(A->coeffs, B->coeffs, d*B->length);
    A->length = B->length;
}

void n_poly_fq_set_n_fq(
    n_poly_t A,
    const mp_limb_t * c,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);

    n_poly_fit_length(A, d);
    _nmod_vec_set(A->coeffs, c, d);
    A->length = 1;
    _n_poly_fq_normalise(A, d);
}

void n_poly_fq_set_fq_nmod(
    n_poly_t A,
    const fq_nmod_t c,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);

    n_poly_fit_length(A, d);
    n_fq_set_fq_nmod(A->coeffs, c, ctx);
    A->length = 1;
    _n_poly_fq_normalise(A, d);
}

void n_poly_fq_shift_right(
    n_poly_t A,
    const n_poly_t B,
    slong n,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);

    if (n < 1)
    {
        n_poly_fq_set(A, B, ctx);
        return;
    }
    else if (B->length <= n)
    {
        n_poly_zero(A);
        return;
    }
    else
    {
        n_poly_fit_length(A, d*(B->length - n));
        flint_mpn_copyi(A->coeffs, B->coeffs + d*n, d*(B->length - n));
        A->length = B->length - n;
    }
}

void n_poly_fq_shift_left(
    n_poly_t A,
    const n_poly_t B,
    slong n,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    n_poly_fit_length(A, d*(B->length + n));
    flint_mpn_copyd(A->coeffs + d*n, B->coeffs, d*B->length);
    flint_mpn_zero(A->coeffs, d*n);
    A->length = B->length + n;
}

void n_poly_fq_truncate(n_poly_t A, slong len, const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    if (A->length > len)
    {
        A->length = len;
        _n_poly_fq_normalise(A, d);
    }
}

void n_poly_fq_evaluate_n_fq(
    mp_limb_t * e,
    const n_poly_t A,
    const mp_limb_t * c,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    mp_limb_t * u = FLINT_ARRAY_ALLOC(d, mp_limb_t);

    _n_fq_zero(e, d);
    for (i = 0; i < A->length; i++)
    {
        n_fq_pow_ui(u, c, i, ctx);
        n_fq_mul(u, u, A->coeffs + d*i, ctx);
        n_fq_add(e, e, u, ctx);
    }

    flint_free(u);
}


void n_poly_fq_set_fq_nmod_poly(
    n_poly_t A,
    const fq_nmod_poly_t B,
    const fq_nmod_ctx_t ctx)
{
    slong i;
    slong d = fq_nmod_ctx_degree(ctx);

    n_poly_fit_length(A, d*B->length);

    for (i = 0; i < B->length; i++)
        n_fq_set_fq_nmod(A->coeffs + d*i, B->coeffs + i, ctx);

    A->length = B->length;
}

void n_poly_fq_get_fq_nmod_poly(
    fq_nmod_poly_t A,
    const n_poly_t B,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;

    FLINT_ASSERT(B->alloc >= d*B->length);

    fq_nmod_poly_fit_length(A, B->length, ctx);

    for (i = 0; i < B->length; i++)
        n_fq_get_fq_nmod(A->coeffs + i, B->coeffs + d*i, ctx);

    A->length = B->length;
}


void n_poly_fq_scalar_mul_n_fq(
    n_poly_t A,
    const n_poly_t B,
    const mp_limb_t * c,
    const fq_nmod_ctx_t ctx)
{
    slong i, d = fq_nmod_ctx_degree(ctx);
    n_poly_fit_length(A, d*B->length);
    for (i = 0; i < B->length; i++)
        n_fq_mul(A->coeffs + d*i, B->coeffs + d*i, c, ctx);
    A->length = B->length;
    _n_poly_fq_normalise(A, d);
}

void n_poly_fq_make_monic(
    n_poly_t A,
    const n_poly_t B,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong itch = FLINT_MAX(N_FQ_MUL_ITCH, N_FQ_INV_ITCH);
    mp_limb_t * tmp, * inv;
    slong i, Blen = B->length;

    if (Blen < 1)
    {
        n_poly_zero(A);
        return;
    }

    n_poly_fit_length(A, d*Blen);

    tmp = FLINT_ARRAY_ALLOC(d*(itch + 1), mp_limb_t);
    inv = tmp + d*itch;

    _n_fq_inv(inv, B->coeffs + d*(Blen - 1), ctx, tmp);

    for (i = 0; i < Blen - 1; i++)
        _n_fq_mul(A->coeffs + d*i, B->coeffs + d*i, inv, ctx, tmp);

    _n_fq_one(A->coeffs + d*(Blen - 1), d);
    A->length = Blen;

    flint_free(tmp);
}


/* multiply A by (x^k - c) */
void n_poly_fq_shift_left_scalar_submul(
    n_poly_t A,
    slong k,
    const mp_limb_t * c,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    mp_limb_t * Acoeffs;
    slong i;
    slong Alen = A->length;
    mp_limb_t * u = FLINT_ARRAY_ALLOC(d, mp_limb_t);

    n_poly_fit_length(A, d*(Alen + k));

    Acoeffs = A->coeffs;

    flint_mpn_copyd(Acoeffs + d*k, Acoeffs, d*Alen);
    flint_mpn_zero(Acoeffs, d*k);

    for (i = 0; i < A->length; i++)
    {
        n_fq_mul(u, c, Acoeffs + d*(i + k), ctx);
        _n_fq_sub(Acoeffs + d*i, Acoeffs + d*i, u, ctx);
    }

    A->length = Alen + k;

    flint_free(u);
}


ulong n_poly_fq_remove(
    n_poly_t f,
    const n_poly_t g,
    const fq_nmod_ctx_t ctx)
{
    n_poly_t q, r;
    ulong i = 0;

    n_poly_init(q);
    n_poly_init(r);

    while (1)
    {
        if (f->length < g->length)
            break;
        n_poly_fq_divrem(q, r, f, g, ctx);
        if (r->length == 0)
            n_poly_swap(q, f);
        else
            break;
        i++;
    }

    n_poly_clear(q);
    n_poly_clear(r);

    return i;
}
