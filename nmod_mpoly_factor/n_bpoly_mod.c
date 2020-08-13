/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


int n_bpoly_mod_is_canonical(const n_bpoly_t A, nmod_t mod)
{
    slong i;

    if (A->length <= 0)
        return A->length == 0;

    for (i = 0; i < A->length; i++)
    {
        if (!n_poly_mod_is_canonical(A->coeffs + i, mod))
            return 0;
    }

    return !n_poly_is_zero(A->coeffs + A->length - 1);
}

void nmod_mpoly_get_bpoly(
    n_bpoly_t A,
    const nmod_mpoly_t B,
    slong varx,
    slong vary,
    const nmod_mpoly_ctx_t ctx)
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
        n_bpoly_set_coeff(A, Bexpx, Bexpy, B->coeffs[j]);
    }
}


void nmod_mpoly_set_bpoly(
    nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const n_bpoly_t B,
    slong varx,
    slong vary,
    const nmod_mpoly_ctx_t ctx)
{
    slong n = ctx->minfo->nvars;
    slong i, j;
    slong NA;
    slong Alen;
    mp_limb_t * Acoeff;
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
    nmod_mpoly_fit_bits(A, Abits, ctx);
    A->bits = Abits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        n_poly_struct * Bc = B->coeffs + i;
        _nmod_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + Bc->length, NA);

        for (j = 0; j < Bc->length; j++)
        {
            if (0 == Bc->coeffs[j])
                continue;

            Aexps[varx] = i;
            Aexps[vary] = j;
            Acoeff[Alen] = Bc->coeffs[j];
            mpoly_set_monomial_ui(Aexp + NA*Alen, Aexps, Abits, ctx->minfo);
            Alen++;
        }
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    A->length = Alen;

    TMP_END;

    nmod_mpoly_sort_terms(A, ctx);
}

void n_bpoly_mod_taylor_shift_var1(n_bpoly_t A, const n_bpoly_t B,
                                                      mp_limb_t c, nmod_t ctx)
{
    slong i;

    n_bpoly_set(A, B);

    for (i = A->length - 1; i >= 0; i--)
        n_poly_mod_taylor_shift(A->coeffs + i, c, ctx);
}

void n_bpoly_mod_sub(
    n_bpoly_t A,
    const n_bpoly_t B,
    const n_bpoly_t C,
    nmod_t mod)
{
    slong i;
    slong Alen = FLINT_MAX(B->length, C->length);

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    n_bpoly_fit_length(A, Alen);

    A->length = 0;
    for (i = 0; i < Alen; i++)
    {
        if (i < B->length)
        {
            if (i < C->length)
            {
                n_poly_mod_sub(A->coeffs + i, B->coeffs + i, C->coeffs + i, mod);
            }
            else
            {
                n_poly_set(A->coeffs + i, B->coeffs + i);
            }
        }
        else
        {
            FLINT_ASSERT(i < C->length);
            n_poly_mod_neg(A->coeffs + i, C->coeffs + i, mod);
        }

        if (!n_poly_is_zero(A->coeffs + i))
            A->length = i + 1;
    }
}


void n_bpoly_mod_add(
    n_bpoly_t A,
    const n_bpoly_t B,
    const n_bpoly_t C,
    nmod_t mod)
{
    slong i;
    slong Alen = FLINT_MAX(B->length, C->length);

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    n_bpoly_fit_length(A, Alen);

    A->length = 0;
    for (i = 0; i < Alen; i++)
    {
        if (i < B->length)
        {
            if (i < C->length)
            {
                n_poly_mod_add(A->coeffs + i, B->coeffs + i, C->coeffs + i, mod);
            }
            else
            {
                n_poly_set(A->coeffs + i, B->coeffs + i);
            }
        }
        else
        {
            FLINT_ASSERT(i < C->length);
            n_poly_set(A->coeffs + i, C->coeffs + i);
        }

        if (!n_poly_is_zero(A->coeffs + i))
            A->length = i + 1;
    }
}


void n_bpoly_mod_make_primitive(n_poly_t g, n_bpoly_t A, nmod_t ctx)
{
    mp_limb_t c = 1;
    slong Alen = A->length;
    slong i;
    n_poly_t q, r;

    n_poly_init(q);
    n_poly_init(r);

    n_poly_zero(g);
    for (i = 0; i < Alen; i++)
    {
        n_poly_mod_gcd(q, g, A->coeffs + i, ctx);
        n_poly_swap(g, q);
    }

    for (i = 0; i < Alen; i++)
    {
        n_poly_mod_divrem(q, r, A->coeffs + i, g, ctx);
        FLINT_ASSERT(n_poly_is_zero(r));
        n_poly_swap(A->coeffs + i, q);
    }

    /* make lc_xy(A) one */
    if (Alen > 0)
    {
        c = A->coeffs[Alen - 1].coeffs[A->coeffs[Alen - 1].length - 1];
        if (c != 1)
        {
            _n_poly_mod_scalar_mul_nmod(g, g, c, ctx);
            c = nmod_inv(c, ctx);
            for (i = 0; i < Alen; i++)
                _n_poly_mod_scalar_mul_nmod(A->coeffs + i,
                                            A->coeffs + i, c, ctx);
        }
    }

    n_poly_clear(q);
    n_poly_clear(r);
}


void n_bpoly_mod_mul(n_bpoly_t A, const n_bpoly_t B, const n_bpoly_t C,
                                                                    nmod_t mod)
{
    slong i, j;
    n_poly_struct * t;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    if (B->length <= 0 || C->length <= 0)
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
            _n_poly_mod_mul(t, B->coeffs + i, C->coeffs + j, mod);
            n_poly_mod_add(A->coeffs + i + j, A->coeffs + i + j, t, mod);
        }
    }

    A->length = B->length + C->length - 1;
    n_bpoly_normalise(A);
}


void n_bpoly_mod_mul_series(n_bpoly_t A, const n_bpoly_t B, const n_bpoly_t C,
                                                       slong order, nmod_t mod)
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
            n_poly_mod_mullow(t, B->coeffs + i, C->coeffs + j, order, mod);
            n_poly_mod_add(A->coeffs + i + j, A->coeffs + i + j, t, mod);
        }
    }

    A->length = B->length + C->length - 1;
    n_bpoly_normalise(A);

    n_poly_clear(t);
}


void n_bpoly_mod_derivative(
    n_bpoly_t A,
    const n_bpoly_t B,
    nmod_t ctx)
{
    slong i;

    FLINT_ASSERT(A != B);

    if (B->length < 2)
    {
        n_bpoly_zero(A);
        return;
    }

    n_bpoly_fit_length(A, B->length - 1);

    for (i = 1; i < B->length; i++)
        n_poly_mod_scalar_mul_ui(A->coeffs + i - 1, B->coeffs + i, i, ctx);

    A->length = B->length - 1;
    n_bpoly_normalise(A);
}


/*
    division in ((Z/nZ)[y]/y^order)[x]
    inputs need not be reduced mod y^order
*/
void n_bpoly_mod_divrem_series(
    n_bpoly_t Q,
    n_bpoly_t R,
    const n_bpoly_t A,
    const n_bpoly_t B,
    slong order,
    nmod_t ctx)
{
    slong i, qoff;
    n_poly_t q, t;

    FLINT_ASSERT(R != A);
    FLINT_ASSERT(R != B);
    FLINT_ASSERT(Q != A);
    FLINT_ASSERT(Q != B);
    FLINT_ASSERT(B->length > 0);

    n_poly_init(q);
    n_poly_init(t);

    n_bpoly_set(R, A);
    for (i = 0; i < R->length; i++)
        n_poly_truncate(R->coeffs + i, order);
    n_bpoly_normalise(R);

    Q->length = 0;

    while (R->length >= B->length)
    {
        n_poly_mod_div_series(q, R->coeffs + R->length - 1,
                                        B->coeffs + B->length - 1, order, ctx);

        for (i = 0; i < B->length; i++)
        {
            n_poly_mod_mullow(t, B->coeffs + i, q, order, ctx);
            n_poly_mod_sub(R->coeffs + i + R->length - B->length,
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

        n_poly_set(Q->coeffs + qoff, q);

        FLINT_ASSERT(n_poly_is_zero(R->coeffs + R->length - 1));

        n_bpoly_normalise(R);
    }

    n_poly_clear(q);
    n_poly_clear(t);
}


int n_bpoly_mod_divides(
    n_bpoly_t Q,
    const n_bpoly_t A,
    const n_bpoly_t B,
    nmod_t ctx)
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
    n_bpoly_set(R, A);

    Q->length = 0;

    while (R->length >= B->length)
    {
        _n_poly_mod_divrem(q, t, R->coeffs + R->length - 1,
                                 B->coeffs + B->length - 1, ctx);
        if (!n_poly_is_zero(t))
        {
            divides = 0;
            goto cleanup;
        }

        for (i = 0; i < B->length; i++)
        {
            _n_poly_mod_mul(t, B->coeffs + i, q, ctx);
            n_poly_mod_sub(R->coeffs + i + R->length - B->length,
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

        n_poly_set(Q->coeffs + qoff, q);

        while (R->length > 0 && n_poly_is_zero(R->coeffs + R->length - 1))
            R->length--;
    }

    divides = (R->length == 0);

cleanup:

    n_poly_clear(q);
    n_poly_clear(t);
    n_bpoly_clear(R);

    return divides;
}

void n_bpoly_mod_taylor_shift_var0(n_bpoly_t A, mp_limb_t c, nmod_t mod)
{
    slong n, i, j;
    n_poly_t t;

    FLINT_ASSERT(c < mod.n);

    if (c == 0)
        return;

    n_poly_init(t);
    n = A->length;

    for (i = n - 2; i >= 0; i--)
    {
        for (j = i; j < n - 1; j++)
        {
            _n_poly_mod_scalar_mul_nmod(t, A->coeffs + j + 1, c, mod);
            n_poly_mod_add(A->coeffs + j, A->coeffs + j, t, mod);
        }
    }

    n_poly_clear(t);
}


