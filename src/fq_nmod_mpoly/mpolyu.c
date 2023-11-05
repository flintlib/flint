/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "fq_nmod_mpoly.h"
#include "fq_nmod_mpoly_factor.h"


int fq_nmod_mpolyu_is_canonical(const fq_nmod_mpolyu_t A,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    for (i = 0; i < A->length; i++)
    {
        if ((slong)(A->exps[i]) < 0)
        {
            return 0;
        }

        if (i > 0 && A->exps[i - 1] <= A->exps[i])
        {
            return 0;
        }

        if (!fq_nmod_mpoly_is_canonical(A->coeffs + i, ctx))
        {
            return 0;
        }

        if (fq_nmod_mpoly_is_zero(A->coeffs + i, ctx))
        {
            return 0;
        }
    }
    return 1;
}

void fq_nmod_mpolyu_init(fq_nmod_mpolyu_t A, flint_bitcnt_t bits,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}


void fq_nmod_mpolyu_clear(fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t uctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fq_nmod_mpoly_clear(A->coeffs + i, uctx);
    flint_free(A->coeffs);
    flint_free(A->exps);
}

void fq_nmod_mpolyu_zero(fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t uctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++) {
        fq_nmod_mpoly_clear(A->coeffs + i, uctx);
        fq_nmod_mpoly_init(A->coeffs + i, uctx);
    }
    A->length = 0;
}

int fq_nmod_mpolyu_is_one(fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t uctx)
{
    if (A->length != 1 || A->exps[0] != UWORD(0))
        return 0;

    return fq_nmod_mpoly_is_one(A->coeffs + 0, uctx);
}

void fq_nmod_mpolyu_print_pretty(const fq_nmod_mpolyu_t poly,
                                const char ** x, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    if (poly->length == 0)
        flint_printf("0");
    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        fq_nmod_mpoly_print_pretty(poly->coeffs + i,x,ctx);
        flint_printf(")*X^%wd",poly->exps[i]);
    }
}

void fq_nmod_mpolyu_fit_length(fq_nmod_mpolyu_t A, slong length,
                                                const fq_nmod_mpoly_ctx_t uctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        A->exps = (ulong *) flint_realloc(A->exps, new_alloc*sizeof(ulong));
        A->coeffs = (fq_nmod_mpoly_struct *) flint_realloc(A->coeffs,
                                       new_alloc*sizeof(fq_nmod_mpoly_struct));

        for (i = old_alloc; i < new_alloc; i++)
            fq_nmod_mpoly_init3(A->coeffs + i, 0, A->bits, uctx);

        A->alloc = new_alloc;
    }
}

void fq_nmod_mpolyu_one(fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t uctx)
{
    fq_nmod_mpolyu_fit_length(A, 1, uctx);
    A->exps[0] = 0;
    fq_nmod_mpoly_one(A->coeffs + 0, uctx);
    A->length = 1;
}


void fq_nmod_mpolyu_degrees_si(
    slong * degs,
    const fq_nmod_mpolyu_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    flint_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    ulong * pmax, mask;
    TMP_INIT;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    if (A->length < 1)
    {
        for (j = 0; j < ctx->minfo->nvars; j++)
            degs[j] = -1;
    }

    TMP_START;

    mask = mpoly_overflow_mask_sp(bits);

    pmax = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(pmax, N);

    for (i = 0; i < A->length; i++)
    {
        ulong * Aiexps = A->coeffs[i].exps;
        FLINT_ASSERT(A->coeffs[i].bits == bits);

        for (j = 0; j < A->coeffs[i].length; j++)
            mpoly_monomial_max(pmax, pmax, Aiexps + N*j, bits, N, mask);
    }

    mpoly_unpack_vec_ui((ulong *) degs, pmax, bits, ctx->minfo->nvars, 1);

    for (i = 0; i < ctx->minfo->nvars/2; i++)
        FLINT_SWAP(slong, degs[i], degs[ctx->minfo->nvars - i - 1]);

    TMP_END;
}

void fq_nmod_mpolyu_repack_bits_inplace(
    fq_nmod_mpolyu_t A,
    flint_bitcnt_t bits,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    if (bits == A->bits)
        return;

    A->bits = bits;

    for (i = 0; i < A->alloc; i++)
        fq_nmod_mpoly_repack_bits_inplace(A->coeffs + i, bits, ctx);
}


/* if the coefficient doesn't exist, a new one is created (and set to zero) */
fq_nmod_mpoly_struct * _fq_nmod_mpolyu_get_coeff(fq_nmod_mpolyu_t A,
                                     ulong pow, const fq_nmod_mpoly_ctx_t uctx)
{
    slong i, j;
    fq_nmod_mpoly_struct * xk;

    for (i = 0; i < A->length && A->exps[i] >= pow; i++)
    {
        if (A->exps[i] == pow)
        {
            return A->coeffs + i;
        }
    }

    fq_nmod_mpolyu_fit_length(A, A->length + 1, uctx);

    for (j = A->length; j > i; j--)
    {
        A->exps[j] = A->exps[j - 1];
        fq_nmod_mpoly_swap(A->coeffs + j, A->coeffs + j - 1, uctx);
    }

    A->length++;
    A->exps[i] = pow;
    xk = A->coeffs + i;
    xk->length = 0;
    FLINT_ASSERT(xk->bits == A->bits);

    return xk;
}



/*
    Convert B to A using the variable permutation perm.
    The uctx should be the context of the coefficients of A.
    The ctx should be the context of B.

    operation on each term:

    for 0 <= k <= m
        l = perm[k]
        Aexp[k] = (Bexp[l] - shift[l])/stride[l]
*/
void fq_nmod_mpoly_to_mpolyu_perm_deflate(
    fq_nmod_mpolyu_t A,
    const fq_nmod_mpoly_ctx_t uctx,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i, j, k, l;
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong NA, NB;
    ulong * uexps;
    ulong * Bexps;
    fq_nmod_mpoly_struct * Ac;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m + 1 <= n);

    TMP_START;

    FLINT_ASSERT(fq_nmod_mpoly_is_canonical(B, ctx));

    uexps = (ulong *) TMP_ALLOC((m + 1)*sizeof(ulong));
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    fq_nmod_mpolyu_zero(A, uctx);

    NA = mpoly_words_per_exp(A->bits, uctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    for (j = 0; j < B->length; j++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*j, B->bits, ctx->minfo);
        for (k = 0; k < m + 1; k++)
        {
            l = perm[k];
            FLINT_ASSERT(stride[l] != UWORD(0));
            FLINT_ASSERT(((Bexps[l] - shift[l]) % stride[l]) == UWORD(0));
            uexps[k] = (Bexps[l] - shift[l]) / stride[l];
        }
        Ac = _fq_nmod_mpolyu_get_coeff(A, uexps[0], uctx);
        FLINT_ASSERT(Ac->bits == A->bits);

        fq_nmod_mpoly_fit_length(Ac, Ac->length + 1, uctx);
        _n_fq_set(Ac->coeffs + d*Ac->length, B->coeffs + d*j, d);
        mpoly_set_monomial_ui(Ac->exps + NA*Ac->length, uexps + 1, A->bits, uctx->minfo);
        Ac->length++;
    }

    for (i = 0; i < A->length; i++)
        fq_nmod_mpoly_sort_terms(A->coeffs + i, uctx);

    TMP_END;

    FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(A, uctx));
}

/*
    Convert B to A using the variable permutation vector perm.
    A must be constructed with bits = Abits.

    operation on each term:

        for 0 <= l < n
            Aexp[l] = shift[l]

        for 0 <= k <= m
            l = perm[k]
            Aexp[l] += scale[l]*Bexp[k]
*/
void fq_nmod_mpoly_from_mpolyu_perm_inflate(
    fq_nmod_mpoly_t A, flint_bitcnt_t Abits,
    const fq_nmod_mpoly_ctx_t ctx,
    const fq_nmod_mpolyu_t B,
    const fq_nmod_mpoly_ctx_t uctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong i, j, k, l;
    slong NA, NB;
    slong Alen;
    mp_limb_t * Acoeff;
    ulong * Aexp;
    ulong * uexps;
    ulong * Aexps;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m + 1 <= n);
    TMP_START;

    uexps = (ulong *) TMP_ALLOC((m + 1)*sizeof(ulong));
    Aexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, uctx->minfo);

    fq_nmod_mpoly_fit_length_reset_bits(A, B->length, Abits, ctx);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        fq_nmod_mpoly_struct * Bc = B->coeffs + i;
        FLINT_ASSERT(Bc->bits == B->bits);

        _fq_nmod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc, d,
                                  &Aexp, &A->exps_alloc, NA, Alen + Bc->length);

        for (j = 0; j < Bc->length; j++)
        {
            _n_fq_set(Acoeff + d*(Alen + j), Bc->coeffs + d*j, d);
            mpoly_get_monomial_ui(uexps + 1, Bc->exps + NB*j, Bc->bits, uctx->minfo);
            uexps[0] = B->exps[i];
            for (l = 0; l < n; l++)
            {
                Aexps[l] = shift[l];
            }
            for (k = 0; k < m + 1; k++)
            {
                l = perm[k];
                Aexps[l] += stride[l]*uexps[k];
            }
            mpoly_set_monomial_ui(Aexp + NA*(Alen + j), Aexps, Abits, ctx->minfo);
        }
        Alen += Bc->length;
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->length = Alen;

    TMP_END;

    fq_nmod_mpoly_sort_terms(A, ctx);
}


/*
    Convert B to A using the variable permutation perm.
    The uctx should be the context of the coefficients of A.
    The ctx should be the context of B.

    operation on each term:

    for 0 <= k < m + 2
        l = perm[k]
        Aexp[k] = (Bexp[l] - shift[l])/stride[l]

    the most significant main variable uses Aexp[0]
    the least significant main variable uses Aexp[1]
    the coefficients of A use variables Aexp[2], ..., Aexp[m + 1]
    maxexps if it exists is supposed to be a degree bound on B
*/
void fq_nmod_mpoly_to_mpolyuu_perm_deflate(
    fq_nmod_mpolyu_t A,
    const fq_nmod_mpoly_ctx_t uctx,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i, j, k, l;
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong NA, NB;
    ulong * uexps;
    ulong * Bexps;
    fq_nmod_mpoly_struct * Ac;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m + 2 <= n);

    fq_nmod_mpolyu_zero(A, uctx);

    TMP_START;

    uexps = (ulong *) TMP_ALLOC((m + 2)*sizeof(ulong));
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));
    NA = mpoly_words_per_exp(A->bits, uctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    for (j = 0; j < B->length; j++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*j, B->bits, ctx->minfo);
        for (k = 0; k < m + 2; k++)
        {
            l = perm[k];
            if (stride[l] == 1)
            {
                uexps[k] = (Bexps[l] - shift[l]);
            }
            else
            {
                FLINT_ASSERT(stride[l] != 0);
                FLINT_ASSERT(((Bexps[l] - shift[l]) % stride[l]) == 0);
                uexps[k] = (Bexps[l] - shift[l]) / stride[l];
            }
        }
        FLINT_ASSERT(FLINT_BIT_COUNT(uexps[0]) < FLINT_BITS/2);
        FLINT_ASSERT(FLINT_BIT_COUNT(uexps[1]) < FLINT_BITS/2);
        Ac = _fq_nmod_mpolyu_get_coeff(A, (uexps[0] << (FLINT_BITS/2)) + uexps[1], uctx);
        FLINT_ASSERT(Ac->bits == A->bits);

        fq_nmod_mpoly_fit_length(Ac, Ac->length + 1, uctx);
        _n_fq_set(Ac->coeffs + d*Ac->length, B->coeffs + d*j, d);
        mpoly_set_monomial_ui(Ac->exps + NA*Ac->length, uexps + 2, A->bits, uctx->minfo);
        Ac->length++;
    }

    for (i = 0; i < A->length; i++)
        fq_nmod_mpoly_sort_terms(A->coeffs + i, uctx);

    TMP_END;
}


/*
    Convert B to A using the variable permutation vector perm.
    A must be constructed with bits = Abits.

    operation on each term:

        for 0 <= l < n
            Aexp[l] = shift[l]

        for 0 <= k < m + 2
            l = perm[k]
            Aexp[l] += scale[l]*Bexp[k]
*/
void fq_nmod_mpoly_from_mpolyuu_perm_inflate( /* only for 2 main vars */
    fq_nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fq_nmod_mpoly_ctx_t ctx,
    const fq_nmod_mpolyu_t B,
    const fq_nmod_mpoly_ctx_t uctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong i, j, k, l;
    slong NA, NB;
    slong Alen;
    mp_limb_t * Acoeff;
    ulong * Aexp;
    ulong * uexps;
    ulong * Aexps;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m + 2 <= n);
    TMP_START;

    uexps = (ulong *) TMP_ALLOC((m + 2)*sizeof(ulong));
    Aexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, uctx->minfo);

    fq_nmod_mpoly_fit_length_reset_bits(A, B->length, Abits, ctx);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        fq_nmod_mpoly_struct * Bc = B->coeffs + i;
        FLINT_ASSERT(Bc->bits == B->bits);

        _fq_nmod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc, d,
                                  &Aexp, &A->exps_alloc, NA, Alen + Bc->length);

        for (j = 0; j < Bc->length; j++)
        {
            _n_fq_set(Acoeff + d*(Alen + j), Bc->coeffs + d*j, d);
            mpoly_get_monomial_ui(uexps + 2, Bc->exps + NB*j, Bc->bits, uctx->minfo);
            uexps[0] = B->exps[i] >> (FLINT_BITS/2);
            uexps[1] = B->exps[i] & ((-UWORD(1)) >> (FLINT_BITS - FLINT_BITS/2));
            for (l = 0; l < n; l++)
            {
                Aexps[l] = shift[l];
            }
            for (k = 0; k < m + 2; k++)
            {
                l = perm[k];
                Aexps[l] += stride[l]*uexps[k];
            }
            mpoly_set_monomial_ui(Aexp + NA*(Alen + j), Aexps, Abits, ctx->minfo);
        }
        Alen += Bc->length;
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    _fq_nmod_mpoly_set_length(A, Alen, ctx);

    fq_nmod_mpoly_sort_terms(A, ctx);
    TMP_END;
}




/** 0 variables *************************************/


void fq_nmod_mpoly_to_mpolyl_perm_deflate(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t lctx,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
    slong d;
    slong j, k, l;
    slong m = lctx->minfo->nvars;
    slong n = ctx->minfo->nvars;
    slong NA, NB;
    ulong * lexps;
    ulong * Bexps;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m <= n);
    TMP_START;

    fq_nmod_mpoly_fit_length(A, B->length, lctx);
    A->length = B->length;

    lexps = (ulong *) TMP_ALLOC(m*sizeof(ulong));
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));
    NA = mpoly_words_per_exp(A->bits, lctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    d = fq_nmod_ctx_degree(ctx->fqctx);
    FLINT_ASSERT(d == fq_nmod_ctx_degree(lctx->fqctx));
    _nmod_vec_set(A->coeffs, B->coeffs, B->length*d);

    for (j = 0; j < B->length; j++)
    {
        mpoly_get_monomial_ui(Bexps, B->exps + NB*j, B->bits, ctx->minfo);
        for (k = 0; k < m; k++)
        {
            l = perm[k];
            if (stride[l] == 1)
            {
                lexps[k] = (Bexps[l] - shift[l]);
            }
            else
            {
                FLINT_ASSERT(stride[l] != 0);
                FLINT_ASSERT(((Bexps[l] - shift[l]) % stride[l]) == 0);
                lexps[k] = (Bexps[l] - shift[l]) / stride[l];
            }
        }

        mpoly_set_monomial_ui(A->exps + NA*j, lexps, A->bits, lctx->minfo);
    }

    TMP_END;

    fq_nmod_mpoly_sort_terms(A, lctx);

    FLINT_ASSERT(fq_nmod_mpoly_is_canonical(A, lctx));
}


/*
    Convert B to A using the variable permutation vector perm.
    A must be constructed with bits = Abits.

    operation on each term:

        for 0 <= l < n
            Aexp[l] = shift[l]

        for 0 <= k < m + 2
            l = perm[k]
            Aexp[l] += scale[l]*Bexp[k]
*/
void fq_nmod_mpoly_from_mpolyl_perm_inflate(
    fq_nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fq_nmod_mpoly_ctx_t ctx,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t lctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
    slong d;
    slong n = ctx->minfo->nvars;
    slong m = lctx->minfo->nvars;
    slong i, k, l;
    slong NA, NB;
    ulong * Bexps;
    ulong * Aexps;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m <= n);
    TMP_START;

    Bexps = (ulong *) TMP_ALLOC(m*sizeof(ulong));
    Aexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, lctx->minfo);

    fq_nmod_mpoly_fit_length_reset_bits(A, B->length, Abits, ctx);
    A->length = B->length;

    d = fq_nmod_ctx_degree(ctx->fqctx);
    FLINT_ASSERT(d == fq_nmod_ctx_degree(lctx->fqctx));
    _nmod_vec_set(A->coeffs, B->coeffs, B->length*d);

    for (i = 0; i < B->length; i++)
    {
	    mpoly_get_monomial_ui(Bexps, B->exps + NB*i, B->bits, lctx->minfo);

        for (l = 0; l < n; l++)
        {
            Aexps[l] = shift[l];
        }
        for (k = 0; k < m; k++)
        {
            l = perm[k];
            Aexps[l] += stride[l]*Bexps[k];
        }

        mpoly_set_monomial_ui(A->exps + NA*i, Aexps, Abits, ctx->minfo);
    }

    TMP_END;

    fq_nmod_mpoly_sort_terms(A, ctx);

    FLINT_ASSERT(fq_nmod_mpoly_is_canonical(A, ctx));
}

void fq_nmod_mpolyu_shift_right(fq_nmod_mpolyu_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(A->exps[i] >= s);
        A->exps[i] -= s;
    }
}

void fq_nmod_mpolyu_shift_left(fq_nmod_mpolyu_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        A->exps[i] += s;
    }
}

void fq_nmod_mpolyu_scalar_mul_fq_nmod(fq_nmod_mpolyu_t A, fq_nmod_t c,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
        fq_nmod_mpoly_scalar_mul_fq_nmod(A->coeffs + i, A->coeffs + i, c, ctx);
}


void fq_nmod_mpolyu_set(fq_nmod_mpolyu_t A, const fq_nmod_mpolyu_t B,
                                                const fq_nmod_mpoly_ctx_t uctx)
{
    slong i;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(fq_nmod_mpolyu_is_canonical(B, uctx));

    fq_nmod_mpolyu_fit_length(A, B->length, uctx);
    A->length = B->length;

    for (i = 0; i < B->length; i++)
    {
        fq_nmod_mpoly_set(A->coeffs + i, B->coeffs + i, uctx);
        A->exps[i] = B->exps[i];
    }
}


void fq_nmod_mpolyu_divexact_mpoly_inplace(
    fq_nmod_mpolyu_t A,
    fq_nmod_mpoly_t c,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    flint_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    ulong * cmpmask;
    fq_nmod_mpoly_t t;
    TMP_INIT;

    FLINT_ASSERT(bits == c->bits);
    FLINT_ASSERT(c->length > 0);

    if (fq_nmod_mpoly_is_fq_nmod(c, ctx))
    {
        slong d = fq_nmod_ctx_degree(ctx->fqctx);
        mp_limb_t * inv;

        if (_n_fq_is_one(c->coeffs + d*0, d))
            return;

        TMP_START;

        inv = (mp_limb_t *) TMP_ALLOC(d*sizeof(mp_limb_t));

        n_fq_inv(inv, c->coeffs + d*0, ctx->fqctx);

        for (i = 0; i < A->length; i++)
        {
            fq_nmod_mpoly_struct * Ai = A->coeffs + i;
            for (j = 0; j < Ai->length; j++)
                n_fq_mul(Ai->coeffs + d*j, Ai->coeffs + d*j, inv, ctx->fqctx);
        }

        TMP_END;

        return;
    }

    fq_nmod_mpoly_init3(t, 0, bits, ctx);

    TMP_START;

    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    for (i = A->length - 1; i >= 0; i--)
    {
        FLINT_ASSERT(A->coeffs[i].bits == bits);
        _fq_nmod_mpoly_divides_monagan_pearce(t, A->coeffs[i].coeffs,
                   A->coeffs[i].exps, A->coeffs[i].length, c->coeffs, c->exps,
                                      c->length, bits, N, cmpmask, ctx->fqctx);
        fq_nmod_mpoly_swap(A->coeffs + i, t, ctx);
        FLINT_ASSERT(A->coeffs[i].length > 0);
    }

    TMP_END;

    fq_nmod_mpoly_clear(t, ctx);
}


/*
    A = B * c and preserve the bit packing
*/
void fq_nmod_mpolyu_mul_mpoly(
    fq_nmod_mpolyu_t A,
    fq_nmod_mpolyu_t B,
    fq_nmod_mpoly_t c,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    flint_bitcnt_t bits = B->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    ulong * cmpmask;
    TMP_INIT;

    TMP_START;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == c->bits);

    fq_nmod_mpolyu_fit_length(A, B->length, ctx);

    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    for (i = 0; i < B->length; i++)
    {
        A->exps[i] = B->exps[i];
        fq_nmod_mpoly_fit_length(A->coeffs + i, B->coeffs[i].length + c->length, ctx);
        _fq_nmod_mpoly_mul_johnson(A->coeffs + i, c->coeffs, c->exps, c->length,
                 B->coeffs[i].coeffs, B->coeffs[i].exps, B->coeffs[i].length,
                                                 bits, N, cmpmask, ctx->fqctx);
    }
    A->length = B->length;

    TMP_END;
}


void fq_nmod_mpolyu_mul_mpoly_inplace(
    fq_nmod_mpolyu_t A,
    fq_nmod_mpoly_t c,
    const fq_nmod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    slong i;
    ulong * cmpmask;
    fq_nmod_mpoly_t t;
    TMP_INIT;

    FLINT_ASSERT(bits == c->bits);
    FLINT_ASSERT(c->length > 0);

    if (fq_nmod_mpoly_is_fq_nmod(c, ctx))
    {
        if (n_fq_is_one(c->coeffs + 0, ctx->fqctx))
            return;

        for (i = 0; i < A->length; i++)
            fq_nmod_mpoly_scalar_mul_n_fq(A->coeffs + i, A->coeffs + i,
                                                           c->coeffs + 0, ctx);
        return;
    }

    fq_nmod_mpoly_init3(t, 0, bits, ctx);

    TMP_START;

    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    for (i = A->length - 1; i >= 0; i--)
    {
        FLINT_ASSERT(A->coeffs[i].bits == bits);
        _fq_nmod_mpoly_mul_johnson(t, A->coeffs[i].coeffs,
                  A->coeffs[i].exps, A->coeffs[i].length, c->coeffs, c->exps,
                                      c->length, bits, N, cmpmask, ctx->fqctx);
        fq_nmod_mpoly_swap(A->coeffs + i, t, ctx);
        FLINT_ASSERT(A->coeffs[i].length > 0);
    }

    TMP_END;

    fq_nmod_mpoly_clear(t, ctx);
}


int fq_nmod_mpolyu_content_mpoly(
    fq_nmod_mpoly_t g,
    const fq_nmod_mpolyu_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;

    success = _fq_nmod_mpoly_vec_content_mpoly(g, A->coeffs, A->length, ctx);
    if (!success)
        fq_nmod_mpoly_zero(g, ctx);

    fq_nmod_mpoly_repack_bits_inplace(g, A->bits, ctx);

    return success;
}
