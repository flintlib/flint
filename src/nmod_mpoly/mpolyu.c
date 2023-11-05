/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "nmod_mpoly_factor.h"

void nmod_mpolyu_init(nmod_mpolyu_t A, flint_bitcnt_t bits,
                                                    const nmod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}


void nmod_mpolyu_clear(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        nmod_mpoly_clear(A->coeffs + i, uctx);
    flint_free(A->coeffs);
    flint_free(A->exps);
}


void nmod_mpolyu_swap(nmod_mpolyu_t A, nmod_mpolyu_t B,
                                                   const nmod_mpoly_ctx_t uctx)
{
   nmod_mpolyu_struct t = *A;
   *A = *B;
   *B = t;
}

void nmod_mpolyu_zero(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx)
{
    A->length = 0;
}

int nmod_mpolyu_is_one(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx)
{
    if (A->length != 1 || A->exps[0] != UWORD(0))
        return 0;

    return nmod_mpoly_is_one(A->coeffs + 0, uctx);
}

void nmod_mpolyu_print_pretty(const nmod_mpolyu_t poly,
                                   const char ** x, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    if (poly->length == 0)
        flint_printf("0");
    for (i = 0; i < poly->length; i++)
    {
        if (i != 0)
            flint_printf(" + ");
        flint_printf("(");
        nmod_mpoly_print_pretty(poly->coeffs + i,x,ctx);
        flint_printf(")*X^%wd",poly->exps[i]);
    }
}

void nmod_mpolyu_fit_length(nmod_mpolyu_t A, slong length,
                                                   const nmod_mpoly_ctx_t uctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        A->exps = (ulong *) flint_realloc(A->exps, new_alloc*sizeof(ulong));
        A->coeffs = (nmod_mpoly_struct *) flint_realloc(A->coeffs,
                                          new_alloc*sizeof(nmod_mpoly_struct));

        for (i = old_alloc; i < new_alloc; i++)
            nmod_mpoly_init3(A->coeffs + i, 0, A->bits, uctx);

        A->alloc = new_alloc;
    }
}

void nmod_mpolyu_one(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx)
{
    nmod_mpolyu_fit_length(A, WORD(1), uctx);
    A->exps[0] = UWORD(0);
    nmod_mpoly_one(A->coeffs + 0, uctx);
    A->length = WORD(1);
}


void nmod_mpolyu_degrees_si(
    slong * degs,
    const nmod_mpolyu_t A,
    const nmod_mpoly_ctx_t ctx)
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

void nmod_mpolyu_repack_bits_inplace(
    nmod_mpolyu_t A,
    flint_bitcnt_t bits,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;

    if (bits == A->bits)
        return;

    A->bits = bits;

    for (i = 0; i < A->alloc; i++)
        nmod_mpoly_repack_bits_inplace(A->coeffs + i, bits, ctx);
}

/* if the coefficient doesn't exist, a new one is created (and set to zero) */
nmod_mpoly_struct * _nmod_mpolyu_get_coeff(nmod_mpolyu_t A,
                                        ulong pow, const nmod_mpoly_ctx_t uctx)
{
    slong i, j;
    nmod_mpoly_struct * xk;

    for (i = 0; i < A->length && A->exps[i] >= pow; i++)
    {
        if (A->exps[i] == pow)
        {
            return A->coeffs + i;
        }
    }

    nmod_mpolyu_fit_length(A, A->length + 1, uctx);

    for (j = A->length; j > i; j--)
    {
        A->exps[j] = A->exps[j - 1];
        nmod_mpoly_swap(A->coeffs + j, A->coeffs + j - 1, uctx);
    }

    A->length++;
    A->exps[i] = pow;
    xk = A->coeffs + i;
    xk->length = 0;
    FLINT_ASSERT(xk->bits == A->bits);

    return xk;
}


typedef struct
{
    volatile slong index;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
    slong length;
    nmod_mpoly_struct * coeffs;
    const nmod_mpoly_ctx_struct * ctx;
}
_sort_arg_struct;

typedef _sort_arg_struct _sort_arg_t[1];


static void _worker_sort(void * varg)
{
    _sort_arg_struct * arg = (_sort_arg_struct *) varg;
    slong i;

get_next_index:

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&arg->mutex);
#endif
    i = arg->index;
    arg->index++;
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&arg->mutex);
#endif

    if (i >= arg->length)
        goto cleanup;

    nmod_mpoly_sort_terms(arg->coeffs + i, arg->ctx);

    goto get_next_index;

cleanup:

    return;
}


/** 1 variables *************************************/


/*
    Convert B to A using the variable permutation perm.
    The uctx (m vars) should be the context of the coefficients of A.
    The ctx (n vars) should be the context of B.

    operation on each term:

    for 0 <= k < m + 1
        l = perm[k]
        Aexp[k] = (Bexp[l] - shift[l])/stride[l]

    the most significant main variable uses k = 0
    the coefficients of A use variables k = 1 ... m
*/
void nmod_mpoly_to_mpolyu_perm_deflate_threaded_pool(
    nmod_mpolyu_t A,
    const nmod_mpoly_ctx_t uctx,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong i, j, k, l;
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong NA, NB;
    ulong * uexps;
    ulong * Bexps;
    nmod_mpoly_struct * Ac;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m + 1 <= n);

    TMP_START;

    uexps = (ulong *) TMP_ALLOC((m + 1)*sizeof(ulong));
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    nmod_mpolyu_zero(A, uctx);

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
        Ac = _nmod_mpolyu_get_coeff(A, uexps[0], uctx);
        FLINT_ASSERT(Ac->bits == A->bits);

        nmod_mpoly_fit_length(Ac, Ac->length + 1, uctx);
        Ac->coeffs[Ac->length] = B->coeffs[j];
        mpoly_set_monomial_ui(Ac->exps + NA*Ac->length, uexps + 1, A->bits, uctx->minfo);
        Ac->length++;
    }

    if (num_handles > 0)
    {
        _sort_arg_t arg;

#if FLINT_USES_PTHREAD
	pthread_mutex_init(&arg->mutex, NULL);
#endif
	arg->index = 0;
        arg->coeffs = A->coeffs;
        arg->length = A->length;
        arg->ctx = uctx;

        for (i = 0; i < num_handles; i++)
        {
            thread_pool_wake(global_thread_pool, handles[i], 0, _worker_sort, arg);
        }
        _worker_sort(arg);
        for (i = 0; i < num_handles; i++)
        {
            thread_pool_wait(global_thread_pool, handles[i]);
        }

#if FLINT_USES_PTHREAD
        pthread_mutex_destroy(&arg->mutex);
#endif
    }
    else
    {
        for (i = 0; i < A->length; i++)
        {
            nmod_mpoly_sort_terms(A->coeffs + i, uctx);
        }
    }

    TMP_END;
}

/*
    Convert B to A using the variable permutation vector perm.
    This function inverts nmod_mpoly_to_mpolyu_perm_deflate.
    A will be constructed with bits = Abits.

    operation on each term:

        for 0 <= l < n
            Aexp[l] = shift[l]

        for 0 <= k < m + 1
            l = perm[k]
            Aexp[l] += scale[l]*Bexp[k]
*/
void nmod_mpoly_from_mpolyu_perm_inflate(
    nmod_mpoly_t A, flint_bitcnt_t Abits,
    const nmod_mpoly_ctx_t ctx,
    const nmod_mpolyu_t B,
    const nmod_mpoly_ctx_t uctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
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

    nmod_mpoly_fit_length_reset_bits(A, B->length, Abits, ctx);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        nmod_mpoly_struct * Bc = B->coeffs + i;
        FLINT_ASSERT(Bc->bits == B->bits);

        _nmod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                               &Aexp, &A->exps_alloc, NA, Alen + Bc->length);

        for (j = 0; j < Bc->length; j++)
        {
            Acoeff[Alen + j] = Bc->coeffs[j];
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
    _nmod_mpoly_set_length(A, Alen, ctx);

    TMP_END;

    nmod_mpoly_sort_terms(A, ctx);
}


/** 0 variables *************************************/


void nmod_mpoly_to_mpolyl_perm_deflate(
    nmod_mpoly_t A,
    const nmod_mpoly_ctx_t lctx,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
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

    nmod_mpoly_fit_length(A, B->length, ctx);
    A->length = B->length;

    lexps = (ulong *) TMP_ALLOC(m*sizeof(ulong));
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));
    NA = mpoly_words_per_exp(A->bits, lctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    for (j = 0; j < B->length; j++)
    {
        A->coeffs[j] = B->coeffs[j];

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

    nmod_mpoly_sort_terms(A, lctx);

    FLINT_ASSERT(nmod_mpoly_is_canonical(A, lctx));
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
void nmod_mpoly_from_mpolyl_perm_inflate(
    nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const nmod_mpoly_ctx_t ctx,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t lctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
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

    nmod_mpoly_fit_length_reset_bits(A, B->length, Abits, ctx);
    A->length = B->length;

    for (i = 0; i < B->length; i++)
    {
	    A->coeffs[i] = B->coeffs[i];
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

    nmod_mpoly_sort_terms(A, ctx);

    FLINT_ASSERT(nmod_mpoly_is_canonical(A, ctx));
}



/** 2 variables *************************************/

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
void nmod_mpoly_to_mpolyuu_perm_deflate_threaded_pool(
    nmod_mpolyu_t A,
    const nmod_mpoly_ctx_t uctx,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride,
    const ulong * maxexps, /* nullptr is ok */
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong i, j, k, l;
    slong n = ctx->minfo->nvars;
    slong m = uctx->minfo->nvars;
    slong NA, NB;
    ulong * uexps;
    ulong * Bexps;
    nmod_mpoly_struct * Ac;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(m + 2 <= n);

    nmod_mpolyu_zero(A, uctx);

    {
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
            Ac = _nmod_mpolyu_get_coeff(A, (uexps[0] << (FLINT_BITS/2)) + uexps[1], uctx);
            FLINT_ASSERT(Ac->bits == A->bits);

            nmod_mpoly_fit_length(Ac, Ac->length + 1, uctx);
            Ac->coeffs[Ac->length] = B->coeffs[j];
            mpoly_set_monomial_ui(Ac->exps + NA*Ac->length, uexps + 2, A->bits, uctx->minfo);
            Ac->length++;
        }


            for (i = 0; i < A->length; i++)
            {
                nmod_mpoly_sort_terms(A->coeffs + i, uctx);
            }

        TMP_END;
    }
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
void nmod_mpoly_from_mpolyuu_perm_inflate( /* only for 2 main vars */
    nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const nmod_mpoly_ctx_t ctx,
    const nmod_mpolyu_t B,
    const nmod_mpoly_ctx_t uctx,
    const slong * perm,
    const ulong * shift,
    const ulong * stride)
{
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

    nmod_mpoly_fit_length_reset_bits(A, B->length, Abits, ctx);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        nmod_mpoly_struct * Bc = B->coeffs + i;
        FLINT_ASSERT(Bc->bits == B->bits);

        _nmod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                               &Aexp, &A->exps_alloc, NA, Alen + Bc->length);

        for (j = 0; j < Bc->length; j++)
        {
            Acoeff[Alen + j] = Bc->coeffs[j];
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
    _nmod_mpoly_set_length(A, Alen, ctx);

    nmod_mpoly_sort_terms(A, ctx);
    TMP_END;
}



void nmod_mpolyu_shift_right(nmod_mpolyu_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(A->exps[i] >= s);
        A->exps[i] -= s;
    }
}

void nmod_mpolyu_shift_left(nmod_mpolyu_t A, ulong s)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        A->exps[i] += s;
    }
}

void nmod_mpolyu_scalar_mul_nmod(nmod_mpolyu_t A, mp_limb_t c,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < (A->coeffs + i)->length; j++)
        {
            A->coeffs[i].coeffs[j] = nmod_mul(A->coeffs[i].coeffs[j], c, ctx->mod);
        }
    }
}


void nmod_mpolyu_set(nmod_mpolyu_t A, const nmod_mpolyu_t B,
                                                   const nmod_mpoly_ctx_t uctx)
{
    slong i;
    nmod_mpoly_struct * Acoeff, * Bcoeff;
    ulong * Aexp, * Bexp;
    slong Alen, Blen;

    Alen = 0;
    Blen = B->length;
    nmod_mpolyu_fit_length(A, Blen, uctx);
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Aexp = A->exps;
    Bexp = B->exps;

    for (i = 0; i < Blen; i++)
    {
        nmod_mpoly_set(Acoeff + Alen, Bcoeff + i, uctx);
        Aexp[Alen++] = Bexp[i];
    }
    Alen = Blen;

    /* demote remaining coefficients */
    for (i = Alen; i < A->length; i++)
    {
        nmod_mpoly_clear(Acoeff + i, uctx);
        nmod_mpoly_init(Acoeff + i, uctx);
    }
    A->length = Alen;
}



void nmod_mpoly_cvtfrom_poly_notmain(nmod_mpoly_t A, nmod_poly_t a,
                                         slong var, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong k;
    ulong * oneexp;
    slong N;
    TMP_INIT;
    TMP_START;

    FLINT_ASSERT(A->bits <= FLINT_BITS);

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);

    oneexp = (ulong *)TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_monomial_sp(oneexp, var, A->bits, ctx->minfo);

    nmod_mpoly_fit_length(A, nmod_poly_length(a), ctx);

    k = 0;
    for (i = nmod_poly_length(a) - 1; i >= 0; i--)
    {
        mp_limb_t c = nmod_poly_get_coeff_ui(a, i);
        if (c != UWORD(0))
        {
            A->coeffs[k] = c;
            mpoly_monomial_mul_ui(A->exps + N*k, oneexp, N, i);
            k++;
        }
    }
    A->length = k;
    TMP_END;
}
/*
    Set "A" to "a" where "a" is a polynomial in a non-main variable "var"
*/
void nmod_mpolyu_cvtfrom_poly_notmain(nmod_mpolyu_t A, nmod_poly_t a,
                                         slong var, const nmod_mpoly_ctx_t ctx)
{
    nmod_mpolyu_fit_length(A, 1, ctx);
    A->exps[0] = 0;
    nmod_mpoly_cvtfrom_poly_notmain(A->coeffs + 0, a, var, ctx);
    A->length = !nmod_mpoly_is_zero(A->coeffs + 0, ctx);
}



/*
    Assuming that "A" depends only on the main variable,
    convert it to a poly "a".
*/
void nmod_mpolyu_cvtto_poly(nmod_poly_t a, nmod_mpolyu_t A,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    nmod_poly_zero(a);
    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT((A->coeffs + i)->length == 1);
        FLINT_ASSERT(mpoly_monomial_is_zero((A->coeffs + i)->exps,
                      mpoly_words_per_exp((A->coeffs + i)->bits, ctx->minfo)));
        nmod_poly_set_coeff_ui(a, A->exps[i], (A->coeffs + i)->coeffs[0]);
    }
}

/*
    Convert a poly "a" to "A" in the main variable,
*/
void nmod_mpolyu_cvtfrom_poly(nmod_mpolyu_t A, nmod_poly_t a,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong k;
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

    nmod_mpolyu_zero(A, ctx);
    k = 0;
    for (i = nmod_poly_length(a) - 1; i >= 0; i--)
    {
        mp_limb_t c = nmod_poly_get_coeff_ui(a, i);
        if (c != UWORD(0))
        {
            nmod_mpolyu_fit_length(A, k + 1, ctx);
            A->exps[k] = i;
            nmod_mpoly_fit_length_reset_bits(A->coeffs + k, 1, A->bits, ctx);
            A->coeffs[k].coeffs[0] = c;
            A->coeffs[k].length = 1;
            mpoly_monomial_zero(A->coeffs[k].exps + N*0, N);
            k++;
        }
    }
    A->length = k;
}


void nmod_mpolyu_msub(nmod_mpolyu_t R, nmod_mpolyu_t A, nmod_mpolyu_t B,
                           nmod_mpoly_t c, slong e, const nmod_mpoly_ctx_t ctx)
{
    slong i, j, k;
    nmod_mpoly_t T;

    nmod_mpolyu_fit_length(R, A->length + B->length, ctx);

    nmod_mpoly_init(T, ctx);

    i = j = k = 0;
    while (i < A->length || j < B->length)
    {
        if (i < A->length && (j >= B->length || A->exps[i] > B->exps[j] + e))
        {
            /* only A ok */
            nmod_mpoly_set(R->coeffs + k, A->coeffs + i, ctx);
            R->exps[k] = A->exps[i];
            k++;
            i++;
        }
        else if (j < B->length && (i >= A->length || B->exps[j] + e > A->exps[i]))
        {
            /* only B ok */
            nmod_mpoly_mul(R->coeffs + k, B->coeffs + j, c, ctx);
            nmod_mpoly_neg(R->coeffs + k, R->coeffs + k, ctx);
            R->exps[k] = B->exps[j] + e;
            k++;
            j++;
        }
        else if (i < A->length && j < B->length && (A->exps[i] == B->exps[j] + e))
        {
            nmod_mpoly_mul(T, B->coeffs + j, c, ctx);
            nmod_mpoly_sub(R->coeffs + k, A->coeffs + i, T, ctx);
            R->exps[k] = A->exps[i];
            k += !nmod_mpoly_is_zero(R->coeffs + k, ctx);
            i++;
            j++;
        } else
        {
            FLINT_ASSERT(0);
        }
    }

    nmod_mpoly_clear(T, ctx);
    R->length = k;
}


void nmod_mpolyu_divexact_mpoly_inplace(nmod_mpolyu_t A, nmod_mpoly_t c,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    flint_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    ulong * cmpmask;
    nmod_mpoly_t t;
    TMP_INIT;

    FLINT_ASSERT(bits == c->bits);
    FLINT_ASSERT(c->length > 0);

    if (nmod_mpoly_is_ui(c, ctx))
    {
        if (c->coeffs[0] == 1)
            return;

        for (i = 0; i < A->length; i++)
        {
            nmod_mpoly_struct * Ai = A->coeffs + i;
            _nmod_vec_scalar_mul_nmod(Ai->coeffs, Ai->coeffs, Ai->length,
                                   nmod_inv(c->coeffs[0], ctx->mod), ctx->mod);
        }

        return;
    }

    nmod_mpoly_init3(t, 0, bits, ctx);

    TMP_START;

    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    for (i = A->length - 1; i >= 0; i--)
    {
        FLINT_ASSERT(A->coeffs[i].bits == bits);
        _nmod_mpoly_divides_monagan_pearce(t, A->coeffs[i].coeffs,
                   A->coeffs[i].exps, A->coeffs[i].length, c->coeffs, c->exps,
                                        c->length, bits, N, cmpmask, ctx->mod);
        nmod_mpoly_swap(A->coeffs + i, t, ctx);
        FLINT_ASSERT(A->coeffs[i].length > 0);
    }

    TMP_END;

    nmod_mpoly_clear(t, ctx);
}


/*
    A = B * c and preserve the bit packing
*/
void nmod_mpolyu_mul_mpoly(
    nmod_mpolyu_t A,
    nmod_mpolyu_t B,
    nmod_mpoly_t c,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    flint_bitcnt_t bits = B->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    ulong * cmpmask;
    TMP_INIT;

    TMP_START;

    bits = B->bits;
    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(A->bits == c->bits);
    FLINT_ASSERT(A != B);

    nmod_mpolyu_fit_length(A, B->length, ctx);

    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    for (i = 0; i < B->length; i++)
    {
        nmod_mpoly_fit_length(A->coeffs + i,
                                     B->coeffs[i].length + c->length + 1, ctx);
        _nmod_mpoly_mul_johnson(A->coeffs + i, B->coeffs[i].coeffs,
                       B->coeffs[i].exps, B->coeffs[i].length, c->coeffs,
                               c->exps, c->length, bits, N, cmpmask, ctx->mod);
        FLINT_ASSERT(A->coeffs[i].length > 0);
        A->exps[i] = B->exps[i];
    }
    A->length = B->length;
    TMP_END;
}


void nmod_mpolyu_mul_mpoly_inplace(nmod_mpolyu_t A, nmod_mpoly_t c,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    flint_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    ulong * cmpmask;
    nmod_mpoly_t t;
    TMP_INIT;

    FLINT_ASSERT(bits == c->bits);
    FLINT_ASSERT(c->length > 0);

    if (nmod_mpoly_is_ui(c, ctx))
    {
        if (c->coeffs[0] == 1)
            return;

        for (i = 0; i < A->length; i++)
        {
            nmod_mpoly_struct * Ai = A->coeffs + i;
            _nmod_vec_scalar_mul_nmod(Ai->coeffs, Ai->coeffs, Ai->length,
                                                       c->coeffs[0], ctx->mod);
        }

        return;
    }

    nmod_mpoly_init3(t, 0, bits, ctx);

    N = mpoly_words_per_exp(bits, ctx->minfo);

    TMP_START;

    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    for (i = A->length - 1; i >= 0; i--)
    {
        FLINT_ASSERT(A->coeffs[i].bits == bits);
        _nmod_mpoly_mul_johnson(t, A->coeffs[i].coeffs,
                   A->coeffs[i].exps, A->coeffs[i].length, c->coeffs, c->exps,
                                        c->length, bits, N, cmpmask, ctx->mod);
        nmod_mpoly_swap(A->coeffs + i, t, ctx);
        FLINT_ASSERT(A->coeffs[i].length > 0);
    }

    TMP_END;

    nmod_mpoly_clear(t, ctx);
}



int nmod_mpolyu_content_mpoly(
    nmod_mpoly_t g,
    const nmod_mpolyu_t A,
    const nmod_mpoly_ctx_t ctx)
{
    int success;

    success = _nmod_mpoly_vec_content_mpoly(g, A->coeffs, A->length, ctx);
    if (!success)
        nmod_mpoly_zero(g, ctx);

    nmod_mpoly_repack_bits_inplace(g, A->bits, ctx);

    return success;
}
