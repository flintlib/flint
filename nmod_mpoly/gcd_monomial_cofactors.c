/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

/* low priority TODO: these both can be sped up when A->bits == B->bits */

/*
    scan to check if the cofactors G/A and G/B are monomials
    If so, compute G and return 1. Otherwise, return 0.
    This special internal version takes info filled in as part of the 
    general gcd routine. It also assumes that A and B have the same length
    and are packed into <= FLINT_BITS.

    The function should pack its answer into bits = Gbits <= FLINT_BITS
*/
int _nmod_mpoly_gcd_monomial_cofactors_sp(nmod_mpoly_t G, flint_bitcnt_t Gbits,
         const nmod_mpoly_t A, const ulong * Amax_exp , const ulong * Amin_exp,
         const nmod_mpoly_t B, const ulong * Bmax_exp , const ulong * Bmin_exp,
                                                    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    slong NA, NB, NG;
    slong nvars = ctx->minfo->nvars;
    ulong * Aexp, * Bexp;
    mp_limb_t a0inv;
    nmod_mpoly_t T;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->length == B->length);
    FLINT_ASSERT(Gbits <= FLINT_BITS);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);

    for (j = 0; j < nvars; j++)
    {
        if (Amax_exp[j] - Amin_exp[j] != Bmax_exp[j] - Bmin_exp[j])
        {
            return 0;
        }
    }

    TMP_START;

    Aexp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
    Bexp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));

    NA = mpoly_words_per_exp(A->bits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    for (i = 0; i < A->length; i++)
    {
        ulong r, p0, p1, q0, q1, d0, d1, d2;
        /* check if A[0]*B[i] == B[0]*A[i] */
        umul_ppmm(p1, p0, A->coeffs[0], B->coeffs[i]);
        umul_ppmm(q1, q0, B->coeffs[0], A->coeffs[i]);
        sub_dddmmmsss(d2, d1, d0, 0, p1, p0, 0, q1, q0);
        FLINT_ASSERT(d2 == UWORD(0) || d2 == -UWORD(1));
        if (d2)
            sub_ddmmss(d1, d0, q1, q0, p1, p0);
        FLINT_ASSERT(d1 < ctx->ffinfo->mod.n);
        NMOD_RED2(r, d1, d0, ctx->ffinfo->mod);
        if (r != UWORD(0))
        {
            success = 0;
            goto cleanup;
        }

        mpoly_get_monomial_ui(Aexp, A->exps + NA*i, A->bits, ctx->minfo);
        mpoly_get_monomial_ui(Bexp, B->exps + NB*i, B->bits, ctx->minfo);
        for (j = 0; j < nvars; j++)
        {
            if (Aexp[j] - Amin_exp[j] != Bexp[j] - Bmin_exp[j])
            {
                success = 0;
                goto cleanup;
            }
        }
    }

    success = 1;

    nmod_mpoly_init3(T, A->length, Gbits, ctx);

    /* put A's cofactor monomial in Bexp */
    for (j = 0; j < nvars; j++)
    {
        Bexp[j] = FLINT_MAX(Amin_exp[j], Bmin_exp[j]) - Bmin_exp[j];
    }

    a0inv = nmod_inv(A->coeffs[0], ctx->ffinfo->mod);

    NG = mpoly_words_per_exp(Gbits, ctx->minfo);

    for (i = 0; i < A->length; i++)
    {
        mpoly_get_monomial_ui(Aexp, A->exps + NA*i, A->bits, ctx->minfo);
        for (j = 0; j < nvars; j++)
        {
            Aexp[j] -= Bexp[j];
        }
        T->coeffs[i] = nmod_mul(A->coeffs[i], a0inv, ctx->ffinfo->mod);
        mpoly_set_monomial_ui(T->exps + NG*i, Aexp, Gbits, ctx->minfo);
    }
    _nmod_mpoly_set_length(T, A->length, ctx);

    nmod_mpoly_swap(T, G, ctx);
    nmod_mpoly_clear(T, ctx);

cleanup:

    TMP_END;

    return success;
}


/*
    General slow version that gets the job done for arbitrary A and B
*/
int _nmod_mpoly_gcd_monomial_cofactors(nmod_mpoly_t G,
                                const nmod_mpoly_t A, const nmod_mpoly_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    slong NA, NB;
    slong nvars = ctx->minfo->nvars;
    fmpz * Aexp, * Bexp, * Aexp0, * Bexp0, * minAexp, * minBexp;
    fmpz_t t1, t2;
    mp_limb_t a0inv;
    nmod_mpoly_t T;
    TMP_INIT;

    if ( A->length != B->length
         || A->length == 0
         || B->length == 0
       )
    {
        return 0;
    }

    TMP_START;

    Aexp = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    Bexp = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    Aexp0 = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    Bexp0 = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    minAexp = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    minBexp = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (j = 0; j < nvars; j++)
    {
        fmpz_init(Aexp + j);
        fmpz_init(Bexp + j);
        fmpz_init(Aexp0 + j);
        fmpz_init(Bexp0 + j);
        fmpz_init(minAexp + j);
        fmpz_init(minBexp + j);
    }

    fmpz_init(t1);
    fmpz_init(t2);

    NA = mpoly_words_per_exp(A->bits, ctx->minfo);
    NB = mpoly_words_per_exp(B->bits, ctx->minfo);

    mpoly_get_monomial_ffmpz(Aexp0, A->exps + NA*0, A->bits, ctx->minfo);
    mpoly_get_monomial_ffmpz(Bexp0, B->exps + NB*0, B->bits, ctx->minfo);
    _fmpz_vec_set(minAexp, Aexp0, nvars);
    _fmpz_vec_set(minBexp, Bexp0, nvars);

    for (i = 0; i < A->length; i++)
    {
        if (   nmod_mul(A->coeffs[0], B->coeffs[i], ctx->ffinfo->mod)
            != nmod_mul(B->coeffs[0], A->coeffs[i], ctx->ffinfo->mod))
        {
            success = 0;
            goto cleanup;
        }

        mpoly_get_monomial_ffmpz(Aexp, A->exps + NA*i, A->bits, ctx->minfo);
        mpoly_get_monomial_ffmpz(Bexp, B->exps + NB*i, B->bits, ctx->minfo);
        _fmpz_vec_min_inplace(minAexp, Aexp, nvars);
        _fmpz_vec_min_inplace(minBexp, Bexp, nvars);
        for (j = 0; j < nvars; j++)
        {
            fmpz_add(t1, Aexp0 + j, Bexp + j);
            fmpz_add(t2, Bexp0 + j, Aexp + j);
            if (!fmpz_equal(t1, t2))
            {
                success = 0;
                goto cleanup;
            }            
        }
    }

    success = 1;

    nmod_mpoly_init3(T, A->length, FLINT_MIN(A->bits, B->bits), ctx);

    /* put A's cofactor monomial in Bexp */
    _fmpz_vec_max(Bexp, minAexp, minBexp, nvars);
    _fmpz_vec_sub(Bexp, Bexp, minBexp, nvars);

    a0inv = nmod_inv(A->coeffs[0], ctx->ffinfo->mod);

    for (i = 0; i < A->length; i++)
    {
        mpoly_get_monomial_ffmpz(Aexp, A->exps + NA*i, A->bits, ctx->minfo);
        _fmpz_vec_sub(Aexp, Aexp, Bexp, nvars);
        _nmod_mpoly_push_exp_ffmpz(T, Aexp, ctx);
        FLINT_ASSERT(T->length == i + 1);
        T->coeffs[i] = nmod_mul(A->coeffs[i], a0inv, ctx->ffinfo->mod);
    }

    nmod_mpoly_swap(T, G, ctx);
    nmod_mpoly_clear(T, ctx);

cleanup:

    for (j = 0; j < nvars; j++)
    {
        fmpz_clear(Aexp + j);
        fmpz_clear(Bexp + j);
        fmpz_clear(Aexp0 + j);
        fmpz_clear(Bexp0 + j);
        fmpz_clear(minAexp + j);
        fmpz_clear(minBexp + j);
    }

    fmpz_clear(t1);
    fmpz_clear(t2);

    TMP_END;

    return success;
}
