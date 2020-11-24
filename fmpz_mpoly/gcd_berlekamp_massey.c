/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"
#include "nmod_mpoly.h"
#include "fmpz_mpoly.h"
#include "fmpz_mod_mpoly.h"
#include "n_poly.h"
#include "nmod_mpoly_factor.h"


#define LOW_HALF_MASK ((-UWORD(1)) >> (FLINT_BITS - FLINT_BITS/2))


void _mpoly_monomial_evals_nmod(
    mp_limb_t * E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const mp_limb_t * alpha,
    slong vstart,
    const mpoly_ctx_t mctx,
    nmod_t fctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(Abits, mctx);
    slong nvars = mctx->nvars;
    slong * LUToffset;
    ulong * LUTmask;
    mp_limb_t * LUTvalue;
    slong LUTlen;
    mp_limb_t xpoweval;
    ulong * inputexpmask;
    TMP_INIT;

    TMP_START;

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (mp_limb_t *) TMP_ALLOC(N*FLINT_BITS*sizeof(mp_limb_t));

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < Alen; i++)
    {
        for (j = 0; j < N; j++)
        {
            inputexpmask[j] |= (Aexps + N*i)[j];
        }
    }

    LUTlen = 0;
    for (j = nvars - 1; j >= vstart; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, Abits, mctx);

        xpoweval = alpha[j]; /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < Abits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            LUTvalue[LUTlen] = xpoweval;
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
                LUTlen++;
            xpoweval = nmod_mul(xpoweval, xpoweval, fctx);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    for (i = 0; i < Alen; i++)
    {
        xpoweval = 1;
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexps + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                xpoweval = nmod_mul(xpoweval, LUTvalue[j], fctx);
            }
        }
        E[i] = xpoweval;
    }

    TMP_END;
}


/*
    If I was formed from evaluations at
        alpha^alphashift, alpha^(alphashift + 1), ...
    construct the corresponding mpoly if possible with coeffs in (-p/2, p/2]
    The substitution degrees and degree bounds in Ictx are used.
*/
int nmod_mpoly_bma_get_fmpz_mpoly(
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx,
    ulong alphashift,
    nmod_berlekamp_massey_t I,
    const mpoly_bma_interpolate_ctx_t Ictx,
    nmod_t fpctx)
{
    slong i, j, t, N;
    int success;
    ulong new_exp, this_exp;
    slong * shifts, * offsets;
    mp_limb_t * values, * roots;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Alen;
    mp_limb_t T, S, V, V0, V1, V2, p0, p1, r;
    TMP_INIT;

    t = nmod_poly_degree(I->V1);
    FLINT_ASSERT(I->points->length >= t);

    if (t <= 0)
    {
        return 0;
    }

    TMP_START;

    /* use the rt member of I as temp space for roots - slightly dirty */
    nmod_poly_fit_length(I->rt, t);
    I->rt->length = t;
    roots = I->rt->coeffs;
    values = I->points->coeffs;

    success = nmod_poly_find_distinct_nonzero_roots(roots, I->V1);
    if (!success)
    {
        goto cleanup;
    }

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    fmpz_mpoly_fit_length(A, t, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;
    A->length = 0;

    shifts = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    offsets = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, ctx->minfo);
    }

    Alen = 0;
    for (i = 0; i < t; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        V0 = V1 = V2 = T = S = 0;
        r = roots[i];
        for (j = t; j > 0; j--)
        {
            T = nmod_add(nmod_mul(r, T, fpctx), I->V1->coeffs[j], fpctx);
            S = nmod_add(nmod_mul(r, S, fpctx), T, fpctx);
            umul_ppmm(p1, p0, values[j - 1], T);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, WORD(0), p1, p0);
        }
        /* roots[i] should be a root of master */
        FLINT_ASSERT(nmod_add(nmod_mul(r, T, fpctx), I->V1->coeffs[0], fpctx) == 0);
        NMOD_RED3(V, V2, V1, V0, fpctx);
        S = nmod_mul(S, nmod_pow_ui(r, alphashift, fpctx), fpctx);
        V0 = nmod_mul(V, nmod_inv(S, fpctx), fpctx);
        if (V0 == 0)
        {
            /* hmmm */
            continue;
        }

        if (fpctx.n - V0 < V0)
            fmpz_neg_ui(Acoeff + Alen, fpctx.n - V0);
        else
            fmpz_set_ui(Acoeff + Alen, V0);

        mpoly_monomial_zero(Aexp + N*Alen, N);
        new_exp = nmod_discrete_log_pohlig_hellman_run(Ictx->dlogenv_sp, roots[i]);
        for (j = ctx->minfo->nvars - 1; j >= 0; j--)
        {
            this_exp = new_exp % Ictx->subdegs[j];
            new_exp = new_exp / Ictx->subdegs[j];
            if (this_exp >= Ictx->degbounds[j])
            {
                success = 0;
                goto cleanup;
            }
            (Aexp + N*Alen)[offsets[j]] |= this_exp << shifts[j];
        }
        if (new_exp != 0)
        {
            success = 0;
            goto cleanup;
        }
        Alen++;
    }
    A->length = Alen;

    fmpz_mpoly_sort_terms(A, ctx);

    success = 1;

cleanup:

    TMP_END;
    return success;
}


/*
    fmpz_mod_mpoly "skeletons" - just the coefficients
    TODO get rid of this
*/

void fmpz_mpolyc_init(fmpz_mpolyc_t A)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

void fmpz_mpolyc_clear(fmpz_mpolyc_t A)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        fmpz_clear(A->coeffs + i);
    }
    if (A->coeffs)
        flint_free(A->coeffs);
}

void fmpz_mpolyc_fit_length(fmpz_mpolyc_t A, slong length)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->coeffs = (fmpz *) flint_malloc(new_alloc*sizeof(fmpz));
        }
        else
        {
            A->coeffs = (fmpz *) flint_realloc(A->coeffs, new_alloc*sizeof(fmpz));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_init(A->coeffs + i);
        }
        A->alloc = new_alloc;
    }
}

void fmpz_mpolycu_init(fmpz_mpolycu_t A)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

void fmpz_mpolycu_clear(fmpz_mpolycu_t A)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        fmpz_mpolyc_clear(A->coeffs + i);
    }
    if (A->coeffs)
        flint_free(A->coeffs);
}

void fmpz_mpolycu_fit_length(fmpz_mpolycu_t A, slong length)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->coeffs = (fmpz_mpolyc_struct *) flint_malloc(
                                         new_alloc*sizeof(fmpz_mpolyc_struct));
        }
        else
        {
            A->coeffs = (fmpz_mpolyc_struct *) flint_realloc(A->coeffs,
                                         new_alloc*sizeof(fmpz_mpolyc_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mpolyc_init(A->coeffs + i);
        }
        A->alloc = new_alloc;
    }
}


/*
    set out to the evaluation of variables after ksub at alpha^w

    out[0]   = alpha ^ (w * subdegs[n-1] * subdegs[n-2] * ... * * subdegs[1])
      ...
    out[n-3] = alpha ^ (w * subdegs[n-1] * subdegs[n-2])
    out[n-2] = alpha ^ (w * subdegs[n-1])
    out[n-1] = alpha ^ (w)

    secret: subdegs[0] is not used
*/

void nmod_mpoly_bma_interpolate_alpha_powers(
    mp_limb_t * out,
    ulong w,
    const mpoly_bma_interpolate_ctx_t Ictx,
    const fmpz_mpoly_ctx_t ctx,
    nmod_t fpctx)
{
    slong j = ctx->minfo->nvars - 1;
    out[j] = nmod_pow_ui(Ictx->dlogenv_sp->alpha, w, fpctx);
    for (; j > 0; j--)
    {
        out[j - 1] = nmod_pow_ui(out[j], Ictx->subdegs[j], fpctx);
    }
}

void fmpz_mod_mpoly_bma_interpolate_alpha_powers(
    fmpz * out,
    const fmpz_t w,
    const mpoly_bma_interpolate_ctx_t Ictx,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong j = ctx->minfo->nvars - 1;
    fmpz_mod_pow_fmpz(out + j, Ictx->dlogenv->alpha, w, fpctx);
    for (; j > 0; j--)
    {
        fmpz_mod_pow_ui(out + j - 1, out + j, Ictx->subdegs[j], fpctx);
    }
}


/*
    Set the coefficients of E to the evaluation of cooresponding monomials of A
    evaluation at x_i = alpha[i]
*/
void fmpz_mod_mpoly_set_skel(
    fmpz_mpolyc_t M,
    const fmpz_mod_mpoly_ctx_t ctx_mp,
    const fmpz_mpoly_t A,
    const fmpz * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    ulong * Aexp;
    fmpz * Mcoeff;
    slong * LUToffset;
    ulong * LUTmask;
    fmpz * LUTvalue;
    slong LUTlen;
    fmpz_t xpoweval;
    ulong * inputexpmask;
    TMP_INIT;

    TMP_START;

    fmpz_init(xpoweval);

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (fmpz *) TMP_ALLOC(N*FLINT_BITS*sizeof(fmpz));
    for (i = 0; i < N*FLINT_BITS; i++)
    {
        fmpz_init(LUTvalue + i);
    }

    fmpz_mpolyc_fit_length(M, A->length);
    M->length = A->length;

    Mcoeff = M->coeffs;
    Aexp = A->exps;

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < N; j++)
        {
            inputexpmask[j] |= (Aexp + N*i)[j];
        }
    }

    LUTlen = 0;
    for (j = nvars - 1; j >= 0; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, A->bits, ctx->minfo);

        fmpz_set(xpoweval, alpha + j); /* xpoweval = alpha[j]^(2^i) */
        for (i = 0; i < A->bits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            fmpz_set(LUTvalue + LUTlen, xpoweval);
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
            {
                LUTlen++;
            }
            fmpz_mod_mul(xpoweval, xpoweval, xpoweval, ctx_mp->ffinfo);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    for (i = 0; i < A->length; i++)
    {
        fmpz_one(xpoweval);
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                fmpz_mod_mul(xpoweval, xpoweval, LUTvalue + j, ctx_mp->ffinfo);
            }
        }
        fmpz_set(Mcoeff + i, xpoweval);
    }

    fmpz_clear(xpoweval);
    for (i = 0; i < N*FLINT_BITS; i++)
    {
        fmpz_clear(LUTvalue + i);
    }
    TMP_END;
}

void fmpz_mod_mpolyu_set_skel(
    fmpz_mpolycu_t M,
    const fmpz_mod_mpoly_ctx_t ctx_mp,
    const fmpz_mpolyu_t A,
    const fmpz * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mpolycu_fit_length(M, A->length);
    M->length = A->length;
    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_mpoly_set_skel(M->coeffs + i, ctx_mp, A->coeffs + i, alpha, ctx);
    }
}


/* M = S */
void fmpz_mod_mpoly_copy_skel(fmpz_mpolyc_t M, const fmpz_mpolyc_t S)
{
    fmpz_mpolyc_fit_length(M, S->length);
    M->length = S->length;
    _fmpz_vec_set(M->coeffs, S->coeffs, S->length);
}

void fmpz_mod_mpolyu_copy_skel(fmpz_mpolycu_t M, const fmpz_mpolycu_t S)
{
    slong i;
    fmpz_mpolycu_fit_length(M, S->length);
    M->length = S->length;
    for (i = 0; i < S->length; i++)
    {
        fmpz_mod_mpoly_copy_skel(M->coeffs + i, S->coeffs + i);
    }
}


/* Ared = A mod p */
void fmpz_mod_mpoly_red_skel(
    fmpz_mpolyc_t Ared,
    const fmpz_mpoly_t A,
    const fmpz_mod_ctx_t fpctx)
{
    fmpz_mpolyc_fit_length(Ared, A->length);
    Ared->length = A->length;
    _fmpz_vec_scalar_mod_fmpz(Ared->coeffs, A->coeffs, A->length,
                                                  fmpz_mod_ctx_modulus(fpctx));
}

void fmpz_mod_mpolyu_red_skel(
    fmpz_mpolycu_t Ared,
    const fmpz_mpolyu_t A,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    fmpz_mpolycu_fit_length(Ared, A->length);
    Ared->length = A->length;
    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_mpoly_red_skel(Ared->coeffs + i, A->coeffs + i, fpctx);
    }
}


/*
    return Ared.Avar and multiply Avar by Ainc
    the coefficients of Ared must be reduced mod fpctx->p
*/
void fmpz_mod_mpoly_use_skel_mul(
    fmpz_t eval,
    fmpz_mpolyc_t Ared,
    fmpz_mpolyc_t Avar,
    const fmpz_mpolyc_t Ainc,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    fmpz_zero(eval);
    FLINT_ASSERT(Avar->length == Ared->length);
    for (i = 0; i < Ared->length; i++)
    {
        fmpz_addmul(eval, Ared->coeffs + i, Avar->coeffs + i);
        fmpz_mod_mul(Avar->coeffs + i, Avar->coeffs + i, Ainc->coeffs + i, fpctx);
    }
    fmpz_mod(eval, eval, fmpz_mod_ctx_modulus(fpctx));
}


/*
    A is in ZZ[X,Y][x_0, ..., x_(n-1)]
    E is in Fp[X][Y]

    Set E to the evaluation of A. The evaluations of A's monomials are in M.
    M is then multiplied by S.
*/
void fmpz_mod_mpolyuu_use_skel_mul(
    fmpz_mod_mpolyn_t E,
    const fmpz_mpolyu_t A,
    const fmpz_mpolycu_t Ared,
    fmpz_mpolycu_t M,
    const fmpz_mpolycu_t S,
    const fmpz_mod_mpoly_ctx_t ctx_mp)
{
    slong xexp, yexp;
    slong i;
    fmpz_t eval;

    FLINT_ASSERT(E->bits == FLINT_BITS/2);
    FLINT_ASSERT(1 == mpoly_words_per_exp_sp(E->bits, ctx_mp->minfo));

    FLINT_ASSERT(A->length == Ared->length);
    FLINT_ASSERT(A->length == M->length);
    FLINT_ASSERT(A->length == S->length);

    fmpz_init(eval);

    E->length = 0;
    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_mpoly_use_skel_mul(eval, Ared->coeffs + i, M->coeffs + i,
                                                S->coeffs + i, ctx_mp->ffinfo);
        if (fmpz_is_zero(eval))
        {
            continue;
        }

        xexp = A->exps[i] >> (FLINT_BITS/2);
        yexp = A->exps[i] & LOW_HALF_MASK;

        if (E->length > 0 && (E->exps[E->length - 1] >> (FLINT_BITS/2)) == xexp)
        {
            fmpz_mod_poly_set_coeff_fmpz(E->coeffs + E->length - 1, yexp, eval, ctx_mp->ffinfo);
        }
        else
        {
            fmpz_mod_mpolyn_fit_length(E, E->length + 1, ctx_mp);
            fmpz_mod_poly_zero(E->coeffs + E->length, ctx_mp->ffinfo);
            fmpz_mod_poly_set_coeff_fmpz(E->coeffs + E->length, yexp, eval, ctx_mp->ffinfo);
            E->exps[E->length] = xexp << (FLINT_BITS/2);
            E->length++;
        }
    }

    fmpz_clear(eval);
}



int fmpz_mod_bma_get_fmpz_mpoly(
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_t alphashift,
    fmpz_mod_berlekamp_massey_t I,
    const mpoly_bma_interpolate_ctx_t Ictx,
    const fmpz_mod_ctx_t fpctx)
{
    slong i, j, t, N;
    int success;
    ulong this_exp;
    fmpz_t new_exp;
    slong * shifts, * offsets;
    fmpz * values, * roots;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Alen;
    fmpz_t T, S, V, temp, halfp;
    TMP_INIT;

    TMP_START;

    fmpz_init(halfp);
    fmpz_init(T);
    fmpz_init(S);
    fmpz_init(V);
    fmpz_init(temp);
    fmpz_init(new_exp);

    fmpz_tdiv_q_2exp(halfp, fmpz_mod_ctx_modulus(fpctx), 1);

    fmpz_mod_berlekamp_massey_reduce(I, fpctx);
    t = fmpz_mod_poly_degree(I->V1, fpctx);
    FLINT_ASSERT(t > 0);
    FLINT_ASSERT(I->points->length >= t);

    fmpz_mod_poly_fit_length(I->rt, t, fpctx);
    I->rt->length = t;
    success = fmpz_mod_poly_find_distinct_nonzero_roots(I->rt->coeffs, I->V1, fpctx);
    if (!success)
    {
        goto cleanup;
    }

    roots = I->rt->coeffs;
    values = I->points->coeffs;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    fmpz_mpoly_fit_length(A, t, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;
    A->length = 0;

    shifts = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    offsets = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, ctx->minfo);
    }

    Alen = 0;
    for (i = 0; i < t; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        fmpz_zero(V);
        fmpz_zero(T);
        fmpz_zero(S);
        for (j = t; j > 0; j--)
        {
            fmpz_mod_mul(temp, roots + i, T, fpctx);
            fmpz_mod_add(T, temp, I->V1->coeffs + j, fpctx);
            fmpz_mod_mul(temp, roots + i, S, fpctx);
            fmpz_mod_add(S, temp, T, fpctx);
            fmpz_mod_mul(temp, values + j - 1, T, fpctx);
            fmpz_mod_add(V, V, temp, fpctx);
        }
        /* roots[i] should be a root of master */
#if FLINT_WANT_ASSERT
        fmpz_mod_mul(temp, roots + i, T, fpctx);
        fmpz_mod_add(temp, temp, I->V1->coeffs + 0, fpctx);
        FLINT_ASSERT(fmpz_is_zero(temp));
#endif
        fmpz_mod_pow_fmpz(temp, roots + i, alphashift, fpctx);
        fmpz_mod_mul(S, S, temp, fpctx);
        fmpz_mod_inv(temp, S, fpctx);
        fmpz_mod_mul(Acoeff + Alen, V, temp, fpctx);
        if (fmpz_is_zero(Acoeff + Alen))
        {
            /* hmmm */
            continue;
        }

        if (fmpz_cmp(Acoeff + Alen, halfp) > 0)
        {
            fmpz_sub(Acoeff + Alen, Acoeff + Alen, fmpz_mod_ctx_modulus(fpctx));
        }


        mpoly_monomial_zero(Aexp + N*Alen, N);
        fmpz_mod_discrete_log_pohlig_hellman_run(new_exp, Ictx->dlogenv, roots + i);
        for (j = ctx->minfo->nvars - 1; j >= 0; j--)
        {
            this_exp = fmpz_fdiv_ui(new_exp, Ictx->subdegs[j]);
            fmpz_fdiv_q_ui(new_exp, new_exp, Ictx->subdegs[j]);
            if (this_exp >= Ictx->degbounds[j])
            {
                success = 0;
                goto cleanup;
            }
            (Aexp + N*Alen)[offsets[j]] |= this_exp << shifts[j];
        }
        if (!fmpz_is_zero(new_exp))
        {
            success = 0;
            goto cleanup;
        }
        Alen++;
    }
    A->length = Alen;

    fmpz_mpoly_sort_terms(A, ctx);

    success = 1;

cleanup:

    fmpz_clear(T);
    fmpz_clear(S);
    fmpz_clear(V);
    fmpz_clear(temp);
    fmpz_clear(halfp);

    TMP_END;
    return success;
}



void nmod_bma_mpoly_init(nmod_bma_mpoly_t A)
{
    A->length = 0;
    A->alloc = 0;
    A->exps = NULL;
    A->coeffs = NULL;
    A->pointcount = 0;
}

void nmod_bma_mpoly_reset_prime(
    nmod_bma_mpoly_t A,
    nmod_t fpctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        nmod_berlekamp_massey_set_prime(A->coeffs + i, fpctx.n);
    }
}


void nmod_bma_mpoly_clear(nmod_bma_mpoly_t A)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        nmod_berlekamp_massey_clear(A->coeffs + i);
    }

    if (A->exps)
        flint_free(A->exps);
    if (A->coeffs)
        flint_free(A->coeffs);
}

void nmod_bma_mpoly_print(const nmod_bma_mpoly_t A)
{
    slong i;
    flint_printf("0");
    for (i = 0; i < A->length; i++)
    {
        flint_printf(" + [");
        nmod_berlekamp_massey_print(A->coeffs + i);
        flint_printf("]*X^%wd*Y^%wd",
                 A->exps[i] >> (FLINT_BITS/2),
                 A->exps[i] & ((-UWORD(1)) >> (FLINT_BITS/2)));
    }
}


void nmod_bma_mpoly_fit_length(
    nmod_bma_mpoly_t A,
    slong length,
    nmod_t fpctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (nmod_berlekamp_massey_struct *) flint_malloc(
                               new_alloc*sizeof(nmod_berlekamp_massey_struct));
        }
        else
        {
            A->exps = (ulong *) flint_realloc(A->exps, new_alloc*sizeof(ulong));
            A->coeffs = (nmod_berlekamp_massey_struct *) flint_realloc(A->coeffs,
                               new_alloc*sizeof(nmod_berlekamp_massey_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            nmod_berlekamp_massey_init(A->coeffs + i, fpctx.n);
        }
        A->alloc = new_alloc;
    }
}

void nmod_bma_mpoly_zero(nmod_bma_mpoly_t L)
{
    L->length = 0;
    L->pointcount = 0;
}


int nmod_bma_mpoly_reduce(nmod_bma_mpoly_t L)
{
    slong i;
    int changed;

    changed = 0;

    for (i = 0; i < L->length; i++)
    {
        FLINT_ASSERT(L->pointcount == nmod_berlekamp_massey_point_count(L->coeffs + i));
        changed |= nmod_berlekamp_massey_reduce(L->coeffs + i);
    }

    return changed;
}


void nmod_bma_mpoly_add_point(
    nmod_bma_mpoly_t L,
    const n_bpoly_t A,
    const nmod_mpoly_ctx_t ctx_sp)
{
    slong j;
    slong Alen = A->length;
    n_poly_struct * Acoeffs = A->coeffs;
    slong Li, Ai, ai;
    nmod_berlekamp_massey_struct * Lcoeff;
    slong Llen;
    ulong * Lexp;

    if (L->length == 0)
    {
        slong tot = 0;
        for (Ai = 0; Ai < Alen; Ai++)
        for (ai = n_poly_degree(Acoeffs + Ai); ai >= 0; ai--)
            tot += (0 != Acoeffs[Ai].coeffs[ai]);
        nmod_bma_mpoly_fit_length(L, tot, ctx_sp->mod);
    }

    Lcoeff = L->coeffs;
    Llen = L->length;
    Lexp = L->exps;

    Li = 0;
    Ai = A->length - 1;
    ai = 0;
    if (Ai >= 0)
        ai = n_poly_degree(Acoeffs + Ai);

    while (Li < Llen || Ai >= 0)
    {
        if (Li < Llen && Ai >= 0 && Lexp[Li] == pack_exp2(Ai, ai))
        {
            /* L term present, A term present */
add_same_exp:
            nmod_berlekamp_massey_add_point(Lcoeff + Li, Acoeffs[Ai].coeffs[ai]);
            Li++;

            do {
                ai--;
            } while (ai >= 0 && Acoeffs[Ai].coeffs[ai] == 0);
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_poly_degree(Acoeffs + Ai);
            }
        }
        else if (Li < Llen && (Ai < 0 || Lexp[Li] > pack_exp2(Ai, ai)))
        {
            /* L term present, A term missing */
            nmod_berlekamp_massey_add_zeros(Lcoeff + Li, 1);
            Li++;
        }
        else
        {
            /* L term missing, A term present */
            FLINT_ASSERT(Ai >= 0 && (Li >= Llen || Lexp[Li] < pack_exp2(Ai, ai)));
            {
                ulong texp;
                nmod_berlekamp_massey_struct tcoeff;

                nmod_bma_mpoly_fit_length(L, Llen + 1, ctx_sp->mod);
                Lcoeff = L->coeffs;
                Lexp = L->exps;

                texp = Lexp[Llen];
                tcoeff = Lcoeff[Llen];
                for (j = Llen - 1; j >= Li; j--)
                {
                    Lexp[j + 1] = Lexp[j];
                    Lcoeff[j + 1] = Lcoeff[j];
                }
                Lexp[Li] = texp;
                Lcoeff[Li] = tcoeff;
            }

            nmod_berlekamp_massey_start_over(Lcoeff + Li);
            nmod_berlekamp_massey_add_zeros(Lcoeff + Li, L->pointcount);
            Lexp[Li] = pack_exp2(Ai, ai);
            Llen++;
            L->length = Llen;

            goto add_same_exp;
        }
    }

    L->pointcount++;
}



int nmod_zip_mpolyuu_add_point(
    nmod_zip_mpolyu_t L,
    const n_bpoly_t A)
{
    n_poly_struct * Acoeffs = A->coeffs;
    slong Li, Ai, ai;
    nmod_zip_struct * Lcoeff;
    slong Llen;
    ulong * Lexp;
    slong pointcount = L->pointcount;

    Lcoeff = L->coeffs;
    Llen = L->length;
    Lexp = L->exps;

    Li = 0;
    Ai = A->length - 1;
    ai = 0;
    if (Ai >= 0)
        ai = n_poly_degree(Acoeffs + Ai);

    for (Li = 0; Li < Llen; Li++)
    {
        nmod_zip_struct * Lc = Lcoeff + Li;

        if (Ai >= 0 && Lexp[Li] == pack_exp2(Ai, ai))
        {
            /* L present A present */
            FLINT_ASSERT(pointcount <= Lc->ealloc);
            Lc->evals[pointcount] = Acoeffs[Ai].coeffs[ai];
            do {
                ai--;
            } while (ai >= 0 && Acoeffs[Ai].coeffs[ai] == 0);
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_poly_degree(Acoeffs + Ai);
            }
        }
        else if (Ai < 0 || Lexp[Li] > pack_exp2(Ai, ai))
        {
            /* L present A missing */
            FLINT_ASSERT(pointcount <= Lc->ealloc);
            Lc->evals[pointcount] = 0;
        }
        else
        {
            /* L missing A present */
            return 0;
        }
    }

    L->pointcount = pointcount + 1;
    return 1;
}




int nmod_bma_mpoly_get_fmpz_mpolyu(
    fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx,
    ulong alphashift,
    const nmod_bma_mpoly_t L,
    const mpoly_bma_interpolate_ctx_t Ictx,
    nmod_t fpctx)
{
    int success;
    slong i;

    fmpz_mpolyu_fit_length(A, L->length, ctx);
    A->length = 0;
    for (i = 0; i < L->length; i++)
    {
        A->exps[A->length] = L->exps[i];
        success = nmod_mpoly_bma_get_fmpz_mpoly(A->coeffs + A->length, ctx,
                                       alphashift, L->coeffs + i, Ictx, fpctx);
        if (!success)
        {
            return 0;
        }
        A->length += !fmpz_mpoly_is_zero(A->coeffs + A->length, ctx);
    }
    return 1;
}



void fmpz_mod_bma_mpoly_init(fmpz_mod_bma_mpoly_t A)
{
    A->length = 0;
    A->alloc = 0;
    A->exps = NULL;
    A->coeffs = NULL;
    A->pointcount = 0;
}


void fmpz_mod_bma_mpoly_clear(
    fmpz_mod_bma_mpoly_t A,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
    {
        fmpz_mod_berlekamp_massey_clear(A->coeffs + i, fpctx);
    }

    if (A->exps)
        flint_free(A->exps);
    if (A->coeffs)
        flint_free(A->coeffs);
}

#if 0
void fmpz_mod_bma_mpoly_print(
    fmpz_mod_bma_mpoly_t A,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    flint_printf("0");
    for (i = 0; i < A->length; i++)
    {
        flint_printf(" + [");
        fmpz_mod_berlekamp_massey_print(A->coeffs + i, fpctx);
        flint_printf("]*X^%wd*Y^%wd",
                 A->exps[i] >> (FLINT_BITS/2),
                 A->exps[i] & ((-UWORD(1)) >> (FLINT_BITS - FLINT_BITS/2)));
    }
}
#endif


void fmpz_mod_bma_mpoly_fit_length(
    fmpz_mod_bma_mpoly_t A,
    slong length,
    const fmpz_mod_ctx_t fpctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mod_berlekamp_massey_struct *) flint_malloc(
                           new_alloc*sizeof(fmpz_mod_berlekamp_massey_struct));
        }
        else
        {
            A->exps = (ulong *) flint_realloc(A->exps, new_alloc*sizeof(ulong));
            A->coeffs = (fmpz_mod_berlekamp_massey_struct *) flint_realloc(A->coeffs,
                           new_alloc*sizeof(fmpz_mod_berlekamp_massey_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mod_berlekamp_massey_init(A->coeffs + i, fpctx);
        }
        A->alloc = new_alloc;
    }
}

void fmpz_mod_bma_mpoly_zero(fmpz_mod_bma_mpoly_t L)
{
    L->length = 0;
    L->pointcount = 0;
}

int fmpz_mod_bma_mpoly_reduce(fmpz_mod_bma_mpoly_t L, const fmpz_mod_ctx_t fpctx)
{
    slong i;
    int changed;

    changed = 0;

    for (i = 0; i < L->length; i++)
    {
        FLINT_ASSERT(L->pointcount == fmpz_mod_berlekamp_massey_point_count(L->coeffs + i));
        changed |= fmpz_mod_berlekamp_massey_reduce(L->coeffs + i, fpctx);
    }

    return changed;
}

void fmpz_mod_bma_mpoly_add_point(
    fmpz_mod_bma_mpoly_t L,
    const fmpz_mod_mpolyn_t A,
    const fmpz_mod_mpoly_ctx_t ctx_mp)
{
    slong j;
    slong Alen = A->length;
    fmpz_mod_poly_struct * Acoeff = A->coeffs;
    slong Li, Ai, ai;
    ulong Aexp;
    fmpz_mod_berlekamp_massey_struct * Lcoeff;
    slong Llen;
    ulong * Lexp;

    if (L->length == 0)
    {
        slong tot = 0;
        for (Ai = 0; Ai < Alen; Ai++)
        {
            for (ai = (Acoeff + Ai)->length - 1; ai >= 0; ai--)
            {
                tot += !fmpz_is_zero((Acoeff + Ai)->coeffs + ai);
            }
        }
        fmpz_mod_bma_mpoly_fit_length(L, tot, ctx_mp->ffinfo);
    }

    Lcoeff = L->coeffs;
    Llen = L->length;
    Lexp = L->exps;

    Li = 0;
    Ai = ai = 0;
    Aexp = 0;
    if (Ai < Alen)
    {
        ai = fmpz_mod_poly_degree(A->coeffs + Ai, ctx_mp->ffinfo);
        FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/2);
        FLINT_ASSERT(0 == (A->exps[Ai] & LOW_HALF_MASK));
        Aexp = A->exps[Ai] + ai;
    }

    while (Li < Llen || Ai < Alen)
    {
        if (Li < Llen && Ai < Alen && Lexp[Li] == Aexp)
        {
            /* L term present, A term present */
add_same_exp:
            fmpz_mod_berlekamp_massey_add_point(Lcoeff + Li,
                                   (Acoeff + Ai)->coeffs + ai, ctx_mp->ffinfo);
            Li++;
            do {
                ai--;
            } while (ai >= 0 && fmpz_is_zero((Acoeff + Ai)->coeffs + ai));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                {
                    ai = fmpz_mod_poly_degree(A->coeffs + Ai, ctx_mp->ffinfo);
                    FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/2);
                    FLINT_ASSERT(0 == (A->exps[Ai] & LOW_HALF_MASK));
                    Aexp = A->exps[Ai] + ai;
                }
            }
            else
            {
                FLINT_ASSERT(Ai < Alen);
                FLINT_ASSERT(FLINT_BIT_COUNT(ai) < FLINT_BITS/2);
                FLINT_ASSERT(0 == (A->exps[Ai] & LOW_HALF_MASK));
                Aexp = A->exps[Ai] + ai;
            }
        }
        else if (Li < Llen && (Ai >= Alen || Lexp[Li] > Aexp))
        {
            /* L term present, A term missing */
            fmpz_mod_berlekamp_massey_add_zeros(Lcoeff + Li, 1, ctx_mp->ffinfo);
            Li++;
        }
        else
        {
            /* L term missing, A term present */
            FLINT_ASSERT(Ai < Alen && (Li >= Llen || Lexp[Li] < Aexp));
            {
                ulong texp;
                fmpz_mod_berlekamp_massey_struct tcoeff;

                fmpz_mod_bma_mpoly_fit_length(L, Llen + 1, ctx_mp->ffinfo);
                Lcoeff = L->coeffs;
                Lexp = L->exps;

                texp = Lexp[Llen];
                tcoeff = Lcoeff[Llen];
                for (j = Llen - 1; j >= Li; j--)
                {
                    Lexp[j + 1] = Lexp[j];
                    Lcoeff[j + 1] = Lcoeff[j];
                }
                Lexp[Li] = texp;
                Lcoeff[Li] = tcoeff;
            }

            fmpz_mod_berlekamp_massey_start_over(Lcoeff + Li, ctx_mp->ffinfo);
            fmpz_mod_berlekamp_massey_add_zeros(Lcoeff + Li, L->pointcount, ctx_mp->ffinfo);
            Lexp[Li] = Aexp;
            Llen++;
            L->length = Llen;
            goto add_same_exp;
        }
    }

    L->pointcount++;
}

int fmpz_mod_bma_mpoly_get_fmpz_mpolyu(
    fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_t alphashift,
    const fmpz_mod_bma_mpoly_t L,
    const mpoly_bma_interpolate_ctx_t Ictx,
    const fmpz_mod_ctx_t fpctx)
{
    int success;
    slong i;

    fmpz_mpolyu_fit_length(A, L->length, ctx);
    A->length = 0;
    for (i = 0; i < L->length; i++)
    {
        A->exps[A->length] = L->exps[i];
        success = fmpz_mod_bma_get_fmpz_mpoly(A->coeffs + A->length, ctx,
                                       alphashift, L->coeffs + i, Ictx, fpctx);
        if (!success)
        {
            return 0;
        }
        A->length += !fmpz_mpoly_is_zero(A->coeffs + A->length, ctx);
    }
    return 1;
}



ulong fmpz_mod_mpolyn_bidegree(const fmpz_mod_mpolyn_t A)
{
    ulong degx, degy;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->bits == FLINT_BITS/2);

    degx = (A->exps + 1*0)[0] >> (FLINT_BITS/2);
    degy = A->coeffs[0].length - 1;

    FLINT_ASSERT(FLINT_BIT_COUNT(degx) < FLINT_BITS/2);
    FLINT_ASSERT(FLINT_BIT_COUNT(degy) < FLINT_BITS/2);

    return (degx << (FLINT_BITS/2)) + degy;
}

ulong nmod_mpolyn_bidegree(const nmod_mpolyn_t A)
{
    ulong degx, degy;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->bits == FLINT_BITS/2);

    degx = A->exps[0] >> (FLINT_BITS/2);
    degy = nmod_poly_degree(A->coeffs + 0);

    FLINT_ASSERT(FLINT_BIT_COUNT(degx) < FLINT_BITS/2);
    FLINT_ASSERT(FLINT_BIT_COUNT(degy) < FLINT_BITS/2);

    return (degx << (FLINT_BITS/2)) + degy;
}


void fmpz_mpoly_eval_fmpz_mod(
    fmpz_t eval,
    const fmpz_mod_ctx_t fpctx,
    const fmpz_mpoly_t A,
    const fmpz * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    ulong * Aexp;
    slong * LUToffset;
    ulong * LUTmask;
    fmpz * LUTvalue;
    slong LUTlen;
    fmpz_t xpoweval;
    ulong * inputexpmask;
    TMP_INIT;

    TMP_START;

    fmpz_init(xpoweval);

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (fmpz *) TMP_ALLOC(N*FLINT_BITS*sizeof(fmpz));
    for (i = 0; i < N*FLINT_BITS; i++)
    {
        fmpz_init(LUTvalue + i);
    }

    Aexp = A->exps;

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < N; j++)
        {
            inputexpmask[j] |= (Aexp + N*i)[j];
        }
    }

    LUTlen = 0;
    for (j = nvars - 1; j >= 0; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, A->bits, ctx->minfo);

        fmpz_set(xpoweval, alpha + j); /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < A->bits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            fmpz_set(LUTvalue + LUTlen, xpoweval);
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
            {
                LUTlen++;
            }
            fmpz_mod_mul(xpoweval, xpoweval, xpoweval, fpctx);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    fmpz_zero(eval);
    for (i = 0; i < A->length; i++)
    {
        fmpz_mod(xpoweval, A->coeffs + i, fmpz_mod_ctx_modulus(fpctx));
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                fmpz_mod_mul(xpoweval, xpoweval, LUTvalue + j, fpctx);
            }
        }
        fmpz_mod_add(eval, eval, xpoweval, fpctx);
    }

    fmpz_clear(xpoweval);
    for (i = 0; i < N*FLINT_BITS; i++)
    {
        fmpz_clear(LUTvalue + i);
    }
    TMP_END;
}


void fmpz_mpolyuu_eval_nmod(
    n_bpoly_t E,
    const nmod_mpoly_ctx_t ctx_sp,
    const fmpz_mpolyu_t A,
    const mp_limb_t * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong xexp, yexp;
    mp_limb_t eval;

    E->length = 0;
    for (i = 0; i < A->length; i++)
    {
        eval = fmpz_mpoly_eval_nmod(ctx_sp->mod, A->coeffs + i, alpha, ctx);
        if (eval == 0)
            continue;

        xexp = A->exps[i] >> (FLINT_BITS/2);
        yexp = A->exps[i] & LOW_HALF_MASK;

        n_bpoly_set_coeff(E, xexp, yexp, eval);
    }
}


void fmpz_mpolyuu_eval_fmpz_mod(
    fmpz_mod_mpolyn_t E,
    const fmpz_mod_mpoly_ctx_t ctx_mp,
    const fmpz_mpolyu_t A,
    const fmpz * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong xexp, yexp;
    fmpz_t eval;

    FLINT_ASSERT(E->bits == FLINT_BITS/2);
    FLINT_ASSERT(1 == mpoly_words_per_exp_sp(E->bits, ctx_mp->minfo));

    fmpz_init(eval);

    E->length = 0;
    for (i = 0; i < A->length; i++)
    {
        fmpz_mpoly_eval_fmpz_mod(eval, ctx_mp->ffinfo, A->coeffs + i, alpha, ctx);
        if (fmpz_is_zero(eval))
        {
            continue;
        }

        xexp = A->exps[i] >> (FLINT_BITS/2);
        yexp = A->exps[i] & LOW_HALF_MASK;

        if (E->length > 0 && (E->exps[E->length - 1] >> (FLINT_BITS/2)) == xexp)
        {
            fmpz_mod_poly_set_coeff_fmpz(E->coeffs + E->length - 1, yexp, eval, ctx_mp->ffinfo);
        }
        else
        {
            fmpz_mod_mpolyn_fit_length(E, E->length + 1, ctx_mp);
            fmpz_mod_poly_zero(E->coeffs + E->length, ctx_mp->ffinfo);
            fmpz_mod_poly_set_coeff_fmpz(E->coeffs + E->length, yexp, eval, ctx_mp->ffinfo);
            E->exps[E->length] = xexp << (FLINT_BITS/2);
            E->length++;
        }
    }

    fmpz_clear(eval);
}


mp_limb_t fmpz_mpoly_eval_nmod(
    nmod_t fpctx,
    const fmpz_mpoly_t A,
    const mp_limb_t * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    mp_limb_t eval;
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    ulong * Aexp;
    slong * LUToffset;
    ulong * LUTmask;
    mp_limb_t * LUTvalue;
    slong LUTlen;
    mp_limb_t xpoweval;
    ulong * inputexpmask;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);

    TMP_START;

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (mp_limb_t *) TMP_ALLOC(N*FLINT_BITS*sizeof(mp_limb_t));

    Aexp = A->exps;

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < N; j++)
        {
            inputexpmask[j] |= (Aexp + N*i)[j];
        }
    }

    LUTlen = 0;
    for (j = nvars - 1; j >= 0; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, A->bits, ctx->minfo);

        xpoweval = alpha[j]; /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < A->bits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            LUTvalue[LUTlen] = xpoweval;
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
            {
                LUTlen++;
            }
            xpoweval = nmod_mul(xpoweval, xpoweval, fpctx);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    eval = 0;
    for (i = 0; i < A->length; i++)
    {
        xpoweval = fmpz_fdiv_ui(A->coeffs + i, fpctx.n);
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                xpoweval = nmod_mul(xpoweval, LUTvalue[j], fpctx);
            }
        }
        eval = nmod_add(eval, xpoweval, fpctx);
    }

    TMP_END;

    return eval;
}



void nmod_zip_init(nmod_zip_t Z)
{
    Z->mlength = 0;
    Z->malloc = 0;
    Z->coeffs = NULL;
    Z->monomials = NULL;
    Z->ealloc = 0;
    Z->evals = NULL;
}

void nmod_zip_clear(nmod_zip_t Z)
{
    if (Z->coeffs)
        flint_free(Z->coeffs);
    if (Z->monomials)
        flint_free(Z->monomials);
    if (Z->evals)
        flint_free(Z->evals);
}

void nmod_zip_set_lengths(nmod_zip_t A, slong mlength, slong elength)
{
    slong old_malloc = A->malloc;
    slong new_malloc = FLINT_MAX(mlength, A->malloc + A->malloc/2);
    slong old_ealloc = A->ealloc;
    slong new_ealloc = FLINT_MAX(elength, A->ealloc + A->ealloc/2);

    FLINT_ASSERT(mlength > 0);
    FLINT_ASSERT(elength > 0);

    if (mlength > old_malloc)
    {
        if (old_malloc == 0)
        {
            A->coeffs    = (mp_limb_t *) flint_malloc(
                                                 new_malloc*sizeof(mp_limb_t));
            A->monomials = (mp_limb_t *) flint_malloc(
                                                 new_malloc*sizeof(mp_limb_t));
        }
        else
        {
            A->coeffs    = (mp_limb_t *) flint_realloc(A->coeffs,
                                                 new_malloc*sizeof(mp_limb_t));
            A->monomials = (mp_limb_t *) flint_realloc(A->monomials,
                                                 new_malloc*sizeof(mp_limb_t));
        }
        A->malloc = new_malloc;
    }

    A->mlength = mlength;

    if (elength > old_ealloc)
    {
        if (old_ealloc == 0)
        {
            A->evals = (mp_limb_t *) flint_malloc(new_ealloc*sizeof(mp_limb_t));
        }
        else
        {
            A->evals = (mp_limb_t *) flint_realloc(A->evals,
                                                 new_ealloc*sizeof(mp_limb_t));
        }
        A->ealloc = new_ealloc;
    }
}

void nmod_zip_print(const nmod_zip_t Z, slong elength)
{
    slong i;

    printf("m ");
    for (i = 0; i < Z->mlength; i++)
    {
        flint_printf("(%wu %wu) ", Z->coeffs[i], Z->monomials[i]);
    }
    printf("e ");
    for (i = 0; i < elength; i++)
    {
        flint_printf("%wu ", Z->evals[i]);
    }
}



void nmod_zip_mpolyu_init(nmod_zip_mpolyu_t Z)
{
    Z->alloc = 0;
    Z->length = 0;
    Z->exps = NULL;
    Z->coeffs = NULL;
    Z->pointcount = 0;
}

void nmod_zip_mpolyu_clear(nmod_zip_mpolyu_t Z)
{
    slong i;

    for (i = 0; i < Z->alloc; i++)
    {
        nmod_zip_clear(Z->coeffs + i);
    }

    if (Z->exps)
        flint_free(Z->exps);
    if (Z->coeffs)
        flint_free(Z->coeffs);
}

void nmod_zip_mpolyu_fit_length(
    nmod_zip_mpolyu_t A,
    slong length)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, A->alloc + A->alloc/2);

    if (length > old_alloc)
    {
        if (old_alloc == 0)
        {
            A->exps = (ulong *) flint_malloc(new_alloc*sizeof(ulong));
            A->coeffs = (nmod_zip_struct *) flint_malloc(
                                            new_alloc*sizeof(nmod_zip_struct));
        }
        else
        {
            A->exps = (ulong *) flint_realloc(A->exps, new_alloc*sizeof(ulong));
            A->coeffs = (nmod_zip_struct *) flint_realloc(A->coeffs,
                                            new_alloc*sizeof(nmod_zip_struct));
        }

        for (i = old_alloc; i < new_alloc; i++)
        {
            nmod_zip_init(A->coeffs + i);
        }
        A->alloc = new_alloc;
    }
}

void nmod_zip_mpolyu_fit_poly(
    nmod_zip_mpolyu_t Z,
    fmpz_mpolyu_t H,
    slong eval_length)
{
    slong i;
    FLINT_ASSERT(H->length > 0);

    nmod_zip_mpolyu_fit_length(Z, H->length);

    for (i = 0; i < H->length; i++)
    {
        Z->exps[i] = H->exps[i];
        nmod_zip_set_lengths(Z->coeffs + i, (H->coeffs + i)->length, eval_length);
    }

    Z->length = H->length;
    Z->pointcount = 0;
}


void nmod_zip_mpolyu_set_skel(
    nmod_zip_mpolyu_t Z,
    const nmod_mpoly_ctx_t ctx_sp,
    const fmpz_mpolyu_t A,
    const mp_limb_t * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    n_poly_t T;

    n_poly_init(T);

    FLINT_ASSERT(Z->length == A->length);
    for (i = 0; i < Z->length; i++)
    {
        nmod_zip_struct * Zc = Z->coeffs + i;

        n_poly_fit_length(T, A->coeffs[i].length);
        T->length = A->coeffs[i].length;
        _mpoly_monomial_evals_nmod(T->coeffs, A->coeffs[i].exps,
                        A->coeffs[i].bits, A->coeffs[i].length, alpha, 0,
                                                      ctx->minfo, ctx_sp->mod);
        Z->exps[i] = A->exps[i];
        FLINT_ASSERT(Zc->mlength == T->length);
        for (j = 0; j < Zc->mlength; j++)
        {
            Zc->coeffs[j] = 0;
            Zc->monomials[j] = T->coeffs[j];
        }
    }
    Z->pointcount = 0;

    n_poly_clear(T);
}

void nmod_zip_mpolyuu_print(const nmod_zip_mpolyu_t A)
{
    slong i;
    flint_printf("0");
    for (i = 0; i < A->length; i++)
    {
        flint_printf(" + [");
        nmod_zip_print(A->coeffs + i, A->pointcount);
        flint_printf("]*X^%wd*Y^%wd",
                 A->exps[i] >> (FLINT_BITS/2),
                 A->exps[i] & ((-UWORD(1)) >> (FLINT_BITS/2)));
    }
}


void nmod_mpolyuu_eval_step2(
    n_bpoly_t E,
    n_bpoly_t Acur,
    const n_polyun_t Ainc,
    const nmod_mpoly_ctx_t ctx_sp)
{
    slong xexp, yexp;
    slong i;
    slong Alen = Acur->length;
    n_polyun_term_struct * Aterms = Ainc->terms;
    mp_limb_t eval;

    FLINT_ASSERT(Alen == Ainc->length);

    E->length = 0;
    for (i = 0; i < Alen; i++)
    {
        eval = n_poly_mod_eval_step2(Acur->coeffs + i, Aterms[i].coeff,
                                                                  ctx_sp->mod);
        if (eval == 0)
        {
            continue;
        }

        xexp = Aterms[i].exp >> (FLINT_BITS/2);
        yexp = Aterms[i].exp & LOW_HALF_MASK;

        n_bpoly_set_coeff(E, xexp, yexp, eval);
    }

    FLINT_ASSERT(n_bpoly_mod_is_canonical(E, ctx_sp->mod));
}



void n_bpoly_set_mpolyn2(
    n_bpoly_t A,
    const nmod_mpolyn_t B,
    const nmod_mpoly_ctx_t ctx_sp)
{
    slong i, x;

    FLINT_ASSERT(B->bits == FLINT_BITS/2);
    FLINT_ASSERT(1 == mpoly_words_per_exp_sp(B->bits, ctx_sp->minfo));

    A->length = 0;

    for (i = 0; i < B->length; i++)
    {
        x = B->exps[i] >> (FLINT_BITS/2);
        n_bpoly_fit_length(A, x + 1);
        while (A->length <= x)
        {
            n_poly_zero(A->coeffs + A->length);
            A->length++;
        }

        n_poly_set_nmod_poly(A->coeffs + x, B->coeffs + i);
    }

    FLINT_ASSERT(n_bpoly_mod_is_canonical(A, ctx_sp->mod));
}


nmod_zip_find_coeffs_ret_t nmod_zip_find_coeffs(
    nmod_zip_t Z,
    nmod_poly_t master,
    slong pointcount,
    nmod_t fpctx)
{
    slong i, j;
    mp_limb_t V, V0, V1, V2, T, S, r, p0, p1;

    FLINT_ASSERT(pointcount >= Z->mlength);

    nmod_poly_product_roots_nmod_vec(master, Z->monomials, Z->mlength);

    for (i = 0; i < Z->mlength; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        V0 = V1 = V2 = T = S = 0;
        r = Z->monomials[i];
        for (j = Z->mlength; j > 0; j--)
        {
            T = nmod_add(nmod_mul(r, T, fpctx), master->coeffs[j], fpctx);
            S = nmod_add(nmod_mul(r, S, fpctx), T, fpctx);
            umul_ppmm(p1, p0, Z->evals[j - 1], T);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, WORD(0), p1, p0);
        }
        /* roots[i] should be a root of master */
        FLINT_ASSERT(nmod_add(nmod_mul(r, T, fpctx), master->coeffs[0], fpctx) == 0);
        NMOD_RED3(V, V2, V1, V0, fpctx);
        S = nmod_mul(S, r, fpctx); /* shift is one */
        if (S == 0)
        {
            return nmod_zip_find_coeffs_non_invertible;
        }
        Z->coeffs[i] = nmod_mul(V, nmod_inv(S, fpctx), fpctx);
    }

    /* use the coefficients of master as temp work space */
    for (i = 0; i < Z->mlength; i++)
    {
        master->coeffs[i] = nmod_pow_ui(Z->monomials[i], Z->mlength, fpctx);
    }

    /* check that the remaining points match */
    for (i = Z->mlength; i < pointcount; i++)
    {
        V0 = V1 = V2 = S = 0;
        for (j = 0; j < Z->mlength; j++)
        {
            master->coeffs[j] = nmod_mul(master->coeffs[j], Z->monomials[j], fpctx);
            umul_ppmm(p1, p0, Z->coeffs[j], master->coeffs[j]);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, WORD(0), p1, p0);
        }
        NMOD_RED3(V, V2, V1, V0, fpctx);
        if (V != Z->evals[i])
        {
            return nmod_zip_find_coeffs_no_match;
        }
    }

    return nmod_zip_find_coeffs_good;
}

nmod_zip_find_coeffs_ret_t nmod_mpolyu_zip_find_coeffs(
    nmod_zip_mpolyu_t Z,
    const nmod_mpoly_ctx_t ctx_sp)
{
    slong i;
    nmod_zip_find_coeffs_ret_t ret;
    nmod_poly_t T;

    nmod_poly_init_mod(T, ctx_sp->mod);

    for (i = 0; i < Z->length; i++)
    {
        ret = nmod_zip_find_coeffs(Z->coeffs + i, T, Z->pointcount, ctx_sp->mod);
        if (ret != nmod_zip_find_coeffs_good)
        {
            goto cleanup;
        }
    }

    ret = nmod_zip_find_coeffs_good;

cleanup:

    nmod_poly_clear(T);

    return ret;
}


int fmpz_mpolyu_addinterp_zip(
    fmpz_mpolyu_t H,
    const fmpz_t Hmodulus,
    const nmod_zip_mpolyu_t Z,
    nmod_t fpctx)
{
    int changed = 0;
    slong i, j;
    fmpz_t t;

    fmpz_init(t);

    FLINT_ASSERT(H->length == Z->length);
    for (i = 0; i < H->length; i++)
    {
        fmpz_mpoly_struct * Hc = H->coeffs + i;
        nmod_zip_struct * Zc = Z->coeffs + i;

        FLINT_ASSERT(Hc->length == Zc->mlength);
        for (j = 0; j < Hc->length; j++)
        {
            fmpz_CRT_ui(t, Hc->coeffs + j, Hmodulus, Zc->coeffs[j], fpctx.n, 1);
            changed |= !fmpz_equal(t, Hc->coeffs + j);
            fmpz_swap(t, Hc->coeffs + j);
        }
    }

    fmpz_clear(t);
    return changed;
}


/*
    A is in ZZ[x_0, ..., x_(n-1)]

    After the substitutions
        x_0     = x ^ (sub[1] * sub[2] * ... * sub[n-1])

        x_(n-2) = x ^ (sub[n-1])
        x_(n-1) = x ^ (1)
    a univariate in ZZ[x] remains. Return the content of this poly.
*/
void fmpz_mpoly_ksub_content(
    fmpz_t content,
    const fmpz_mpoly_t A,
    const ulong * subdegs,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    fmpz_mpoly_t T;
    fmpz_mpoly_ctx_t Tctx;
    fmpz_t e;
    fmpz * Acoeff = A->coeffs;
    ulong * Aexp = A->exps;
    ulong mask;
    slong * offsets, * shifts;
    TMP_INIT;

    TMP_START;
    fmpz_init(e);

    fmpz_mpoly_ctx_init(Tctx, 1, ORD_LEX);
    fmpz_mpoly_init(T, Tctx);

    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    offsets = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    shifts = (slong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(slong));
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, ctx->minfo);
    }

    for (i = 0; i < A->length; i++)
    {
        fmpz_zero(e);
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            fmpz_mul_ui(e, e, subdegs[j]);
            fmpz_add_ui(e, e, ((Aexp + N*i)[offsets[j]]>>shifts[j])&mask);
        }
        _fmpz_mpoly_push_exp_ffmpz(T, e, Tctx);
        fmpz_set(T->coeffs + T->length - 1, Acoeff + i);
    }
    fmpz_mpoly_sort_terms(T, Tctx);
    fmpz_mpoly_combine_like_terms(T, Tctx);
    _fmpz_vec_content(content, T->coeffs, T->length);
    fmpz_mpoly_clear(T, Tctx);
    fmpz_mpoly_ctx_clear(Tctx);

    fmpz_clear(e);
    TMP_END;
}




/*
    return 1: good (or at least not bad)
           0: no match
          -1: better degree bound is in GevaldegXY
*/

int static _random_check_sp(
    ulong * GevaldegXY,
    ulong GdegboundXY,
    int which_check,
    n_bpoly_t Aeval_sp,
    n_bpoly_t Beval_sp,
    n_bpoly_t Geval_sp,
    n_bpoly_t Abareval_sp,
    n_bpoly_t Bbareval_sp,
    mp_limb_t * alphas_sp,
    const fmpz_mpolyu_t H,
    const fmpz_mpolyu_t A,
    const fmpz_mpolyu_t B,
    const fmpz_mpoly_t Gamma,
    const fmpz_mpoly_ctx_t ctx,
    const nmod_mpoly_ctx_t ctx_sp,
    flint_rand_t randstate,
    n_poly_bpoly_stack_t Sp)
{
    mp_limb_t Gammaeval_sp;
    int success;
    int point_try_count;
    slong i;

    for (point_try_count = 0; point_try_count < 10; point_try_count++)
    {
        for (i = 0; i < ctx->minfo->nvars; i++)
            alphas_sp[i] = n_urandint(randstate, ctx_sp->mod.n);
        fmpz_mpolyuu_eval_nmod(Aeval_sp, ctx_sp, A, alphas_sp, ctx);
        fmpz_mpolyuu_eval_nmod(Beval_sp, ctx_sp, B, alphas_sp, ctx);

        if (Aeval_sp->length < 1 || Beval_sp->length < 1 ||
            n_bpoly_bidegree(Aeval_sp) != A->exps[0] ||
            n_bpoly_bidegree(Beval_sp) != B->exps[0])
        {
            continue;
        }

        Gammaeval_sp = fmpz_mpoly_eval_nmod(ctx_sp->mod, Gamma, alphas_sp, ctx);
        FLINT_ASSERT(Gammaeval_sp != 0);

        success = n_bpoly_mod_gcd_brown_smprime(Geval_sp, Abareval_sp,
                             Bbareval_sp, Aeval_sp, Beval_sp, ctx_sp->mod, Sp);
        if (!success)
            continue;

        n_bpoly_scalar_mul_nmod(Geval_sp, Gammaeval_sp, ctx_sp->mod);

        FLINT_ASSERT(Geval_sp->length > 0);
        *GevaldegXY = n_bpoly_bidegree(Geval_sp);

        if (GdegboundXY < *GevaldegXY)
        {
            continue;
        }
        else if (GdegboundXY > *GevaldegXY)
        {
            return -1;
        }

        if (which_check == 0)
        {
            fmpz_mpolyuu_eval_nmod(Bbareval_sp, ctx_sp, H, alphas_sp, ctx);
            return n_bpoly_equal(Bbareval_sp, Geval_sp);
        }
        else
        {
            fmpz_mpolyuu_eval_nmod(Geval_sp, ctx_sp, H, alphas_sp, ctx);
            return n_bpoly_equal(Geval_sp,
                                 which_check == 1 ? Abareval_sp : Bbareval_sp);
        }
    }

    return 1; /* Hmm */
}

int static _random_check_mp(
    ulong * GevaldegXY,
    ulong GdegboundXY,
    int which_check,
    fmpz_mod_mpolyn_t Aeval,
    fmpz_mod_mpolyn_t Beval,
    fmpz_mod_mpolyn_t Geval,
    fmpz_mod_mpolyn_t Abareval,
    fmpz_mod_mpolyn_t Bbareval,
    fmpz_t Gammaeval,
    fmpz * alphas,
    const fmpz_mpolyu_t H,
    const fmpz_mpolyu_t A,
    const fmpz_mpolyu_t B,
    const fmpz_mpoly_t Gamma,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_mpoly_ctx_t ctx_mp,
    flint_rand_t randstate)
{
    int success;
    int point_try_count;
    slong i;

    /* try to test H at a random evaluation point */
    for (point_try_count = 0; point_try_count < 10; point_try_count++)
    {
        for (i = 0; i < ctx->minfo->nvars; i++)
            fmpz_randm(alphas + i, randstate, fmpz_mod_ctx_modulus(ctx_mp->ffinfo));
        fmpz_mpolyuu_eval_fmpz_mod(Aeval, ctx_mp, A, alphas, ctx);
        fmpz_mpolyuu_eval_fmpz_mod(Beval, ctx_mp, B, alphas, ctx);

        if (Aeval->length == 0 || Beval->length == 0 ||
            fmpz_mod_mpolyn_bidegree(Aeval) != A->exps[0] ||
            fmpz_mod_mpolyn_bidegree(Beval) != B->exps[0])
        {
            continue;
        }

        fmpz_mpoly_eval_fmpz_mod(Gammaeval, ctx_mp->ffinfo, Gamma, alphas, ctx);
        FLINT_ASSERT(!fmpz_is_zero(Gammaeval));

        success = fmpz_mod_mpolyn_gcd_brown_bivar(Geval, Abareval, Bbareval,
                                                         Aeval, Beval, ctx_mp);
        if (!success)
            continue;

        fmpz_mod_mpolyn_scalar_mul_fmpz_mod(Geval, Gammaeval, ctx_mp);

        FLINT_ASSERT(Geval->length > 0);
        *GevaldegXY = fmpz_mod_mpolyn_bidegree(Geval);

        if (GdegboundXY < *GevaldegXY)
        {
            continue;
        }
        else if (GdegboundXY > *GevaldegXY)
        {
            return -1;
        }

        if (which_check == 0)
        {
            fmpz_mpolyuu_eval_fmpz_mod(Bbareval, ctx_mp, H, alphas, ctx);
            return fmpz_mod_mpolyn_equal(Bbareval, Geval, ctx_mp);
        }
        else
        {
            fmpz_mpolyuu_eval_fmpz_mod(Geval, ctx_mp, H, alphas, ctx);
            return fmpz_mod_mpolyn_equal(Geval,
                               which_check == 1 ? Abareval : Bbareval, ctx_mp);
        }
    }

    return 1; /* Hmm */
}



static void nmod_mpoly_get_eval_helper(
    n_poly_t Acur,
    n_poly_t Ainc,
    const nmod_mpoly_ctx_t ctx_sp,
    const fmpz_mpoly_t A,
    const mp_limb_t * alphas,
    const fmpz_mpoly_ctx_t ctx)
{
    mp_limb_t * p, * q;
    const fmpz * Acoeffs = A->coeffs;
    slong j, Alen = A->length;

    n_poly_fit_length(Acur, Alen);
    Acur->length = Alen;

    n_poly_fit_length(Ainc, 2*Alen);
    Ainc->length = 2*Alen;

    p = Acur->coeffs;
    _mpoly_monomial_evals_nmod(p, A->exps, A->bits, Alen, alphas, 0,
                                                      ctx->minfo, ctx_sp->mod);
    q = Ainc->coeffs;
    for (j = 0; j < Alen; j++)
    {
        q[2*j + 0] = fmpz_fdiv_ui(Acoeffs + j, ctx_sp->mod.n);
        q[2*j + 1] = p[j];
    }
}

static void nmod_mpolyu_get_eval_helper(
    n_bpoly_t Acur,
    n_polyun_t Ainc,
    const nmod_mpoly_ctx_t ctx_sp,
    const fmpz_mpolyu_t A,
    const mp_limb_t * alphas,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    n_bpoly_fit_length(Acur, A->length);
    Acur->length = A->length;

    n_polyun_fit_length(Ainc, A->length);
    Ainc->length = A->length;

    for (i = 0; i < A->length; i++)
    {
        Ainc->terms[i].exp = A->exps[i];
        nmod_mpoly_get_eval_helper(Acur->coeffs + i, Ainc->terms[i].coeff,
                                           ctx_sp, A->coeffs + i, alphas, ctx);
    }
}


int fmpz_mpolyuu_gcd_berlekamp_massey(
    fmpz_mpolyu_t G,
    fmpz_mpolyu_t Abar,
    fmpz_mpolyu_t Bbar,
    fmpz_mpolyu_t A,
    fmpz_mpolyu_t B,
    const fmpz_mpoly_t Gamma,
    const fmpz_mpoly_ctx_t ctx)
{
    int which_check, changed, success, point_try_count;
    flint_bitcnt_t bits = A->bits;
    mpoly_bma_interpolate_ctx_t Ictx;
    nmod_mpoly_ctx_t ctx_sp;
    fmpz_mod_mpoly_ctx_t ctx_mp;
    fmpz_mpolyu_t H;
    fmpz_mpoly_t Hcontent;
    /* multi precision workspace */
    fmpz_mod_bma_mpoly_t GLambda_mp, AbarLambda_mp, BbarLambda_mp;
    fmpz_mod_mpolyn_t Aeval_mp, Beval_mp, Geval_mp, Abareval_mp, Bbareval_mp;
    fmpz_mpolycu_t Ainc_mp, Acur_mp, Binc_mp, Bcur_mp, Ared_mp, Bred_mp;
    fmpz_mpolyc_t Gammainc_mp, Gammacur_mp, Gammared_mp;
    fmpz_t p, pm1, sshift_mp, last_unlucky_sshift_plus_1_mp, image_count_mp;
    fmpz_t Gammaeval_mp;
    fmpz * alphas_mp;
    /* single precision workspace */
    n_poly_bpoly_stack_t Sp;
    nmod_bma_mpoly_t GLambda_sp, AbarLambda_sp, BbarLambda_sp;
    n_bpoly_t Aeval_sp, Beval_sp, Geval_sp, Abareval_sp, Bbareval_sp;
    n_poly_t Gammacur_sp;
    n_bpoly_t Acur_sp, Bcur_sp;
    n_poly_t Gammainc_sp;
    n_polyun_t Ainc_sp, Binc_sp;
    mp_limb_t p_sp, sshift_sp, last_unlucky_sshift_plus_1_sp, image_count_sp;
    mp_limb_t Gammaeval_sp;
    mp_limb_t * alphas_sp;
    /* misc */
    slong i, j;
    ulong GdegboundXY, GevaldegXY;
    slong * Adegs, * Bdegs;
    flint_rand_t randstate;
    fmpz_t subprod, cAksub, cBksub;
    int unlucky_count;
    fmpz_t Hmodulus;
    nmod_zip_mpolyu_t newZ;
    slong zip_evals;
    ulong ABtotal_length;

    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == G->bits);
    FLINT_ASSERT(bits == Abar->bits);
    FLINT_ASSERT(bits == Bbar->bits);
    FLINT_ASSERT(bits == Gamma->bits);

    /* let's initialize everything at once to avoid complicated cleanup */

    ABtotal_length = 0;
    for (i = 0; i < A->length; i++)
        ABtotal_length += A->coeffs[i].length;
    for (i = 0; i < B->length; i++)
        ABtotal_length += B->coeffs[i].length;

    ABtotal_length = FLINT_MAX(ABtotal_length, 100);

    flint_randinit(randstate);
    fmpz_init(p);
    fmpz_init(pm1); /* p - 1 */
    fmpz_init(image_count_mp);
    fmpz_init(subprod);
    fmpz_init(cAksub);
    fmpz_init(cBksub);
    fmpz_init(sshift_mp);
    fmpz_init(last_unlucky_sshift_plus_1_mp);
    fmpz_init(Hmodulus);
    fmpz_mpolyu_init(H, bits, ctx);
    fmpz_mpoly_init3(Hcontent, 0, bits, ctx);

    mpoly_bma_interpolate_ctx_init(Ictx, ctx->minfo->nvars);

    /* multiprecision workspace */
    fmpz_set_ui(p, 2);    /* modulus no care */
    fmpz_mod_mpoly_ctx_init(ctx_mp, 2, ORD_LEX, p); /* modulus no care */
    fmpz_mod_bma_mpoly_init(GLambda_mp);
    fmpz_mod_bma_mpoly_init(AbarLambda_mp);
    fmpz_mod_bma_mpoly_init(BbarLambda_mp);
    fmpz_init(Gammaeval_mp);
    fmpz_mod_mpolyn_init(Aeval_mp, FLINT_BITS/2, ctx_mp);
    fmpz_mod_mpolyn_init(Beval_mp, FLINT_BITS/2, ctx_mp);
    fmpz_mod_mpolyn_init(Geval_mp, FLINT_BITS/2, ctx_mp);
    fmpz_mod_mpolyn_init(Abareval_mp, FLINT_BITS/2, ctx_mp);
    fmpz_mod_mpolyn_init(Bbareval_mp, FLINT_BITS/2, ctx_mp);
    fmpz_mpolyc_init(Gammainc_mp);
    fmpz_mpolyc_init(Gammacur_mp);
    fmpz_mpolyc_init(Gammared_mp);
    fmpz_mpolycu_init(Ainc_mp);
    fmpz_mpolycu_init(Acur_mp);
    fmpz_mpolycu_init(Ared_mp);
    fmpz_mpolycu_init(Binc_mp);
    fmpz_mpolycu_init(Bcur_mp);
    fmpz_mpolycu_init(Bred_mp);
    alphas_mp = _fmpz_vec_init(ctx->minfo->nvars);

    /* machine precision workspace "sp" */
    nmod_mpoly_ctx_init(ctx_sp, 2, ORD_LEX, 2); /* modulus no care */

    n_poly_stack_init(Sp->poly_stack);
    n_bpoly_stack_init(Sp->bpoly_stack);

    nmod_bma_mpoly_init(GLambda_sp);
    nmod_bma_mpoly_init(AbarLambda_sp);
    nmod_bma_mpoly_init(BbarLambda_sp);

    n_bpoly_init(Aeval_sp);
    n_bpoly_init(Beval_sp);
    n_bpoly_init(Geval_sp);
    n_bpoly_init(Abareval_sp);
    n_bpoly_init(Bbareval_sp);

    n_poly_init(Gammacur_sp);
    n_bpoly_init(Acur_sp);
    n_bpoly_init(Bcur_sp);
    n_poly_init(Gammainc_sp);
    n_polyun_init(Ainc_sp);
    n_polyun_init(Binc_sp);

    alphas_sp = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));

    nmod_zip_mpolyu_init(newZ);

    Adegs      = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    Bdegs      = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));

    GdegboundXY = FLINT_MIN(A->exps[0], B->exps[0]);
    p_sp = UWORD(1) << (FLINT_BITS - 2);
    for (point_try_count = 0; point_try_count < 10; point_try_count++)
    {
        p_sp = n_nextprime(p_sp, 1);
        nmod_mpoly_ctx_set_modulus(ctx_sp, p_sp);

        for (i = 0; i < ctx->minfo->nvars; i++)
            alphas_sp[i] = n_urandint(randstate, p_sp);
        fmpz_mpolyuu_eval_nmod(Aeval_sp, ctx_sp, A, alphas_sp, ctx);
        fmpz_mpolyuu_eval_nmod(Beval_sp, ctx_sp, B, alphas_sp, ctx);

        if (Aeval_sp->length < 1 || Beval_sp->length < 1 ||
            n_bpoly_bidegree(Aeval_sp) != A->exps[0] ||
            n_bpoly_bidegree(Beval_sp) != B->exps[0])
        {
            continue;
        }
        success = n_bpoly_mod_gcd_brown_smprime(Geval_sp, Abareval_sp,
                             Bbareval_sp, Aeval_sp, Beval_sp, ctx_sp->mod, Sp);
        if (success)
        {
            FLINT_ASSERT(Geval_sp->length > 0);
            GdegboundXY = n_bpoly_bidegree(Geval_sp);
            break;
        }
    }

    fmpz_mpolyu_inner_degrees_si(Adegs, A, ctx);
    fmpz_mpolyu_inner_degrees_si(Bdegs, B, ctx);
    for (i = 0; i < ctx->minfo->nvars; i++)
        Ictx->degbounds[i] = 1 + FLINT_MAX(Adegs[i], Bdegs[i]);

    /* initialization done! */

    if (GdegboundXY == 0)
    {
        fmpz_mpolyu_one(G, ctx);
        fmpz_mpolyu_swap(Abar, A, ctx);
        fmpz_mpolyu_swap(Bbar, B, ctx);
        success = 1;
        goto cleanup;
    }

    for (i = 0; i < ctx->minfo->nvars; i++)
        Ictx->subdegs[i] = Ictx->degbounds[i];

    goto got_ksub;

pick_ksub:

    if (ctx->minfo->nvars > 1)
    {
        /* just increment the smallest subdegs[j] */
        j = 1;
        for (i = 2; i < ctx->minfo->nvars; i++)
        {
            if (Ictx->subdegs[i] < Ictx->subdegs[j])
            {
                j = i;
            }
        }
        Ictx->subdegs[j]++;
    }

got_ksub:

    fmpz_one(subprod);
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        if ((slong)(Ictx->subdegs[i]) < 0)
        {
            success = 0;
            goto cleanup;
        }
        fmpz_mul_ui(subprod, subprod, Ictx->subdegs[i]);
    }

    /* see if the ksub killed either lc(A) or lc(B) */
    fmpz_mpoly_ksub_content(cAksub, A->coeffs + 0, Ictx->subdegs, ctx);
    fmpz_mpoly_ksub_content(cBksub, B->coeffs + 0, Ictx->subdegs, ctx);
    if (fmpz_is_zero(cAksub) || fmpz_is_zero(cBksub))
        goto pick_ksub;

pick_bma_prime:

    if (fmpz_cmp_ui(p, ABtotal_length) < 0)
        fmpz_set_ui(p, ABtotal_length);

    if (fmpz_cmp(p, subprod) < 0)
        fmpz_set(p, subprod);

    success = fmpz_next_smooth_prime(p, p);
    fmpz_sub_ui(pm1, p, 1);
    if (!success)
    {
        /* ran out of smooth primes - absolute falure */
        success = 0;
        goto cleanup;
    }

    /* make sure reduction does not kill either leading coeff after ksub */
    if (fmpz_divisible(cAksub, p) || fmpz_divisible(cBksub, p))
    {
        goto pick_bma_prime;
    }

    /* make sure p does not divide any coefficient of Gamma */
    for (i = 0; i < Gamma->length; i++)
    {
        if (fmpz_divisible(Gamma->coeffs + i, p))
        {
            goto pick_bma_prime;
        }
    }

    if (fmpz_abs_fits_ui(p))
    {
        p_sp = fmpz_get_ui(p);
        sshift_sp = 1;

        unlucky_count = 0;
        last_unlucky_sshift_plus_1_sp = 0;

        nmod_mpoly_ctx_set_modulus(ctx_sp, p_sp);
        nmod_discrete_log_pohlig_hellman_precompute_prime(Ictx->dlogenv_sp, p_sp);

        nmod_bma_mpoly_reset_prime(GLambda_sp, ctx_sp->mod);
        nmod_bma_mpoly_reset_prime(AbarLambda_sp, ctx_sp->mod);
        nmod_bma_mpoly_reset_prime(BbarLambda_sp, ctx_sp->mod);
        nmod_bma_mpoly_zero(GLambda_sp);
        nmod_bma_mpoly_zero(AbarLambda_sp);
        nmod_bma_mpoly_zero(BbarLambda_sp);

        FLINT_ASSERT(sshift_sp == 1);
        nmod_mpoly_bma_interpolate_alpha_powers(alphas_sp, sshift_sp,
                                                       Ictx, ctx, ctx_sp->mod);

        nmod_mpoly_get_eval_helper(Gammacur_sp, Gammainc_sp, ctx_sp, Gamma, alphas_sp, ctx);
        nmod_mpolyu_get_eval_helper(Acur_sp, Ainc_sp, ctx_sp, A, alphas_sp, ctx);
        nmod_mpolyu_get_eval_helper(Bcur_sp, Binc_sp, ctx_sp, B, alphas_sp, ctx);

        image_count_sp = 0;

    next_bma_image_sp:

        /* image count is also the current power of alpha we are evaluating */
        image_count_sp++;

        FLINT_ASSERT(GLambda_sp->pointcount == AbarLambda_sp->pointcount);
        FLINT_ASSERT(GLambda_sp->pointcount == BbarLambda_sp->pointcount);
        FLINT_ASSERT(sshift_sp + GLambda_sp->pointcount == image_count_sp);

        if (image_count_sp >= p_sp - 1)
            goto pick_bma_prime;

        Gammaeval_sp = n_poly_mod_eval_step2(Gammacur_sp, Gammainc_sp, ctx_sp->mod);
        nmod_mpolyuu_eval_step2(Aeval_sp, Acur_sp, Ainc_sp, ctx_sp);
        nmod_mpolyuu_eval_step2(Beval_sp, Bcur_sp, Binc_sp, ctx_sp);

        if (Aeval_sp->length < 1 || Beval_sp->length < 1 ||
            n_bpoly_bidegree(Aeval_sp) != A->exps[0] ||
            n_bpoly_bidegree(Beval_sp) != B->exps[0])
        {
            sshift_sp += GLambda_sp->pointcount + 1;
            nmod_bma_mpoly_zero(GLambda_sp);
            nmod_bma_mpoly_zero(AbarLambda_sp);
            nmod_bma_mpoly_zero(BbarLambda_sp);
            goto next_bma_image_sp;
        }

        /* the evaluation killed neither lc(A) nor lc(B) */
        FLINT_ASSERT(Gammaeval_sp != 0);

        success = n_bpoly_mod_gcd_brown_smprime(Geval_sp, Abareval_sp,
                             Bbareval_sp, Aeval_sp, Beval_sp, ctx_sp->mod, Sp);
        if (!success)
        {
            sshift_sp += GLambda_sp->pointcount + 1;
            nmod_bma_mpoly_zero(GLambda_sp);
            nmod_bma_mpoly_zero(AbarLambda_sp);
            nmod_bma_mpoly_zero(BbarLambda_sp);
            goto next_bma_image_sp;
        }

        FLINT_ASSERT(Geval_sp->length > 0);
        GevaldegXY = n_bpoly_bidegree(Geval_sp);
        n_bpoly_scalar_mul_nmod(Geval_sp, Gammaeval_sp, ctx_sp->mod);

        if (GdegboundXY < GevaldegXY)
        {
            if (sshift_sp == last_unlucky_sshift_plus_1_sp)
                goto pick_ksub;

            if (++unlucky_count > 2)
                goto pick_bma_prime;

            last_unlucky_sshift_plus_1_sp = sshift_sp + 1;
            sshift_sp += GLambda_sp->pointcount + 1;
            nmod_bma_mpoly_zero(GLambda_sp);
            nmod_bma_mpoly_zero(AbarLambda_sp);
            nmod_bma_mpoly_zero(BbarLambda_sp);
            goto next_bma_image_sp;        
        }
        else if (GdegboundXY > GevaldegXY)
        {
            /* new bound on deg_XY(G) */
            GdegboundXY = GevaldegXY;
            if (GdegboundXY == 0)
            {
                fmpz_mpolyu_one(G, ctx);
                fmpz_mpolyu_swap(Abar, A, ctx);
                fmpz_mpolyu_swap(Bbar, B, ctx);
                success = 1;
                goto cleanup;
            }
            sshift_sp += GLambda_sp->pointcount;
            nmod_bma_mpoly_zero(GLambda_sp);
            nmod_bma_mpoly_zero(AbarLambda_sp);
            nmod_bma_mpoly_zero(BbarLambda_sp);
            nmod_bma_mpoly_add_point(GLambda_sp, Geval_sp, ctx_sp);
            nmod_bma_mpoly_add_point(AbarLambda_sp, Abareval_sp, ctx_sp);
            nmod_bma_mpoly_add_point(BbarLambda_sp, Bbareval_sp, ctx_sp);
            goto next_bma_image_sp;
        }

        nmod_bma_mpoly_add_point(GLambda_sp, Geval_sp, ctx_sp);
        nmod_bma_mpoly_add_point(AbarLambda_sp, Abareval_sp, ctx_sp);
        nmod_bma_mpoly_add_point(BbarLambda_sp, Bbareval_sp, ctx_sp);

        if ((GLambda_sp->pointcount & 7) != 0)
        {
            goto next_bma_image_sp;
        }

        if (GLambda_sp->pointcount/2 >= Gamma->length &&
            !nmod_bma_mpoly_reduce(GLambda_sp) &&
            nmod_bma_mpoly_get_fmpz_mpolyu(H, ctx, sshift_sp,
                                           GLambda_sp, Ictx, ctx_sp->mod) &&
            H->length > 0 && H->coeffs[0].length == Gamma->length)
        {
            which_check = 0;
            goto check_sp;
        }

        if (AbarLambda_sp->pointcount/2 >= A->coeffs[0].length &&
            !nmod_bma_mpoly_reduce(AbarLambda_sp) &&
            nmod_bma_mpoly_get_fmpz_mpolyu(H, ctx, sshift_sp,
                                           AbarLambda_sp, Ictx, ctx_sp->mod) &&
            H->length > 0 && H->coeffs[0].length == A->coeffs[0].length)
        {
            which_check = 1;
            goto check_sp;
        }

        if (BbarLambda_sp->pointcount/2 >= B->coeffs[0].length &&
            !nmod_bma_mpoly_reduce(BbarLambda_sp) &&
            nmod_bma_mpoly_get_fmpz_mpolyu(H, ctx, sshift_sp,
                                           BbarLambda_sp, Ictx, ctx_sp->mod) &&
            H->length > 0 && H->coeffs[0].length == B->coeffs[0].length)
        {
            which_check = 2;
            goto check_sp;
        }

        if (GLambda_sp->pointcount/2 > ABtotal_length)
        {
            success = 0;
            goto cleanup;
        }

        goto next_bma_image_sp;

    check_sp:

        FLINT_ASSERT(0 <= which_check && which_check <= 2);
        FLINT_ASSERT(which_check != 0 || GdegboundXY == H->exps[0]);

        success =  _random_check_sp(&GevaldegXY, GdegboundXY, which_check,
                Aeval_sp, Beval_sp, Geval_sp, Abareval_sp, Bbareval_sp,
                alphas_sp, H, A, B, Gamma, ctx, ctx_sp, randstate, Sp);
        if (success < 1)
        {
            if (success == 0)
                goto next_bma_image_sp;

            FLINT_ASSERT(GdegboundXY > GevaldegXY);
            GdegboundXY = GevaldegXY;
            if (GdegboundXY == 0)
            {
                fmpz_mpolyu_one(G, ctx);
                fmpz_mpolyu_swap(Abar, A, ctx);
                fmpz_mpolyu_swap(Bbar, B, ctx);
                success = 1;
                goto cleanup;
            }
            sshift_sp += GLambda_sp->pointcount;
            nmod_bma_mpoly_zero(GLambda_sp);
            nmod_bma_mpoly_zero(AbarLambda_sp);
            nmod_bma_mpoly_zero(BbarLambda_sp);
            goto next_bma_image_sp;
        }
    }
    else
    {
        fmpz_one(sshift_mp);

        unlucky_count = 0;
        fmpz_zero(last_unlucky_sshift_plus_1_mp);

        fmpz_mod_ctx_set_modulus(ctx_mp->ffinfo, p);
        fmpz_mod_discrete_log_pohlig_hellman_precompute_prime(Ictx->dlogenv, p);
        fmpz_mod_bma_mpoly_zero(GLambda_mp);
        fmpz_mod_bma_mpoly_zero(AbarLambda_mp);
        fmpz_mod_bma_mpoly_zero(BbarLambda_mp);

        FLINT_ASSERT(fmpz_is_one(sshift_mp));
        fmpz_mod_mpoly_bma_interpolate_alpha_powers(alphas_mp, sshift_mp,
                                                    Ictx, ctx, ctx_mp->ffinfo);

        fmpz_mod_mpoly_set_skel(Gammainc_mp, ctx_mp, Gamma, alphas_mp, ctx);
        fmpz_mod_mpolyu_set_skel(Ainc_mp, ctx_mp, A, alphas_mp, ctx);
        fmpz_mod_mpolyu_set_skel(Binc_mp, ctx_mp, B, alphas_mp, ctx);

        fmpz_mod_mpoly_red_skel(Gammared_mp, Gamma, ctx_mp->ffinfo);
        fmpz_mod_mpolyu_red_skel(Ared_mp, A, ctx_mp->ffinfo);
        fmpz_mod_mpolyu_red_skel(Bred_mp, B, ctx_mp->ffinfo);

        fmpz_mod_mpoly_copy_skel(Gammacur_mp, Gammainc_mp);
        fmpz_mod_mpolyu_copy_skel(Acur_mp, Ainc_mp);
        fmpz_mod_mpolyu_copy_skel(Bcur_mp, Binc_mp);

        fmpz_zero(image_count_mp);

    next_bma_image_mp:

        /* image count is also the current power of alpha we are evaluating */
        fmpz_add_ui(image_count_mp, image_count_mp, 1);

        FLINT_ASSERT(GLambda_sp->pointcount == AbarLambda_sp->pointcount);
        FLINT_ASSERT(GLambda_sp->pointcount == BbarLambda_sp->pointcount);
    #if FLINT_WANT_ASSERT
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_add_ui(t, sshift_mp, GLambda_mp->pointcount);
            FLINT_ASSERT(fmpz_equal(t, image_count_mp));
            fmpz_clear(t);
        }
    #endif

        if (fmpz_cmp(image_count_mp, pm1) >= 0)
            goto pick_bma_prime;

        fmpz_mod_mpoly_use_skel_mul(Gammaeval_mp, Gammared_mp, Gammacur_mp, Gammainc_mp,
                                                               ctx_mp->ffinfo);
        fmpz_mod_mpolyuu_use_skel_mul(Aeval_mp, A, Ared_mp, Acur_mp, Ainc_mp, ctx_mp);
        fmpz_mod_mpolyuu_use_skel_mul(Beval_mp, B, Bred_mp, Bcur_mp, Binc_mp, ctx_mp);
        if (Aeval_mp->length == 0 || Beval_mp->length == 0 ||
            fmpz_mod_mpolyn_bidegree(Aeval_mp) != A->exps[0] ||
            fmpz_mod_mpolyn_bidegree(Beval_mp) != B->exps[0])
        {
            fmpz_add_ui(sshift_mp, sshift_mp, GLambda_mp->pointcount + 1);
            fmpz_mod_bma_mpoly_zero(GLambda_mp);
            fmpz_mod_bma_mpoly_zero(AbarLambda_mp);
            fmpz_mod_bma_mpoly_zero(BbarLambda_mp);
            goto next_bma_image_mp;
        }

        /* the evaluation killed neither lc(A) nor lc(B) */
        FLINT_ASSERT(!fmpz_is_zero(Gammaeval_mp));

        success = fmpz_mod_mpolyn_gcd_brown_bivar(Geval_mp,
                         Abareval_mp, Bbareval_mp, Aeval_mp, Beval_mp, ctx_mp);
        if (!success)
        {
            fmpz_add_ui(sshift_mp, sshift_mp, GLambda_mp->pointcount + 1);
            fmpz_mod_bma_mpoly_zero(GLambda_mp);
            fmpz_mod_bma_mpoly_zero(AbarLambda_mp);
            fmpz_mod_bma_mpoly_zero(BbarLambda_mp);
            goto next_bma_image_mp;
        }

        FLINT_ASSERT(Geval_mp->length > 0);
        GevaldegXY = fmpz_mod_mpolyn_bidegree(Geval_mp);
        fmpz_mod_mpolyn_scalar_mul_fmpz_mod(Geval_mp, Gammaeval_mp, ctx_mp);

        FLINT_ASSERT(fmpz_equal(Gammaeval_mp,
                        fmpz_mod_mpolyn_leadcoeff_last_ref(Geval_mp, ctx_mp)));

        if (GdegboundXY < GevaldegXY)
        {
            if (fmpz_equal(sshift_mp, last_unlucky_sshift_plus_1_mp))
                goto pick_ksub;

            if (++unlucky_count > 2)
                goto pick_bma_prime;

            fmpz_add_ui(last_unlucky_sshift_plus_1_mp, sshift_mp, 1);
            fmpz_add_ui(sshift_mp, sshift_mp, GLambda_mp->pointcount + 1);
            fmpz_mod_bma_mpoly_zero(GLambda_mp);
            fmpz_mod_bma_mpoly_zero(AbarLambda_mp);
            fmpz_mod_bma_mpoly_zero(BbarLambda_mp);
            goto next_bma_image_mp;
        }
        else if (GdegboundXY > GevaldegXY)
        {
            GdegboundXY = GevaldegXY;
            if (GdegboundXY == 0)
            {
                fmpz_mpolyu_one(G, ctx);
                fmpz_mpolyu_swap(Abar, A, ctx);
                fmpz_mpolyu_swap(Bbar, B, ctx);
                success = 1;
                goto cleanup;
            }
            fmpz_add_ui(sshift_mp, sshift_mp, GLambda_mp->pointcount);
            fmpz_mod_bma_mpoly_zero(GLambda_mp);
            fmpz_mod_bma_mpoly_zero(AbarLambda_mp);
            fmpz_mod_bma_mpoly_zero(BbarLambda_mp);
            fmpz_mod_bma_mpoly_add_point(GLambda_mp, Geval_mp, ctx_mp);
            fmpz_mod_bma_mpoly_add_point(AbarLambda_mp, Abareval_mp, ctx_mp);
            fmpz_mod_bma_mpoly_add_point(BbarLambda_mp, Bbareval_mp, ctx_mp);
            goto next_bma_image_mp;
        }

        fmpz_mod_bma_mpoly_add_point(GLambda_mp, Geval_mp, ctx_mp);
        fmpz_mod_bma_mpoly_add_point(AbarLambda_mp, Abareval_mp, ctx_mp);
        fmpz_mod_bma_mpoly_add_point(BbarLambda_mp, Bbareval_mp, ctx_mp);

        if ((GLambda_mp->pointcount & 7) != 0)
        {
            goto next_bma_image_mp;
        }

        if (GLambda_mp->pointcount/2 >= Gamma->length &&
            !fmpz_mod_bma_mpoly_reduce(GLambda_mp, ctx_mp->ffinfo) &&
            fmpz_mod_bma_mpoly_get_fmpz_mpolyu(H, ctx, sshift_mp, GLambda_mp,
                                                      Ictx, ctx_mp->ffinfo) &&
            H->length > 0 && H->coeffs[0].length == Gamma->length)
        {
            which_check = 0;
            goto check_mp;
        }

        if (AbarLambda_mp->pointcount/2 >= A->coeffs[0].length &&
            !fmpz_mod_bma_mpoly_reduce(AbarLambda_mp, ctx_mp->ffinfo) &&
            fmpz_mod_bma_mpoly_get_fmpz_mpolyu(H, ctx, sshift_mp, AbarLambda_mp,
                                                      Ictx, ctx_mp->ffinfo) &&
            H->length > 0 && H->coeffs[0].length == A->coeffs[0].length)
        {
            which_check = 1;
            goto check_mp;
        }

        if (BbarLambda_mp->pointcount/2 >= B->coeffs[0].length &&
            !fmpz_mod_bma_mpoly_reduce(BbarLambda_mp, ctx_mp->ffinfo) &&
            fmpz_mod_bma_mpoly_get_fmpz_mpolyu(H, ctx, sshift_mp, BbarLambda_mp,
                                                      Ictx, ctx_mp->ffinfo) &&
            H->length > 0 && H->coeffs[0].length == B->coeffs[0].length)
        {
            which_check = 2;
            goto check_mp;
        }

        if (GLambda_mp->pointcount/2 > ABtotal_length)
        {
            success = 0;
            goto cleanup;
        }

        goto next_bma_image_mp;

    check_mp:

        FLINT_ASSERT(0 <= which_check && which_check <= 2);
        FLINT_ASSERT(which_check != 0 || GdegboundXY == H->exps[0]);

        success = _random_check_mp(&GevaldegXY, GdegboundXY, which_check,
                    Aeval_mp, Beval_mp, Geval_mp, Abareval_mp, Bbareval_mp, Gammaeval_mp,
                        alphas_mp, H, A, B, Gamma, ctx, ctx_mp, randstate);
        if (success < 1)
        {
            if (success == 0)
                goto next_bma_image_mp;

            FLINT_ASSERT(GdegboundXY > GevaldegXY);
            GdegboundXY = GevaldegXY;
            if (GdegboundXY == 0)
            {
                fmpz_mpolyu_one(G, ctx);
                fmpz_mpolyu_swap(Abar, A, ctx);
                fmpz_mpolyu_swap(Bbar, B, ctx);
                success = 1;
                goto cleanup;
            }
            fmpz_add_ui(sshift_mp, sshift_mp, GLambda_mp->pointcount);
            fmpz_mod_bma_mpoly_zero(GLambda_mp);
            fmpz_mod_bma_mpoly_zero(AbarLambda_mp);
            fmpz_mod_bma_mpoly_zero(BbarLambda_mp);
            goto next_bma_image_mp;
        }
    }

    /* assume that H is correct modulo Hmodulus = p */
    FLINT_ASSERT(0 <= which_check && which_check <= 2);
    fmpz_set(Hmodulus, p);

    FLINT_ASSERT(H->length > 0);
    zip_evals = H->coeffs[0].length;
    for (i = 1; i < H->length; i++)
        zip_evals = FLINT_MAX(zip_evals, H->coeffs[i].length);
    zip_evals += 1;
    nmod_zip_mpolyu_fit_poly(newZ, H, zip_evals);

    p_sp = UWORD(1) << (FLINT_BITS - 2);

pick_zip_prime:

    if (p_sp >= UWORD_MAX_PRIME)
    {
        success = 0;
        goto cleanup;
    }
    p_sp = n_nextprime(p_sp, 1);

    if (0 == fmpz_fdiv_ui(Hmodulus, p_sp))
        goto pick_zip_prime;

    nmod_mpoly_ctx_set_modulus(ctx_sp, p_sp);

    FLINT_ASSERT(p_sp > 3);
    for (i = 0; i < ctx->minfo->nvars; i++)
        alphas_sp[i] = n_urandint(randstate, p_sp - 3) + 2;

    nmod_zip_mpolyu_set_skel(newZ, ctx_sp, H, alphas_sp, ctx);
    nmod_mpoly_get_eval_helper(Gammacur_sp, Gammainc_sp, ctx_sp, Gamma, alphas_sp, ctx);
    nmod_mpolyu_get_eval_helper(Acur_sp, Ainc_sp, ctx_sp, A, alphas_sp, ctx);
    nmod_mpolyu_get_eval_helper(Bcur_sp, Binc_sp, ctx_sp, B, alphas_sp, ctx);

next_zip_image:

    Gammaeval_sp = n_poly_mod_eval_step2(Gammacur_sp, Gammainc_sp, ctx_sp->mod);
    nmod_mpolyuu_eval_step2(Aeval_sp, Acur_sp, Ainc_sp, ctx_sp);
    nmod_mpolyuu_eval_step2(Beval_sp, Bcur_sp, Binc_sp, ctx_sp);

    if (Aeval_sp->length < 1 || Beval_sp->length < 1 ||
        n_bpoly_bidegree(Aeval_sp) != A->exps[0] ||
        n_bpoly_bidegree(Beval_sp) != B->exps[0])
    {
        goto pick_zip_prime;
    }

    /* the evaluation killed neither lc(A) nor lc(B) */
    FLINT_ASSERT(Gammaeval_sp != 0);

    success = n_bpoly_mod_gcd_brown_smprime(Geval_sp, Abareval_sp, Bbareval_sp,
                                  Aeval_sp, Beval_sp, ctx_sp->mod, Sp);
    if (!success)
        goto pick_zip_prime;

    FLINT_ASSERT(Geval_sp->length > 0);
    GevaldegXY = n_bpoly_bidegree(Geval_sp);

    if (GevaldegXY > GdegboundXY)
    {
        goto pick_zip_prime;        
    }
    else if (GevaldegXY < GdegboundXY)
    {
        GdegboundXY = GevaldegXY;
        if (GdegboundXY == 0)
        {
            fmpz_mpolyu_one(G, ctx);
            fmpz_mpolyu_swap(Abar, A, ctx);
            fmpz_mpolyu_swap(Bbar, B, ctx);
            success = 1;
            goto cleanup;
        }
        goto pick_bma_prime;
    }

    n_bpoly_scalar_mul_nmod(Geval_sp, Gammaeval_sp, ctx_sp->mod);

    success = nmod_zip_mpolyuu_add_point(newZ,
                            which_check == 1 ? Abareval_sp :
                            which_check == 2 ? Bbareval_sp : Geval_sp);
    if (!success)
        goto pick_bma_prime;

    if (newZ->pointcount < zip_evals)
        goto next_zip_image;

    switch (nmod_mpolyu_zip_find_coeffs(newZ, ctx_sp))
    {
        default:
            FLINT_ASSERT(0);
        case nmod_zip_find_coeffs_no_match:
            goto pick_bma_prime;
        case nmod_zip_find_coeffs_non_invertible:
            goto pick_zip_prime;
        case nmod_zip_find_coeffs_good:
            break;
    }

    FLINT_ASSERT(bits == H->bits);
    changed = fmpz_mpolyu_addinterp_zip(H, Hmodulus, newZ, ctx_sp->mod);
    fmpz_mul_ui(Hmodulus, Hmodulus, ctx_sp->mod.n);

    if (changed)
    {
        /* TODO if the coefficients of H are getting to large? */
        goto pick_zip_prime;
    }

    success = fmpz_mpolyu_content_mpoly_threaded_pool(Hcontent, H, ctx, NULL, 0);
    FLINT_ASSERT(Hcontent->bits == bits);
    if (!success)
    {
        success = 0;
        goto cleanup;
    }

    if (which_check == 1)
    {
        fmpz_mpolyu_divexact_mpoly(Abar, H, 1, Hcontent, ctx);
        if (!fmpz_mpolyuu_divides(G, A, Abar, 2, ctx) ||
            !fmpz_mpolyuu_divides(Bbar, B, G, 2, ctx))
        {
            goto pick_zip_prime;
        }
    }
    else if (which_check == 2)
    {
        fmpz_mpolyu_divexact_mpoly(Bbar, H, 1, Hcontent, ctx);
        if (!fmpz_mpolyuu_divides(G, B, Bbar, 2, ctx) ||
            !fmpz_mpolyuu_divides(Abar, A, G, 2, ctx))
        {
            goto pick_zip_prime;
        }
    }
    else
    {
        FLINT_ASSERT(which_check == 0);
        fmpz_mpolyu_divexact_mpoly(G, H, 1, Hcontent, ctx);
        if (!fmpz_mpolyuu_divides(Abar, A, G, 2, ctx) ||
            !fmpz_mpolyuu_divides(Bbar, B, G, 2, ctx))
        {
            goto pick_zip_prime;
        }
    }

    success = 1;

cleanup:

    /* machine precision workspace */
    flint_free(alphas_sp);
    n_bpoly_clear(Aeval_sp);
    n_bpoly_clear(Beval_sp);
    n_bpoly_clear(Geval_sp);
    n_bpoly_clear(Abareval_sp);
    n_bpoly_clear(Bbareval_sp);
    n_poly_clear(Gammacur_sp);
    n_bpoly_clear(Acur_sp);
    n_bpoly_clear(Bcur_sp);
    n_poly_clear(Gammainc_sp);
    n_polyun_clear(Ainc_sp);
    n_polyun_clear(Binc_sp);
    nmod_bma_mpoly_clear(GLambda_sp);
    nmod_bma_mpoly_clear(AbarLambda_sp);
    nmod_bma_mpoly_clear(BbarLambda_sp);
    n_poly_stack_clear(Sp->poly_stack);
    n_bpoly_stack_clear(Sp->bpoly_stack);
    nmod_mpoly_ctx_clear(ctx_sp);

    /* multiprecision workspace */
    _fmpz_vec_clear(alphas_mp, ctx->minfo->nvars);
    fmpz_clear(Gammaeval_mp);
    fmpz_mod_mpolyn_clear(Aeval_mp, ctx_mp);
    fmpz_mod_mpolyn_clear(Beval_mp, ctx_mp);
    fmpz_mod_mpolyn_clear(Geval_mp, ctx_mp);
    fmpz_mod_mpolyn_clear(Abareval_mp, ctx_mp);
    fmpz_mod_mpolyn_clear(Bbareval_mp, ctx_mp);
    fmpz_mpolyc_clear(Gammainc_mp);
    fmpz_mpolyc_clear(Gammacur_mp);
    fmpz_mpolyc_clear(Gammared_mp);
    fmpz_mpolycu_clear(Ainc_mp);
    fmpz_mpolycu_clear(Acur_mp);
    fmpz_mpolycu_clear(Ared_mp);
    fmpz_mpolycu_clear(Binc_mp);
    fmpz_mpolycu_clear(Bcur_mp);
    fmpz_mpolycu_clear(Bred_mp);
    fmpz_mod_bma_mpoly_clear(GLambda_mp, ctx_mp->ffinfo);
    fmpz_mod_bma_mpoly_clear(AbarLambda_mp, ctx_mp->ffinfo);
    fmpz_mod_bma_mpoly_clear(BbarLambda_mp, ctx_mp->ffinfo);
    fmpz_mod_mpoly_ctx_clear(ctx_mp);

    fmpz_mpoly_clear(Hcontent, ctx);
    fmpz_mpolyu_clear(H, ctx);
    flint_free(Adegs);
    flint_free(Bdegs);
    nmod_zip_mpolyu_clear(newZ);

    mpoly_bma_interpolate_ctx_clear(Ictx);

    fmpz_clear(Hmodulus);
    fmpz_clear(last_unlucky_sshift_plus_1_mp);
    fmpz_clear(sshift_mp);
    fmpz_clear(cBksub);
    fmpz_clear(cAksub);
    fmpz_clear(subprod);
    fmpz_clear(image_count_mp);
    fmpz_clear(pm1);
    fmpz_clear(p);
    flint_randclear(randstate);

    FLINT_ASSERT(G->bits == bits);
    FLINT_ASSERT(Abar->bits == bits);
    FLINT_ASSERT(Bbar->bits == bits);

    return success;
}


int fmpz_mpoly_gcd_berlekamp_massey(
    fmpz_mpoly_t G,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    flint_bitcnt_t wbits;
    int success = 0;
    fmpz_mpoly_ctx_t uctx;
    fmpz_mpolyu_t Auu, Buu, Guu, Abaruu, Bbaruu;
    fmpz_mpoly_t Ac, Bc, Gc, Gamma;
    slong * Adegs, * Bdegs, * perm;
    ulong * shift, * stride;
    ulong max_main_degree, max_minor_degree;

    if (fmpz_mpoly_is_zero(A, ctx))
    {
        if (fmpz_mpoly_is_zero(B, ctx))
        {
            fmpz_mpoly_zero(G, ctx);
            return 1;
        }
        if (fmpz_sgn(B->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, B, ctx);
            return 1;
        }
        else
        {
            fmpz_mpoly_set(G, B, ctx);
            return 1;
        }
    }

    if (fmpz_mpoly_is_zero(B, ctx))
    {
        if (fmpz_sgn(A->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, A, ctx);
            return 1;
        }
        else
        {
            fmpz_mpoly_set(G, A, ctx);
            return 1;
        }
    }

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
    {
        return 0;
    }

    if (ctx->minfo->nvars < 3)
    {
        return fmpz_mpoly_gcd_zippel(G, A, B, ctx);
    }

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->nvars >= 3);
    FLINT_ASSERT(!fmpz_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!fmpz_mpoly_is_zero(B, ctx));

    Adegs = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    Bdegs = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    perm = (slong *) flint_malloc((ctx->minfo->nvars)*sizeof(slong));
    shift = (ulong *) flint_malloc((ctx->minfo->nvars)*sizeof(ulong));
    stride = (ulong *) flint_malloc((ctx->minfo->nvars)*sizeof(ulong));

    mpoly_degrees_si(Adegs, A->exps, A->length, A->bits, ctx->minfo);
    mpoly_degrees_si(Bdegs, B->exps, B->length, B->bits, ctx->minfo);

    max_main_degree = 0;
    max_minor_degree = 0;
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        perm[i] = i;
        shift[i] = 0;
        stride[i] = 1;
        FLINT_ASSERT(Adegs[i] >= 0);
        FLINT_ASSERT(Bdegs[i] >= 0);
        if (i < 2)
        {
            max_main_degree = FLINT_MAX(max_main_degree, Adegs[i]);
            max_main_degree = FLINT_MAX(max_main_degree, Bdegs[i]);
        }
        else
        {
            max_minor_degree = FLINT_MAX(max_minor_degree, Adegs[i]);
            max_minor_degree = FLINT_MAX(max_minor_degree, Bdegs[i]);
        }
    }

    fmpz_mpoly_ctx_init(uctx, ctx->minfo->nvars - 2, ORD_LEX);

    /* wbits is bits for intermediates in ZZ[x_0,x_1][x_2,...,x_(n-1)] */
    wbits = 1 + FLINT_BIT_COUNT(max_minor_degree);
    wbits = FLINT_MAX(MPOLY_MIN_BITS, wbits);
    wbits = mpoly_fix_bits(wbits, uctx->minfo);
    FLINT_ASSERT(wbits <= FLINT_BITS);

    fmpz_mpolyu_init(Auu, wbits, uctx);
    fmpz_mpolyu_init(Buu, wbits, uctx);
    fmpz_mpolyu_init(Guu, wbits, uctx);
    fmpz_mpolyu_init(Abaruu, wbits, uctx);
    fmpz_mpolyu_init(Bbaruu, wbits, uctx);
    fmpz_mpoly_init3(Ac, 0, wbits, uctx);
    fmpz_mpoly_init3(Bc, 0, wbits, uctx);
    fmpz_mpoly_init3(Gc, 0, wbits, uctx);
    fmpz_mpoly_init3(Gamma, 0, wbits, uctx);

    /* two main variables must be packed into bits = FLINT_BITS/2 */
    if (FLINT_BIT_COUNT(max_main_degree) >= FLINT_BITS/2)
    {
        success = 0;
        goto cleanup;
    }

    fmpz_mpoly_to_mpolyuu_perm_deflate_threaded_pool(Auu, uctx, A, ctx,
                                           perm, shift, stride, NULL, NULL, 0);
    fmpz_mpoly_to_mpolyuu_perm_deflate_threaded_pool(Buu, uctx, B, ctx,
                                           perm, shift, stride, NULL, NULL, 0);

    success = fmpz_mpolyu_content_mpoly_threaded_pool(Ac, Auu, uctx, NULL, 0);
    success = success &&
              fmpz_mpolyu_content_mpoly_threaded_pool(Bc, Buu, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    fmpz_mpolyu_divexact_mpoly_inplace(Auu, Ac, uctx);
    fmpz_mpolyu_divexact_mpoly_inplace(Buu, Bc, uctx);

    success = _fmpz_mpoly_gcd_threaded_pool(Gamma, wbits, Auu->coeffs + 0,
                                            Buu->coeffs + 0, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    success = fmpz_mpolyuu_gcd_berlekamp_massey(Guu, Abaruu, Bbaruu,
                                                        Auu, Buu, Gamma, uctx);
    if (!success)
        goto cleanup;

    success = _fmpz_mpoly_gcd_threaded_pool(Gc, wbits, Ac, Bc, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    fmpz_mpolyu_mul_mpoly_inplace(Guu, Gc, uctx);

    fmpz_mpoly_from_mpolyuu_perm_inflate(G, FLINT_MIN(A->bits, B->bits), ctx,
                                               Guu, uctx, perm, shift, stride);
    if (fmpz_sgn(G->coeffs + 0) < 0)
        fmpz_mpoly_neg(G, G, ctx);

    success = 1;

cleanup:

    flint_free(Adegs);
    flint_free(Bdegs);
    flint_free(perm);
    flint_free(shift);
    flint_free(stride);

    fmpz_mpolyu_clear(Auu, uctx);
    fmpz_mpolyu_clear(Buu, uctx);
    fmpz_mpolyu_clear(Guu, uctx);
    fmpz_mpolyu_clear(Abaruu, uctx);
    fmpz_mpolyu_clear(Bbaruu, uctx);
    fmpz_mpoly_clear(Ac, uctx);
    fmpz_mpoly_clear(Bc, uctx);
    fmpz_mpoly_clear(Gc, uctx);
    fmpz_mpoly_clear(Gamma, uctx);

    fmpz_mpoly_ctx_clear(uctx);

    return success;
}

