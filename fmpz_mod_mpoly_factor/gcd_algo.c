/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_vec.h"
#include "fmpz_mod_mpoly.h"
#include "fmpz_mod_mpoly_factor.h"

/*
    For each j, set out[j] to the evaluation of A at x_i = alpha[i] (i != j)
    i.e. if nvars = 3
        out[0] = A(x, alpha[1], alpha[2])
        out[1] = A(alpha[0], x, alpha[2])
        out[2] = A(alpha[0], alpha[1], x)

    If ignore[j] is nonzero, then out[j] need not be calculated, probably
    because we shouldn't calculate it in dense form.
*/
static void fmpz_mod_mpoly_evals(
    slong * Atdeg,  /* total degree of deflated A, or -1 for overflow */
    fmpz_mod_poly_struct * out,
    const int * ignore,
    const fmpz_mod_mpoly_t A,
    ulong * Amin_exp,
    ulong * Amax_exp,
    ulong * Astride,
    const fmpz * alphas,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong nvars = ctx->minfo->nvars;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    slong * offsets, * shifts;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    ulong * varexps;
    ulong varexp;
    slong total_degree, lo, hi;
    fmpz_t meval, t, t1;

    fmpz_init(meval);
    fmpz_init(t);
    fmpz_init(t1);

    offsets = FLINT_ARRAY_ALLOC(2*nvars, slong);
    shifts = offsets + nvars;
    varexps = FLINT_ARRAY_ALLOC(nvars, ulong);
    for (j = 0; j < nvars; j++)
    {
        fmpz_mod_poly_zero(out + j, ctx->ffinfo);
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, ctx->minfo);
    }

    total_degree = 0;
    for (i = 0; i < A->length; i++)
    {
        fmpz * s = A->coeffs + i;

        hi = lo = 0;
        for (j = 0; j < nvars; j++)
        {
            varexp = ((A->exps + N*i)[offsets[j]]>>shifts[j])&mask;

            FLINT_ASSERT((Astride[j] == 0 && varexp == Amin_exp[j]) ||
                                     (varexp - Amin_exp[j]) % Astride[j] == 0);

            varexps[j] = Astride[j] < 2 ? varexp - Amin_exp[j] :
                                         (varexp - Amin_exp[j])/Astride[j];

            add_ssaaaa(hi, lo, hi, lo, 0, varexps[j]);

            fmpz_mod_pow_ui(t1, alphas + j, varexps[j], ctx->ffinfo);
            fmpz_mod_mul(meval, s, t1, ctx->ffinfo);
            s = meval;
        }

        if (hi == 0 && FLINT_SIGN_EXT(lo) == 0 && total_degree >= 0)
            total_degree = FLINT_MAX(total_degree, lo);
        else
            total_degree = -1;

        for (j = 0; j < nvars; j++)
        {
            varexp = varexps[j];

            if (ignore[j])
                continue;

            fmpz_mod_poly_fit_length(out + j, varexp + 1, ctx->ffinfo);

            while (out[j].length <= varexp)
            {
                fmpz_zero(out[j].coeffs + out[j].length);
                out[j].length++;
            }

            fmpz_mod_inv(t1, alphas + j, ctx->ffinfo);
            fmpz_mod_pow_ui(t1, t1, varexps[j], ctx->ffinfo);
            fmpz_mod_mul(t, meval, t1, ctx->ffinfo);
            fmpz_mod_add(out[j].coeffs + varexp, out[j].coeffs + varexp, t,
                                                                  ctx->ffinfo);
        }
    }

    *Atdeg = total_degree;

    for (j = 0; j < nvars; j++)
        _fmpz_mod_poly_normalise(out + j);

    flint_free(offsets);
    flint_free(varexps);

    fmpz_clear(meval);
    fmpz_clear(t);
    fmpz_clear(t1);
}


static void _set_estimates(
    mpoly_gcd_info_t I,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int try_count = 0;
    slong nvars = ctx->minfo->nvars;
    slong i, j;
    fmpz_mod_poly_t Geval;
    fmpz_mod_poly_struct * Aevals, * Bevals;
    fmpz * alphas;
    flint_rand_t state;
    slong ignore_limit;
    int * ignore;

    flint_randinit(state);

    ignore = FLINT_ARRAY_ALLOC(nvars, int);
    alphas  = _fmpz_vec_init(nvars);
    Aevals = FLINT_ARRAY_ALLOC(nvars, fmpz_mod_poly_struct);
    Bevals = FLINT_ARRAY_ALLOC(nvars, fmpz_mod_poly_struct);

    fmpz_mod_poly_init(Geval, ctx->ffinfo);
    for (j = 0; j < nvars; j++)
    {
        fmpz_mod_poly_init(Aevals + j, ctx->ffinfo);
        fmpz_mod_poly_init(Bevals + j, ctx->ffinfo);
    }

    ignore_limit = (A->length + B->length)/4096;
    ignore_limit = FLINT_MAX(WORD(9999), ignore_limit);
    I->Gdeflate_deg_bounds_are_nice = 1;
    for (j = 0; j < nvars; j++)
    {
        if (I->Adeflate_deg[j] > ignore_limit ||
            I->Bdeflate_deg[j] > ignore_limit)
        {
            ignore[j] = 1;
            I->Gdeflate_deg_bounds_are_nice = 0;
        }
        else
        {
            ignore[j] = 0;
        }
    }

try_again:

    if (++try_count > 10)
    {
        I->Gdeflate_deg_bounds_are_nice = 0;
        for (j = 0; j < nvars; j++)
        {
            I->Gdeflate_deg_bound[j] = FLINT_MIN(I->Adeflate_deg[j],
                                                 I->Bdeflate_deg[j]);
            I->Gterm_count_est[j] = 1 + I->Gdeflate_deg_bound[j]/2;
        }

        goto cleanup;
    }

    for (j = 0; j < nvars; j++)
        fmpz_mod_rand_not_zero(alphas + j, state, ctx->ffinfo);

    fmpz_mod_mpoly_evals(&I->Adeflate_tdeg, Aevals, ignore, A,
                            I->Amin_exp, I->Amax_exp, I->Gstride, alphas, ctx);
    fmpz_mod_mpoly_evals(&I->Bdeflate_tdeg, Bevals, ignore, B,
                            I->Bmin_exp, I->Bmax_exp, I->Gstride, alphas, ctx);

    for (j = 0; j < nvars; j++)
    {
        if (ignore[j])
        {
            I->Gdeflate_deg_bound[j] = FLINT_MIN(I->Adeflate_deg[j],
                                                 I->Bdeflate_deg[j]);
            I->Gterm_count_est[j] = 1 + I->Gdeflate_deg_bound[j]/2;
        }
        else
        {
            if (I->Adeflate_deg[j] != fmpz_mod_poly_degree(Aevals + j, ctx->ffinfo) ||
                I->Bdeflate_deg[j] != fmpz_mod_poly_degree(Bevals + j, ctx->ffinfo))
            {
                goto try_again;
            }

            fmpz_mod_poly_gcd(Geval, Aevals + j, Bevals + j, ctx->ffinfo);

            I->Gterm_count_est[j] = 0;
            I->Gdeflate_deg_bound[j] = fmpz_mod_poly_degree(Geval, ctx->ffinfo);
            for (i = I->Gdeflate_deg_bound[j]; i >= 0; i--)
                I->Gterm_count_est[j] += (Geval->coeffs[i] != 0);
        }
    }

cleanup:

    fmpz_mod_poly_clear(Geval, ctx->ffinfo);
    for (j = 0; j < nvars; j++)
    {
        fmpz_mod_poly_clear(Aevals + j, ctx->ffinfo);
        fmpz_mod_poly_clear(Bevals + j, ctx->ffinfo);
    }

    flint_free(ignore);
    _fmpz_vec_clear(alphas, nvars);
    flint_free(Aevals);
    flint_free(Bevals);

    flint_randclear(state);

    return;
}


/* (Abar, Bbar) = (A, B) */
static void _parallel_set(
    fmpz_mod_mpoly_t Abar, /* could be NULL */
    fmpz_mod_mpoly_t Bbar, /* could be NULL */
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    if (Abar == B && Bbar == A)
    {
        FLINT_ASSERT(Abar != NULL && Bbar != NULL);
        fmpz_mod_mpoly_set(Abar, B, ctx);
        fmpz_mod_mpoly_set(Bbar, A, ctx);
        fmpz_mod_mpoly_swap(Abar, Bbar, ctx);
    }
    else if (Abar == B && Bbar != A)
    {
        FLINT_ASSERT(Abar != NULL);
        if (Bbar != NULL)
            fmpz_mod_mpoly_set(Bbar, B, ctx);
        fmpz_mod_mpoly_set(Abar, A, ctx);
    }
    else
    {
        if (Abar != NULL)
            fmpz_mod_mpoly_set(Abar, A, ctx);
        if (Bbar != NULL)
            fmpz_mod_mpoly_set(Bbar, B, ctx);
    }
}


/* The variables in ess(A) and ess(B) are disjoint. gcd is trivial to compute */
static void _do_trivial(
    fmpz_mod_mpoly_t G,
    fmpz_mod_mpoly_t Abar,  /* could be NULL */
    fmpz_mod_mpoly_t Bbar,  /* could be NULL */
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const mpoly_gcd_info_t I,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    _parallel_set(Abar, Bbar, A, B, ctx);

    if (Abar != NULL)
        mpoly_monomials_shift_right_ui(Abar->exps, Abar->bits, Abar->length,
                                                      I->Gmin_exp, ctx->minfo);

    if (Bbar != NULL)
        mpoly_monomials_shift_right_ui(Bbar->exps, Bbar->bits, Bbar->length,
                                                      I->Gmin_exp, ctx->minfo);

    fmpz_mod_mpoly_fit_length_reset_bits(G, 1, I->Gbits, ctx);
    mpoly_set_monomial_ui(G->exps, I->Gmin_exp, I->Gbits, ctx->minfo);
    fmpz_one(G->coeffs + 0);
    _fmpz_mod_mpoly_set_length(G, 1, ctx);
}

/*********************** Easy when B is a monomial ***************************/
static int _do_monomial_gcd(
    fmpz_mod_mpoly_t G,
    fmpz_mod_mpoly_t Abar,  /* could be NULL */
    fmpz_mod_mpoly_t Bbar,  /* could be NULL */
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    flint_bitcnt_t Gbits = FLINT_MIN(A->bits, B->bits);
    fmpz * minAfields, * minAdegs, * minBdegs;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length == 1);

    TMP_START;

    /* get the field-wise minimum of A */
    minAfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_init(minAfields + i);
    mpoly_min_fields_fmpz(minAfields, A->exps, A->length, A->bits, ctx->minfo);

    /* unpack to get the min degrees of each variable in A */
    minAdegs = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init(minAdegs + i);
    mpoly_get_monomial_ffmpz_unpacked_ffmpz(minAdegs, minAfields, ctx->minfo);

    /* get the degree of each variable in B */
    minBdegs = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init(minBdegs + i);
    mpoly_get_monomial_ffmpz(minBdegs, B->exps, B->bits, ctx->minfo);

    /* compute the degree of each variable in G */
    _fmpz_vec_min_inplace(minBdegs, minAdegs, ctx->minfo->nvars);

    _parallel_set(Abar, Bbar, A, B, ctx);

    if (Abar != NULL)
        mpoly_monomials_shift_right_ffmpz(Abar->exps, Abar->bits, Abar->length,
                                                         minBdegs, ctx->minfo);

    if (Bbar != NULL)
        mpoly_monomials_shift_right_ffmpz(Bbar->exps, Bbar->bits, Bbar->length,
                                                         minBdegs, ctx->minfo);

    fmpz_mod_mpoly_fit_length_reset_bits(G, 1, Gbits, ctx);
    mpoly_set_monomial_ffmpz(G->exps, minBdegs, Gbits, ctx->minfo);
    fmpz_one(G->coeffs + 0);
    _fmpz_mod_mpoly_set_length(G, 1, ctx);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(minAfields + i);
    }
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        fmpz_clear(minAdegs + i);
        fmpz_clear(minBdegs + i);
    }

    TMP_END;

    return 1;
}


/********************** See if cofactors are monomials ***********************/
static int _try_monomial_cofactors(
    fmpz_mod_mpoly_t G,
    fmpz_mod_mpoly_t Abar,  /* could be NULL */
    fmpz_mod_mpoly_t Bbar,  /* could be NULL */
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    slong NA, NG;
    slong nvars = ctx->minfo->nvars;
    fmpz * Abarexps, * Bbarexps, * Texps;
    fmpz_t a0, b0, t1, t2;
    fmpz_mod_mpoly_t T;
    flint_bitcnt_t Gbits = FLINT_MIN(A->bits, B->bits);
    flint_bitcnt_t Abarbits = A->bits;
    flint_bitcnt_t Bbarbits = B->bits;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (A->length != B->length)
        return 0;

    fmpz_init(t1);
    fmpz_init(t2);
    fmpz_init(a0);
    fmpz_init(b0);

    for (i = A->length - 1; i > 0; i--)
    {
        fmpz_mod_mul(t1, A->coeffs + 0, B->coeffs + i, ctx->ffinfo);
        fmpz_mod_mul(t2, B->coeffs + 0, A->coeffs + i, ctx->ffinfo);
        success = fmpz_equal(t1, t2);
        if (!success)
            goto cleanup_less;
    }

    fmpz_set(a0, A->coeffs + 0);
    fmpz_set(b0, B->coeffs + 0);

    TMP_START;

    Abarexps = (fmpz *) TMP_ALLOC(3*nvars*sizeof(fmpz));
    Bbarexps = Abarexps + 1*nvars;
    Texps    = Abarexps + 2*nvars;
    for (j = 0; j < nvars; j++)
    {
        fmpz_init(Abarexps + j);
        fmpz_init(Bbarexps + j);
        fmpz_init(Texps + j);
    }

    success = mpoly_monomial_cofactors(Abarexps, Bbarexps, A->exps, A->bits,
                                      B->exps, B->bits, A->length, ctx->minfo);
    if (!success)
        goto cleanup_more;

    fmpz_mod_mpoly_init3(T, A->length, Gbits, ctx);
    NG = mpoly_words_per_exp(Gbits, ctx->minfo);
    NA = mpoly_words_per_exp(A->bits, ctx->minfo);
    fmpz_mod_inv(t1, A->coeffs + 0, ctx->ffinfo);
    T->length = A->length;
    for (i = 0; i < A->length; i++)
    {
        mpoly_get_monomial_ffmpz(Texps, A->exps + NA*i, A->bits, ctx->minfo);
        _fmpz_vec_sub(Texps, Texps, Abarexps, nvars);
        mpoly_set_monomial_ffmpz(T->exps + NG*i, Texps, Gbits, ctx->minfo);
        fmpz_mod_mul(T->coeffs + i, A->coeffs + i, t1, ctx->ffinfo);
    }
    fmpz_mod_mpoly_swap(G, T, ctx);
    fmpz_mod_mpoly_clear(T, ctx);

    if (Abar != NULL)
    {
        fmpz_mod_mpoly_fit_length_reset_bits(Abar, 1, Abarbits, ctx);
        mpoly_set_monomial_ffmpz(Abar->exps, Abarexps, Abarbits, ctx->minfo);
        fmpz_set(Abar->coeffs + 0, a0);
        _fmpz_mod_mpoly_set_length(Abar, 1, ctx);
    }

    if (Bbar != NULL)
    {
        fmpz_mod_mpoly_fit_length_reset_bits(Bbar, 1, Bbarbits, ctx);
        mpoly_set_monomial_ffmpz(Bbar->exps, Bbarexps, Bbarbits, ctx->minfo);
        fmpz_set(Bbar->coeffs + 0, b0);
        _fmpz_mod_mpoly_set_length(Bbar, 1, ctx);
    }

    success = 1;

cleanup_more:

    for (j = 0; j < nvars; j++)
    {
        fmpz_clear(Abarexps + j);
        fmpz_clear(Bbarexps + j);
        fmpz_clear(Texps + j);
    }

    TMP_END;

cleanup_less:

    fmpz_clear(t1);
    fmpz_clear(t2);
    fmpz_clear(a0);
    fmpz_clear(b0);

    return success;
}


/***  ess(A) and ess(B) depend on only one variable v_in_both ****************/
static int _do_univar(
    fmpz_mod_mpoly_t G,
    fmpz_mod_mpoly_t Abar,
    fmpz_mod_mpoly_t Bbar,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    slong v_in_both,
    const mpoly_gcd_info_t I,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_poly_t a, b, g, t, r;

    fmpz_mod_poly_init(a, ctx->ffinfo);
    fmpz_mod_poly_init(b, ctx->ffinfo);
    fmpz_mod_poly_init(g, ctx->ffinfo);
    fmpz_mod_poly_init(t, ctx->ffinfo);
    fmpz_mod_poly_init(r, ctx->ffinfo);

    _fmpz_mod_mpoly_to_fmpz_mod_poly_deflate(a, A, v_in_both, I->Amin_exp, I->Gstride, ctx);
    _fmpz_mod_mpoly_to_fmpz_mod_poly_deflate(b, B, v_in_both, I->Bmin_exp, I->Gstride, ctx);
    fmpz_mod_poly_gcd(g, a, b, ctx->ffinfo);
    _fmpz_mod_mpoly_from_fmpz_mod_poly_inflate(G, I->Gbits, g, v_in_both,
                                                 I->Gmin_exp, I->Gstride, ctx);
    if (Abar != NULL)
    {
        fmpz_mod_poly_divrem(t, r, a, g, ctx->ffinfo);
        _fmpz_mod_mpoly_from_fmpz_mod_poly_inflate(Abar, I->Abarbits, t, v_in_both,
                                              I->Abarmin_exp, I->Gstride, ctx);
    }

    if (Bbar != NULL)
    {
        fmpz_mod_poly_divrem(t, r, b, g, ctx->ffinfo);
        _fmpz_mod_mpoly_from_fmpz_mod_poly_inflate(Bbar, I->Bbarbits, t, v_in_both,
                                              I->Bbarmin_exp, I->Gstride, ctx);
    }

    fmpz_mod_poly_clear(a, ctx->ffinfo);
    fmpz_mod_poly_clear(b, ctx->ffinfo);
    fmpz_mod_poly_clear(g, ctx->ffinfo);
    fmpz_mod_poly_clear(t, ctx->ffinfo);
    fmpz_mod_poly_clear(r, ctx->ffinfo);

    return 1;
}


/********* Assume B has length one when converted to univar format ***********/
static int _try_missing_var(
    fmpz_mod_mpoly_t G, flint_bitcnt_t Gbits,
    fmpz_mod_mpoly_t Abar,
    fmpz_mod_mpoly_t Bbar,
    slong var,
    const fmpz_mod_mpoly_t A, ulong Ashift,
    const fmpz_mod_mpoly_t B, ulong Bshift,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;
    fmpz_mod_mpoly_univar_t Au;

    fmpz_mod_mpoly_univar_init(Au, ctx);

#if FLINT_WANT_ASSERT
    fmpz_mod_mpoly_to_univar(Au, B, var, ctx);
    FLINT_ASSERT(Au->length == 1);
#endif

    fmpz_mod_mpoly_to_univar(Au, A, var, ctx);

    fmpz_mod_mpoly_univar_fit_length(Au, Au->length + 1, ctx);
    fmpz_mod_mpoly_set(Au->coeffs + Au->length, B, ctx);
    Au->length++;

    if (Abar == NULL && Bbar == NULL)
    {
        success = _fmpz_mod_mpoly_vec_content_mpoly(G, Au->coeffs, Au->length, ctx);
        if (!success)
            goto cleanup;

        fmpz_mod_mpoly_repack_bits_inplace(G, Gbits, ctx);
        _mpoly_gen_shift_left(G->exps, G->bits, G->length,
                                   var, FLINT_MIN(Ashift, Bshift), ctx->minfo);
    }
    else
    {
        fmpz_mod_mpoly_t tG, tAbar, tBbar;

        fmpz_mod_mpoly_init(tG, ctx);
        fmpz_mod_mpoly_init(tAbar, ctx);
        fmpz_mod_mpoly_init(tBbar, ctx);

        success = _fmpz_mod_mpoly_vec_content_mpoly(tG, Au->coeffs, Au->length, ctx);
        if (!success)
            goto cleanup;

        fmpz_mod_mpoly_repack_bits_inplace(tG, Gbits, ctx);
        _mpoly_gen_shift_left(tG->exps, tG->bits, tG->length,
                                   var, FLINT_MIN(Ashift, Bshift), ctx->minfo);

        if (Abar != NULL)
        {
            success = fmpz_mod_mpoly_divides(tAbar, A, tG, ctx);
            FLINT_ASSERT(success);
        }

        if (Bbar != NULL)
        {
            success = fmpz_mod_mpoly_divides(tBbar, B, tG, ctx);
            FLINT_ASSERT(success);
        }

        fmpz_mod_mpoly_swap(G, tG, ctx);

        if (Abar != NULL)
            fmpz_mod_mpoly_swap(Abar, tAbar, ctx);

        if (Bbar != NULL)
            fmpz_mod_mpoly_swap(Bbar, tBbar, ctx);

        fmpz_mod_mpoly_clear(tG, ctx);
        fmpz_mod_mpoly_clear(tAbar, ctx);
        fmpz_mod_mpoly_clear(tBbar, ctx);        
    }

    success = 1;

cleanup:

    fmpz_mod_mpoly_univar_clear(Au, ctx);

    return success;
}


/************************ See if B divides A ********************************/
static int _try_divides(
    fmpz_mod_mpoly_t G,
    fmpz_mod_mpoly_t Abar,
    fmpz_mod_mpoly_t Bbar,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t BB,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success = 0;
    fmpz_mod_mpoly_t Q, B, M;

    fmpz_mod_mpoly_init(Q, ctx);
    fmpz_mod_mpoly_init(B, ctx);
    fmpz_mod_mpoly_init(M, ctx);

    /* BB = M*B */
    fmpz_mod_mpoly_term_content(M, BB, ctx);
    fmpz_mod_mpoly_divides(B, BB, M, ctx);

    success = fmpz_mod_mpoly_divides(Q, A, B, ctx);
    if (success)
    {
        _do_monomial_gcd(G, Abar, Bbar, Q, M, ctx);
        fmpz_mod_mpoly_mul(G, G, B, ctx);
    }

    fmpz_mod_mpoly_clear(Q, ctx);
    fmpz_mod_mpoly_clear(B, ctx);
    fmpz_mod_mpoly_clear(M, ctx);

    return success;
}


/********************* Hit A and B with a prs ********************************/

static int _try_prs(
    fmpz_mod_mpoly_t G,
    fmpz_mod_mpoly_t Abar,
    fmpz_mod_mpoly_t Bbar,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const mpoly_gcd_info_t I,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;
    slong j, var = -WORD(1);
    fmpz_mod_mpoly_t Ac, Bc, Gc, s, t;
    fmpz_mod_mpoly_univar_t Ax, Bx, Gx;

    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        if (I->Amax_exp[j] <= I->Amin_exp[j] ||
            I->Bmax_exp[j] <= I->Bmin_exp[j] /* not needed */)
        {
            FLINT_ASSERT(I->Amax_exp[j] == I->Amin_exp[j]);
            FLINT_ASSERT(I->Bmax_exp[j] == I->Bmin_exp[j]);
            continue;
        }

        FLINT_ASSERT(I->Gstride[j] != UWORD(0));
        FLINT_ASSERT((I->Amax_exp[j] - I->Amin_exp[j]) % I->Gstride[j] == 0);
        FLINT_ASSERT((I->Bmax_exp[j] - I->Bmin_exp[j]) % I->Gstride[j] == 0);

        var = j;
        break;
    }

    if (var < 0)
        return 0;

    fmpz_mod_mpoly_init(Ac, ctx);
    fmpz_mod_mpoly_init(Bc, ctx);
    fmpz_mod_mpoly_init(Gc, ctx);
    fmpz_mod_mpoly_init(s, ctx);
    fmpz_mod_mpoly_init(t, ctx);
    fmpz_mod_mpoly_univar_init(Ax, ctx);
    fmpz_mod_mpoly_univar_init(Bx, ctx);
    fmpz_mod_mpoly_univar_init(Gx, ctx);

    fmpz_mod_mpoly_to_univar(Ax, A, var, ctx);
    fmpz_mod_mpoly_to_univar(Bx, B, var, ctx);

    success = _fmpz_mod_mpoly_vec_content_mpoly(Ac, Ax->coeffs, Ax->length, ctx) &&
              _fmpz_mod_mpoly_vec_content_mpoly(Bc, Bx->coeffs, Bx->length, ctx) &&
              fmpz_mod_mpoly_gcd(Gc, Ac, Bc, ctx);
    if (!success)
        goto cleanup;

    _fmpz_mod_mpoly_vec_divexact_mpoly(Ax->coeffs, Ax->length, Ac, ctx);
    _fmpz_mod_mpoly_vec_divexact_mpoly(Bx->coeffs, Bx->length, Bc, ctx);

    success = fmpz_cmp(Ax->exps + 0, Bx->exps + 0) > 0 ?
                            fmpz_mod_mpoly_univar_pseudo_gcd(Gx, Ax, Bx, ctx) :
                            fmpz_mod_mpoly_univar_pseudo_gcd(Gx, Bx, Ax, ctx);
    if (!success)
        goto cleanup;

    if (fmpz_mod_mpoly_gcd(t, Ax->coeffs + 0, Bx->coeffs + 0, ctx) &&
                                                                t->length == 1)
    {
        fmpz_mod_mpoly_term_content(s, Gx->coeffs + 0, ctx);
        fmpz_mod_mpoly_divexact(t, Gx->coeffs + 0, s, ctx);
        _fmpz_mod_mpoly_vec_divexact_mpoly(Gx->coeffs, Gx->length, t, ctx);
    }
    else if (fmpz_mod_mpoly_gcd(t, Ax->coeffs + Ax->length - 1,
                           Bx->coeffs + Bx->length - 1, ctx) && t->length == 1)
    {
        fmpz_mod_mpoly_term_content(s, Gx->coeffs + Gx->length - 1, ctx);
        fmpz_mod_mpoly_divexact(t, Gx->coeffs + Gx->length - 1, s, ctx);
        _fmpz_mod_mpoly_vec_divexact_mpoly(Gx->coeffs, Gx->length, t, ctx);
    }

    success = _fmpz_mod_mpoly_vec_content_mpoly(t, Gx->coeffs, Gx->length, ctx);
    if (!success)
        goto cleanup;

    _fmpz_mod_mpoly_vec_divexact_mpoly(Gx->coeffs, Gx->length, t, ctx);
    _fmpz_mod_mpoly_vec_mul_mpoly(Gx->coeffs, Gx->length, Gc, ctx);
    _fmpz_mod_mpoly_from_univar(Gc, I->Gbits, Gx, var, ctx);

    if (Abar != NULL)
        fmpz_mod_mpoly_divexact(s, A, Gc, ctx);

    if (Bbar != NULL)
        fmpz_mod_mpoly_divexact(t, B, Gc, ctx);

    fmpz_mod_mpoly_swap(G, Gc, ctx);

    if (Abar != NULL)
        fmpz_mod_mpoly_swap(Abar, s, ctx);

    if (Bbar != NULL)
        fmpz_mod_mpoly_swap(Bbar, t, ctx);

    success = 1;

cleanup:

    fmpz_mod_mpoly_clear(Ac, ctx);
    fmpz_mod_mpoly_clear(Bc, ctx);
    fmpz_mod_mpoly_clear(Gc, ctx);
    fmpz_mod_mpoly_clear(s, ctx);
    fmpz_mod_mpoly_clear(t, ctx);
    fmpz_mod_mpoly_univar_clear(Ax, ctx);
    fmpz_mod_mpoly_univar_clear(Bx, ctx);
    fmpz_mod_mpoly_univar_clear(Gx, ctx);

    return success;
}


/********************** Hit A and B with zippel ******************************/

static int _try_zippel(
    fmpz_mod_mpoly_t G,
    fmpz_mod_mpoly_t Abar,
    fmpz_mod_mpoly_t Bbar,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const mpoly_gcd_info_t I,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, k;
    slong m = I->mvars;
    int success;
    flint_bitcnt_t wbits;
    fmpz_mod_mpoly_ctx_t lctx;
    fmpz_mod_mpoly_t Al, Bl, Gl, Abarl, Bbarl;
    fmpz_mod_mpoly_t Ac, Bc, Gc, Abarc, Bbarc;
    slong max_deg;
    flint_rand_t state;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (!(I->can_use & MPOLY_GCD_USE_ZIPPEL))
        return 0;

    FLINT_ASSERT(m >= 2);

    flint_randinit(state);

    fmpz_mod_mpoly_ctx_init(lctx, m, ORD_LEX, fmpz_mod_ctx_modulus(ctx->ffinfo));

    max_deg = 0;
    for (i = 0; i < m; i++)
    {
        k = I->zippel_perm[i];
        max_deg = FLINT_MAX(max_deg, I->Adeflate_deg[k]);
        max_deg = FLINT_MAX(max_deg, I->Bdeflate_deg[k]);
    }

    wbits = 1 + FLINT_BIT_COUNT(max_deg);
    wbits = mpoly_fix_bits(wbits, lctx->minfo);
    FLINT_ASSERT(wbits <= FLINT_BITS);

    fmpz_mod_mpoly_init3(Al, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Bl, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Gl, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Abarl, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Bbarl, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Ac, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Bc, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Gc, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Abarc, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Bbarc, 0, wbits, lctx);

    fmpz_mod_mpoly_to_mpolyl_perm_deflate(Al, lctx, A, ctx,
                                      I->zippel_perm, I->Amin_exp, I->Gstride);
    fmpz_mod_mpoly_to_mpolyl_perm_deflate(Bl, lctx, B, ctx,
                                      I->zippel_perm, I->Bmin_exp, I->Gstride);

    success = fmpz_mod_mpolyl_content(Ac, Al, 1, lctx) &&
              fmpz_mod_mpolyl_content(Bc, Bl, 1, lctx);
    if (!success)
        goto cleanup;

    success = _fmpz_mod_mpoly_gcd_algo(Gc, Abar == NULL ? NULL : Abarc,
                 Bbar == NULL ? NULL : Bbarc, Ac, Bc, lctx, MPOLY_GCD_USE_ALL);
    if (!success)
        goto cleanup;

    success = fmpz_mod_mpoly_divides(Al, Al, Ac, lctx);
    FLINT_ASSERT(success);

    success = fmpz_mod_mpoly_divides(Bl, Bl, Bc, lctx);
    FLINT_ASSERT(success);

    fmpz_mod_mpoly_repack_bits_inplace(Al, wbits, lctx);
    fmpz_mod_mpoly_repack_bits_inplace(Bl, wbits, lctx);

    success = fmpz_mod_mpolyl_gcdp_zippel(Gl, Abarl, Bbarl, Al, Bl,
                                                           m - 1, lctx, state);
    if (!success)
        goto cleanup;

    fmpz_mod_mpoly_mul(Gl, Gl, Gc, lctx);
    fmpz_mod_mpoly_from_mpolyl_perm_inflate(G, I->Gbits, ctx, Gl, lctx,
                                      I->zippel_perm, I->Gmin_exp, I->Gstride);
    if (Abar != NULL)
    {
        fmpz_mod_mpoly_mul(Abarl, Abarl, Abarc, lctx);
        fmpz_mod_mpoly_from_mpolyl_perm_inflate(Abar, I->Abarbits, ctx, Abarl, lctx,
                                   I->zippel_perm, I->Abarmin_exp, I->Gstride);
    }

    if (Bbar != NULL)
    {
        fmpz_mod_mpoly_mul(Bbarl, Bbarl, Bbarc, lctx);
        fmpz_mod_mpoly_from_mpolyl_perm_inflate(Bbar, I->Bbarbits, ctx, Bbarl, lctx,
                                   I->zippel_perm, I->Bbarmin_exp, I->Gstride);
    }

    success = 1;

cleanup:

    fmpz_mod_mpoly_clear(Al, lctx);
    fmpz_mod_mpoly_clear(Bl, lctx);
    fmpz_mod_mpoly_clear(Gl, lctx);
    fmpz_mod_mpoly_clear(Abarl, lctx);
    fmpz_mod_mpoly_clear(Bbarl, lctx);
    fmpz_mod_mpoly_clear(Ac, lctx);
    fmpz_mod_mpoly_clear(Bc, lctx);
    fmpz_mod_mpoly_clear(Gc, lctx);
    fmpz_mod_mpoly_clear(Abarc, lctx);
    fmpz_mod_mpoly_clear(Bbarc, lctx);

    fmpz_mod_mpoly_ctx_clear(lctx);
    flint_randclear(state);

    return success;
}

static int _try_zippel2(
    fmpz_mod_mpoly_t G,
    fmpz_mod_mpoly_t Abar,
    fmpz_mod_mpoly_t Bbar,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const mpoly_gcd_info_t I,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, k;
    slong m = I->mvars;
    int success;
    flint_bitcnt_t wbits;
    fmpz_mod_mpoly_ctx_t lctx;
    fmpz_mod_mpoly_t Al, Bl, Gl, Abarl, Bbarl;
    fmpz_mod_mpoly_t Al_lc, Bl_lc, Ac, Bc, Gc, Abarc, Bbarc, Gamma;
    slong * tmp, * Gl_degs, * Al_degs, * Bl_degs, * Gamma_degs, * Gguess;
    slong max_degree;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (!(I->can_use & MPOLY_GCD_USE_ZIPPEL2))
        return 0;

    FLINT_ASSERT(m >= 3);

    tmp = FLINT_ARRAY_ALLOC(5*m, slong);
    Al_degs   = tmp + 1*m;
    Bl_degs   = tmp + 2*m;
    Gl_degs   = tmp + 3*m;
    Gamma_degs = tmp + 4*m;

    fmpz_mod_mpoly_ctx_init(lctx, m, ORD_LEX, fmpz_mod_ctx_modulus(ctx->ffinfo));

    max_degree = 0;
    for (i = 0; i < m; i++)
    {
        k = I->zippel2_perm[i];

        Gl_degs[i] = I->Gdeflate_deg_bound[k];

        Al_degs[i] = I->Adeflate_deg[k];
        max_degree = FLINT_MAX(max_degree, Al_degs[i]);

        Bl_degs[i] = I->Bdeflate_deg[k];
        max_degree = FLINT_MAX(max_degree, Bl_degs[i]);
    }

    wbits = 1 + FLINT_BIT_COUNT(max_degree);
    wbits = mpoly_fix_bits(wbits, lctx->minfo);
    FLINT_ASSERT(wbits <= FLINT_BITS);

    fmpz_mod_mpoly_init3(Al, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Bl, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Gl, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Abarl, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Bbarl, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Ac, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Bc, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Gc, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Abarc, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Bbarc, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Gamma, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Al_lc, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Bl_lc, 0, wbits, lctx);

    fmpz_mod_mpoly_to_mpolyl_perm_deflate(Al, lctx, A, ctx,
                                     I->zippel2_perm, I->Amin_exp, I->Gstride);
    fmpz_mod_mpoly_to_mpolyl_perm_deflate(Bl, lctx, B, ctx,
                                     I->zippel2_perm, I->Bmin_exp, I->Gstride);

    success = fmpz_mod_mpolyl_content(Ac, Al, 2, lctx) &&
              fmpz_mod_mpolyl_content(Bc, Bl, 2, lctx);
    if (!success)
        goto cleanup;

    success = _fmpz_mod_mpoly_gcd_algo(Gc, Abar == NULL ? NULL : Abarc,
                 Bbar == NULL ? NULL : Bbarc, Ac, Bc, lctx, MPOLY_GCD_USE_ALL);
    if (!success)
        goto cleanup;

    fmpz_mod_mpoly_degrees_si(tmp, Ac, lctx);
    for (i = 0; i < m; i++)
        Al_degs[i] -= tmp[i];

    success = fmpz_mod_mpoly_divides(Al, Al, Ac, lctx);
    FLINT_ASSERT(success);

    fmpz_mod_mpoly_degrees_si(tmp, Bc, lctx);
    for (i = 0; i < m; i++)
        Bl_degs[i] -= tmp[i];

    success = fmpz_mod_mpoly_divides(Bl, Bl, Bc, lctx);
    FLINT_ASSERT(success);

    fmpz_mod_mpoly_degrees_si(tmp, Gc, lctx);
    for (i = 0; i < m; i++)
        Gl_degs[i] -= tmp[i];

    fmpz_mod_mpoly_repack_bits_inplace(Al, wbits, lctx);
    fmpz_mod_mpoly_repack_bits_inplace(Bl, wbits, lctx);
    fmpz_mod_mpolyl_lead_coeff(Al_lc, Al, 2, lctx);
    fmpz_mod_mpolyl_lead_coeff(Bl_lc, Bl, 2, lctx);
    success = fmpz_mod_mpoly_gcd(Gamma, Al_lc, Bl_lc, lctx);
    if (!success)
        goto cleanup;
    fmpz_mod_mpoly_repack_bits_inplace(Gamma, wbits, lctx);

    fmpz_mod_mpoly_degrees_si(Gamma_degs, Gamma, lctx);

    Gguess = I->Gdeflate_deg_bounds_are_nice ? Gl_degs : NULL;

    success = fmpz_mod_mpolyl_gcd_zippel2_smprime(Gl, Gguess, Abarl, Bbarl,
                            Al, Al_degs, Bl, Bl_degs, Gamma, Gamma_degs, lctx);
    if (!success)
        goto cleanup;

    fmpz_mod_mpoly_mul(Gl, Gl, Gc, lctx);
    fmpz_mod_mpoly_from_mpolyl_perm_inflate(G, I->Gbits, ctx, Gl, lctx,
                                     I->zippel2_perm, I->Gmin_exp, I->Gstride);
    if (Abar != NULL)
    {
        fmpz_mod_mpoly_mul(Abarl, Abarl, Abarc, lctx);
        fmpz_mod_mpoly_from_mpolyl_perm_inflate(Abar, I->Abarbits, ctx, Abarl, lctx,
                                  I->zippel2_perm, I->Abarmin_exp, I->Gstride);
    }

    if (Bbar != NULL)
    {
        fmpz_mod_mpoly_mul(Bbarl, Bbarl, Bbarc, lctx);
        fmpz_mod_mpoly_from_mpolyl_perm_inflate(Bbar, I->Bbarbits, ctx, Bbarl, lctx,
                                  I->zippel2_perm, I->Bbarmin_exp, I->Gstride);
    }

    success = 1;

cleanup:

    fmpz_mod_mpoly_clear(Al, lctx);
    fmpz_mod_mpoly_clear(Bl, lctx);
    fmpz_mod_mpoly_clear(Gl, lctx);
    fmpz_mod_mpoly_clear(Abarl, lctx);
    fmpz_mod_mpoly_clear(Bbarl, lctx);
    fmpz_mod_mpoly_clear(Ac, lctx);
    fmpz_mod_mpoly_clear(Bc, lctx);
    fmpz_mod_mpoly_clear(Gc, lctx);
    fmpz_mod_mpoly_clear(Abarc, lctx);
    fmpz_mod_mpoly_clear(Bbarc, lctx);
    fmpz_mod_mpoly_clear(Gamma, lctx);
    fmpz_mod_mpoly_clear(Al_lc, lctx);
    fmpz_mod_mpoly_clear(Bl_lc, lctx);

    fmpz_mod_mpoly_ctx_clear(lctx);

    flint_free(tmp);

    return success;
}


/******************** Hit A and B with hensel lifting ************************/
static int _try_hensel(
    fmpz_mod_mpoly_t G,
    fmpz_mod_mpoly_t Abar,
    fmpz_mod_mpoly_t Bbar,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const mpoly_gcd_info_t I,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, k;
    slong m = I->mvars;
    int success;
    flint_bitcnt_t wbits;
    fmpz_mod_mpoly_ctx_t lctx;
    fmpz_mod_mpoly_t Al, Bl, Gl, Abarl, Bbarl;
    fmpz_mod_mpoly_t Ac, Bc, Gc, Abarc, Bbarc;
    slong max_deg;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (!(I->can_use & MPOLY_GCD_USE_HENSEL))
        return 0;

    FLINT_ASSERT(m >= 2);

    fmpz_mod_mpoly_ctx_init(lctx, m, ORD_LEX, fmpz_mod_ctx_modulus(ctx->ffinfo));

    max_deg = 0;
    for (i = 0; i < m; i++)
    {
        k = I->hensel_perm[i];
        max_deg = FLINT_MAX(max_deg, I->Adeflate_deg[k]);
        max_deg = FLINT_MAX(max_deg, I->Bdeflate_deg[k]);
    }

    wbits = 1 + FLINT_BIT_COUNT(max_deg);
    wbits = mpoly_fix_bits(wbits, lctx->minfo);
    FLINT_ASSERT(wbits <= FLINT_BITS);

    fmpz_mod_mpoly_init3(Al, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Bl, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Gl, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Abarl, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Bbarl, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Ac, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Bc, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Gc, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Abarc, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Bbarc, 0, wbits, lctx);

    fmpz_mod_mpoly_to_mpolyl_perm_deflate(Al, lctx, A, ctx,
                                      I->hensel_perm, I->Amin_exp, I->Gstride);
    fmpz_mod_mpoly_to_mpolyl_perm_deflate(Bl, lctx, B, ctx,
                                      I->hensel_perm, I->Bmin_exp, I->Gstride);

    success = fmpz_mod_mpolyl_content(Ac, Al, 1, lctx) &&
              fmpz_mod_mpolyl_content(Bc, Bl, 1, lctx);
    if (!success)
        goto cleanup;

    success = _fmpz_mod_mpoly_gcd_algo(Gc, Abar == NULL ? NULL : Abarc,
                 Bbar == NULL ? NULL : Bbarc, Ac, Bc, lctx, MPOLY_GCD_USE_ALL);
    if (!success)
        goto cleanup;

    success = fmpz_mod_mpoly_divides(Al, Al, Ac, lctx);
    FLINT_ASSERT(success);

    success = fmpz_mod_mpoly_divides(Bl, Bl, Bc, lctx);
    FLINT_ASSERT(success);

    fmpz_mod_mpoly_repack_bits_inplace(Al, wbits, lctx);
    fmpz_mod_mpoly_repack_bits_inplace(Bl, wbits, lctx);

    max_deg = I->Gdeflate_deg_bound[I->hensel_perm[0]];
    success = fmpz_mod_mpolyl_gcd_hensel_smprime(Gl, max_deg, Abarl, Bbarl, Al, Bl, lctx);
    if (!success)
        goto cleanup;

    fmpz_mod_mpoly_mul(Gl, Gl, Gc, lctx);
    fmpz_mod_mpoly_from_mpolyl_perm_inflate(G, I->Gbits, ctx, Gl, lctx,
                                      I->hensel_perm, I->Gmin_exp, I->Gstride);
    if (Abar != NULL)
    {
        fmpz_mod_mpoly_mul(Abarl, Abarl, Abarc, lctx);
        fmpz_mod_mpoly_from_mpolyl_perm_inflate(Abar, I->Abarbits, ctx, Abarl, lctx,
                                   I->hensel_perm, I->Abarmin_exp, I->Gstride);
    }

    if (Bbar != NULL)
    {
        fmpz_mod_mpoly_mul(Bbarl, Bbarl, Bbarc, lctx);
        fmpz_mod_mpoly_from_mpolyl_perm_inflate(Bbar, I->Bbarbits, ctx, Bbarl, lctx,
                                   I->hensel_perm, I->Bbarmin_exp, I->Gstride);
    }

    success = 1;

cleanup:

    fmpz_mod_mpoly_clear(Al, lctx);
    fmpz_mod_mpoly_clear(Bl, lctx);
    fmpz_mod_mpoly_clear(Gl, lctx);
    fmpz_mod_mpoly_clear(Abarl, lctx);
    fmpz_mod_mpoly_clear(Bbarl, lctx);
    fmpz_mod_mpoly_clear(Ac, lctx);
    fmpz_mod_mpoly_clear(Bc, lctx);
    fmpz_mod_mpoly_clear(Gc, lctx);
    fmpz_mod_mpoly_clear(Abarc, lctx);
    fmpz_mod_mpoly_clear(Bbarc, lctx);

    fmpz_mod_mpoly_ctx_clear(lctx);

    return success;
}


/*********************** Hit A and B with brown ******************************/

static int _try_brown(
    fmpz_mod_mpoly_t G,
    fmpz_mod_mpoly_t Abar,
    fmpz_mod_mpoly_t Bbar,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    mpoly_gcd_info_t I,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;
    slong m = I->mvars;
    flint_bitcnt_t wbits;
    fmpz_mod_mpoly_ctx_t nctx;
    fmpz_mod_mpolyn_t An, Bn, Gn, Abarn, Bbarn;
    fmpz_mod_poly_polyun_mpolyn_stack_t St;

    if (!(I->can_use & MPOLY_GCD_USE_BROWN))
        return 0;

    FLINT_ASSERT(m >= 2);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    wbits = FLINT_MAX(A->bits, B->bits);

    fmpz_mod_mpoly_ctx_init(nctx, m, ORD_LEX, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_mpolyn_init(An, wbits, nctx);
    fmpz_mod_mpolyn_init(Bn, wbits, nctx);
    fmpz_mod_mpolyn_init(Gn, wbits, nctx);
    fmpz_mod_mpolyn_init(Abarn, wbits, nctx);
    fmpz_mod_mpolyn_init(Bbarn, wbits, nctx);
    fmpz_mod_poly_stack_init(St->poly_stack);
    fmpz_mod_polyun_stack_init(St->polyun_stack);
    fmpz_mod_mpolyn_stack_init(St->mpolyn_stack, wbits, nctx);

    fmpz_mod_mpoly_to_mpolyn_perm_deflate(An, nctx, A, ctx,
                                       I->brown_perm, I->Amin_exp, I->Gstride);
    fmpz_mod_mpoly_to_mpolyn_perm_deflate(Bn, nctx, B, ctx,
                                       I->brown_perm, I->Bmin_exp, I->Gstride);

    FLINT_ASSERT(An->bits == wbits);
    FLINT_ASSERT(Bn->bits == wbits);
    FLINT_ASSERT(An->length > 0);
    FLINT_ASSERT(Bn->length > 0);

    success = fmpz_mod_mpolyn_gcd_brown_smprime(Gn, Abarn, Bbarn, An, Bn,
                                                           m - 1, nctx, I, St);
    if (!success)
        goto cleanup;

    fmpz_mod_mpoly_from_mpolyn_perm_inflate(G, I->Gbits, ctx, Gn, nctx,
                                       I->brown_perm, I->Gmin_exp, I->Gstride);

    if (Abar != NULL)
        fmpz_mod_mpoly_from_mpolyn_perm_inflate(Abar, I->Abarbits, ctx, Abarn, nctx,
                                    I->brown_perm, I->Abarmin_exp, I->Gstride);

    if (Bbar != NULL)
        fmpz_mod_mpoly_from_mpolyn_perm_inflate(Bbar, I->Bbarbits, ctx, Bbarn, nctx,
                                    I->brown_perm, I->Bbarmin_exp, I->Gstride);

    success = 1;

cleanup:

    fmpz_mod_poly_stack_clear(St->poly_stack);
    fmpz_mod_polyun_stack_clear(St->polyun_stack);
    fmpz_mod_mpolyn_stack_clear(St->mpolyn_stack, nctx);
    fmpz_mod_mpolyn_clear(An, nctx);
    fmpz_mod_mpolyn_clear(Bn, nctx);
    fmpz_mod_mpolyn_clear(Gn, nctx);
    fmpz_mod_mpolyn_clear(Abarn, nctx);
    fmpz_mod_mpolyn_clear(Bbarn, nctx);
    fmpz_mod_mpoly_ctx_clear(nctx);

    return success;
}

static int _try_nmod(
    fmpz_mod_mpoly_t G,
    fmpz_mod_mpoly_t Abar,
    fmpz_mod_mpoly_t Bbar,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    nmod_mpoly_ctx_t nctx;
    nmod_mpoly_t nG, nAbar, nBbar, nA, nB;

    FLINT_ASSERT(fmpz_abs_fits_ui(fmpz_mod_ctx_modulus(ctx->ffinfo)));

    *nctx->minfo = *ctx->minfo;
    nmod_init(&nctx->mod, fmpz_get_ui(fmpz_mod_ctx_modulus(ctx->ffinfo)));
    nmod_mpoly_init(nG, nctx);
    nmod_mpoly_init(nAbar, nctx);
    nmod_mpoly_init(nBbar, nctx);
    nmod_mpoly_init(nA, nctx);
    nmod_mpoly_init(nB, nctx);

    _fmpz_mod_mpoly_get_nmod_mpoly(nA, nctx, A, ctx);
    _fmpz_mod_mpoly_get_nmod_mpoly(nB, nctx, B, ctx);

    success = _nmod_mpoly_gcd_algo_small(nG, Abar == NULL ? NULL : nAbar,
                              Bbar == NULL ? NULL : nBbar, nA, nB, nctx, algo);

    _fmpz_mod_mpoly_set_nmod_mpoly(G, ctx, nG, nctx);
    if (Abar != NULL)
        _fmpz_mod_mpoly_set_nmod_mpoly(Abar, ctx, nAbar, nctx);
    if (Bbar != NULL)
        _fmpz_mod_mpoly_set_nmod_mpoly(Bbar, ctx, nBbar, nctx);

    nmod_mpoly_clear(nG, nctx);
    nmod_mpoly_clear(nAbar, nctx);
    nmod_mpoly_clear(nBbar, nctx);
    nmod_mpoly_clear(nA, nctx);
    nmod_mpoly_clear(nB, nctx);

    return success;
}


/*
    Both A and B have to be packed into bits <= FLINT_BITS
    return is 1 for success, 0 for failure.
*/
static int _fmpz_mod_mpoly_gcd_algo_small(
    fmpz_mod_mpoly_t G,
    fmpz_mod_mpoly_t Abar, /* could be NULL */
    fmpz_mod_mpoly_t Bbar, /* could be NULL */
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong v_in_both;
    slong v_in_either;
    slong v_in_A_only;
    slong v_in_B_only;
    slong j;
    slong nvars = ctx->minfo->nvars;
    mpoly_gcd_info_t I;
#if FLINT_WANT_ASSERT
    fmpz_mod_mpoly_t T, Asave, Bsave;
#endif

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);

    if (A->length == 1)
        return _do_monomial_gcd(G, Bbar, Abar, B, A, ctx);

    if (B->length == 1)
        return _do_monomial_gcd(G, Abar, Bbar, A, B, ctx);

#if FLINT_WANT_ASSERT
    fmpz_mod_mpoly_init(T, ctx);
    fmpz_mod_mpoly_init(Asave, ctx);
    fmpz_mod_mpoly_init(Bsave, ctx);
    fmpz_mod_mpoly_set(Asave, A, ctx);
    fmpz_mod_mpoly_set(Bsave, B, ctx);
#endif

    mpoly_gcd_info_init(I, nvars);

    /* entries of I are all now invalid */

    I->Gbits = FLINT_MIN(A->bits, B->bits);
    I->Abarbits = A->bits;
    I->Bbarbits = B->bits;

    mpoly_gcd_info_limits(I->Amax_exp, I->Amin_exp, I->Alead_count,
                      I->Atail_count, A->exps, A->bits, A->length, ctx->minfo);
    mpoly_gcd_info_limits(I->Bmax_exp, I->Bmin_exp, I->Blead_count,
                      I->Btail_count, B->exps, B->bits, B->length, ctx->minfo);

    /* set ess(p) := p/term_content(p) */

    /* check if the cofactors could be monomials, i.e. ess(A) == ess(B) */
    for (j = 0; j < nvars; j++)
    {
        if (I->Amax_exp[j] - I->Amin_exp[j] != I->Bmax_exp[j] - I->Bmin_exp[j])
            goto skip_monomial_cofactors;
    }

    if (_try_monomial_cofactors(G, Abar, Bbar, A, B, ctx))
        goto successful;

skip_monomial_cofactors:

    if (fmpz_abs_fits_ui(fmpz_mod_ctx_modulus(ctx->ffinfo)))
    {
        success = _try_nmod(G, Abar, Bbar, A, B, ctx, algo);
        goto cleanup;
    }

    mpoly_gcd_info_stride(I->Gstride,
            A->exps, A->bits, A->length, I->Amax_exp, I->Amin_exp,
            B->exps, B->bits, B->length, I->Bmax_exp, I->Bmin_exp, ctx->minfo);

    for (j = 0; j < nvars; j++)
    {
        ulong t = I->Gstride[j];

        if (t == 0)
        {
            FLINT_ASSERT(I->Amax_exp[j] == I->Amin_exp[j] ||
                         I->Bmax_exp[j] == I->Bmin_exp[j]);
        }
        else
        {
            FLINT_ASSERT((I->Amax_exp[j] - I->Amin_exp[j]) % t == 0);
            FLINT_ASSERT((I->Bmax_exp[j] - I->Bmin_exp[j]) % t == 0);
        }

        I->Adeflate_deg[j] = t == 0 ? 0 : (I->Amax_exp[j] - I->Amin_exp[j])/t;
        I->Bdeflate_deg[j] = t == 0 ? 0 : (I->Bmax_exp[j] - I->Bmin_exp[j])/t;

        t = FLINT_MIN(I->Amin_exp[j], I->Bmin_exp[j]);
        I->Gmin_exp[j] = t;
        I->Abarmin_exp[j] = I->Amin_exp[j] - t;
        I->Bbarmin_exp[j] = I->Bmin_exp[j] - t;
    }

    /*
        The following are now valid:
            I->Amax_exp, I->Amin_exp, I->Alead_count, I->Atail_count,
            I->Bmax_exp, I->Bmin_exp, I->Blead_count, I->Btail_count,
            I->Gstride
            I->Adeflate_deg
            I->Bdeflate_deg
            I->Gmin_exp
    */

    /* check if ess(A) and ess(B) have a variable v_in_both in common */
    v_in_both = -WORD(1);
    for (j = 0; j < nvars; j++)
    {
        if (I->Amax_exp[j] > I->Amin_exp[j] && I->Bmax_exp[j] > I->Bmin_exp[j])
        {
            v_in_both = j;
            break;
        }
    }
    if (v_in_both == -WORD(1))
    {
        _do_trivial(G, Abar, Bbar, A, B, I, ctx);
        goto successful;
    }

    /* check if ess(A) and ess(B) depend on another variable v_in_either */
    FLINT_ASSERT(0 <= v_in_both && v_in_both < nvars);

    v_in_either = -WORD(1);
    for (j = 0; j < nvars; j++)
    {
        if (j == v_in_both)
            continue;

        if (I->Amax_exp[j] > I->Amin_exp[j] || I->Bmax_exp[j] > I->Bmin_exp[j])
        {
            v_in_either = j;
            break;
        }
    }

    if (v_in_either == -WORD(1))
    {
        _do_univar(G, Abar, Bbar, A, B, v_in_both, I, ctx);
        goto successful;
    }

    /* check if there is a variable in ess(A) that is not in ess(B) */
    v_in_A_only = -WORD(1);
    v_in_B_only = -WORD(1);
    for (j = 0; j < nvars; j++)
    {
        if (I->Amax_exp[j] > I->Amin_exp[j] && I->Bmax_exp[j] == I->Bmin_exp[j])
        {
            v_in_A_only = j;
            break;
        }
        if (I->Bmax_exp[j] > I->Bmin_exp[j] && I->Amax_exp[j] == I->Amin_exp[j])
        {
            v_in_B_only = j;
            break;
        }
    }
    if (v_in_A_only != -WORD(1))
    {
        success = _try_missing_var(G, I->Gbits, Abar, Bbar, v_in_A_only,
                                   A, I->Amin_exp[v_in_A_only],
                                   B, I->Bmin_exp[v_in_A_only],
                                   ctx);
        goto cleanup;
    }
    if (v_in_B_only != -WORD(1))
    {
        success = _try_missing_var(G, I->Gbits, Bbar, Abar, v_in_B_only,
                                   B, I->Bmin_exp[v_in_B_only],
                                   A, I->Amin_exp[v_in_B_only],
                                   ctx);
        goto cleanup;
    }

    /*
        all variable are now either
            missing from both ess(A) and ess(B), or
            present in both ess(A) and ess(B)
        and there are at least two in the latter case
    */

    mpoly_gcd_info_set_perm(I, A->length, B->length, ctx->minfo);

    /* _set_estimates will probably calculate the correct total degrees */
    I->Adeflate_tdeg = I->Bdeflate_tdeg = -1;

    _set_estimates(I, A, B, ctx);

    /* everything in I is valid now */

    /* check divisibility A/B and B/A */
    {
        int gcd_is_trivial = 1;
        int try_a = I->Gdeflate_deg_bounds_are_nice;
        int try_b = I->Gdeflate_deg_bounds_are_nice;

        for (j = 0; j < nvars; j++)
        {
            if (I->Gdeflate_deg_bound[j] != 0)
                gcd_is_trivial = 0;

            if (I->Adeflate_deg[j] != I->Gdeflate_deg_bound[j])
                try_a = 0;

            if (I->Bdeflate_deg[j] != I->Gdeflate_deg_bound[j])
                try_b = 0;
        }

        if (gcd_is_trivial)
        {
            _do_trivial(G, Abar, Bbar, A, B, I, ctx);
            goto successful;
        }

        if (try_a && _try_divides(G, Bbar, Abar, B, A, ctx))
            goto successful;

        if (try_b && _try_divides(G, Abar, Bbar, A, B, ctx))
            goto successful;
    }

    if (algo == MPOLY_GCD_USE_PRS)
    {
        success = _try_prs(G, Abar, Bbar, A, B, I, ctx);
        goto cleanup;
    }

    if (I->mvars < 3)
    {
        mpoly_gcd_info_measure_brown(I, A->length, B->length, ctx->minfo);
        mpoly_gcd_info_measure_hensel(I, A->length, B->length, ctx->minfo);

        algo &= (MPOLY_GCD_USE_BROWN | MPOLY_GCD_USE_HENSEL);

        if (algo == MPOLY_GCD_USE_BROWN)
        {
            success = _try_brown(G, Abar, Bbar, A, B, I, ctx);
        }
        else if (algo == MPOLY_GCD_USE_HENSEL)
        {
            success = _try_hensel(G, Abar, Bbar, A, B, I, ctx);
        }
        else
        {
            if (I->Adensity + I->Bdensity > 0.05)
            {
                success = _try_brown(G, Abar, Bbar, A, B, I, ctx) ||
                          _try_hensel(G, Abar, Bbar, A, B, I, ctx);
            }
            else
            {
                success = _try_hensel(G, Abar, Bbar, A, B, I, ctx) ||
                          _try_brown(G, Abar, Bbar, A, B, I, ctx);
            }
        }

        goto cleanup;
    }
    else if (algo == MPOLY_GCD_USE_HENSEL)
    {
        mpoly_gcd_info_measure_hensel(I, A->length, B->length, ctx->minfo);
        success = _try_hensel(G, Abar, Bbar, A, B, I, ctx);
        goto cleanup;
    }
    else if (algo == MPOLY_GCD_USE_BROWN)
    {
        mpoly_gcd_info_measure_brown(I, A->length, B->length, ctx->minfo);
        success = _try_brown(G, Abar, Bbar, A, B, I, ctx);
        goto cleanup;
    }
    else if (algo == MPOLY_GCD_USE_ZIPPEL)
    {
        mpoly_gcd_info_measure_zippel(I, A->length, B->length, ctx->minfo);
        success = _try_zippel(G, Abar, Bbar, A, B, I, ctx);
        goto cleanup;
    }
    else if (algo == MPOLY_GCD_USE_ZIPPEL2)
    {
        mpoly_gcd_info_measure_zippel2(I, A->length, B->length, ctx->minfo);
        success = _try_zippel2(G, Abar, Bbar, A, B, I, ctx);
        goto cleanup;
    }
    else
    {
        slong k;
        int deg_is_small = 1;
        double density = I->Adensity + I->Bdensity;

        /*
            mpoly gcd case.
            Only rule is that measure_X must be called before
                try_X is called or I->X_perm is accessed.
        */

        mpoly_gcd_info_measure_hensel(I, A->length, B->length, ctx->minfo);
        mpoly_gcd_info_measure_brown(I, A->length, B->length, ctx->minfo);
        mpoly_gcd_info_measure_zippel(I, A->length, B->length, ctx->minfo);
        mpoly_gcd_info_measure_zippel2(I, A->length, B->length, ctx->minfo);

        if (density > 0.08)
        {
            if (!deg_is_small && _try_hensel(G, Abar, Bbar, A, B, I, ctx))
                goto successful;

            if (_try_brown(G, Abar, Bbar, A, B, I, ctx))
                goto successful;
        }

        if (I->Adeflate_tdeg > 0 && I->Bdeflate_tdeg > 0)
        {
            fmpz_t x;
            double tdensity;
            fmpz_init(x);
            fmpz_bin_uiui(x, (ulong)I->Adeflate_tdeg + I->mvars, I->mvars);
            tdensity = A->length/fmpz_get_d(x);
            fmpz_bin_uiui(x, (ulong)I->Bdeflate_tdeg + I->mvars, I->mvars);
            tdensity += B->length/fmpz_get_d(x);
            density = FLINT_MAX(density, tdensity);
            fmpz_clear(x);
        }

        if (density > (deg_is_small ? 0.05 : 0.001))
        {
            if (_try_hensel(G, Abar, Bbar, A, B, I, ctx))
                goto successful;                
        }

        k = I->zippel2_perm[1];
        k = FLINT_MAX(I->Adeflate_deg[k], I->Bdeflate_deg[k]);
        if ((A->length + B->length)/64 < k)
        {
            if (_try_zippel(G, Abar, Bbar, A, B, I, ctx))
                goto successful;
            if (_try_zippel2(G, Abar, Bbar, A, B, I, ctx))
                goto successful;
        }
        else
        {
            if (_try_zippel2(G, Abar, Bbar, A, B, I, ctx))
                goto successful;
            if (_try_zippel(G, Abar, Bbar, A, B, I, ctx))
                goto successful;
        }

        success = _try_hensel(G, Abar, Bbar, A, B, I, ctx) ||
                  _try_brown(G, Abar, Bbar, A, B, I, ctx);
        goto cleanup;
    }

    success = 0;
    goto cleanup;

successful:

    success = 1;

cleanup:

    mpoly_gcd_info_clear(I);

    if (success)
    {
        FLINT_ASSERT(G->length > 0);

        if (!fmpz_is_one(G->coeffs + 0))
        {
            if (Abar != NULL)
                _fmpz_mod_vec_scalar_mul_fmpz_mod(Abar->coeffs, Abar->coeffs,
                                     Abar->length, G->coeffs + 0, ctx->ffinfo);
            if (Bbar != NULL)
                _fmpz_mod_vec_scalar_mul_fmpz_mod(Bbar->coeffs, Bbar->coeffs,
                                     Bbar->length, G->coeffs + 0, ctx->ffinfo);

            _fmpz_mod_vec_scalar_div_fmpz_mod(G->coeffs, G->coeffs, G->length,
                                                   G->coeffs + 0, ctx->ffinfo);
        }

        FLINT_ASSERT(fmpz_mod_mpoly_divides(T, Asave, G, ctx));
        FLINT_ASSERT(Abar == NULL || fmpz_mod_mpoly_equal(T, Abar, ctx));
            
        FLINT_ASSERT(fmpz_mod_mpoly_divides(T, Bsave, G, ctx));
        FLINT_ASSERT(Bbar == NULL || fmpz_mod_mpoly_equal(T, Bbar, ctx));
    }

#if FLINT_WANT_ASSERT
    fmpz_mod_mpoly_clear(T, ctx);
    fmpz_mod_mpoly_clear(Asave, ctx);
    fmpz_mod_mpoly_clear(Bsave, ctx);
#endif

    return success;
}


/*
    The gcd calculation is unusual.
    First see if both inputs fit into FLINT_BITS.
    Then, try deflation as a last resort.
*/
static int _fmpz_mod_mpoly_gcd_algo_large(
    fmpz_mod_mpoly_t G,
    fmpz_mod_mpoly_t Abar,
    fmpz_mod_mpoly_t Bbar,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong k;
    fmpz * Ashift, * Astride;
    fmpz * Bshift, * Bstride;
    fmpz * Gshift, * Gstride;
    fmpz_mod_mpoly_t Anew, Bnew;
    const fmpz_mod_mpoly_struct * Ause, * Buse;

    if (A->length == 1)
        return _do_monomial_gcd(G, Bbar, Abar, B, A, ctx);

    if (B->length == 1)
        return _do_monomial_gcd(G, Abar, Bbar, A, B, ctx);

    if (_try_monomial_cofactors(G, Abar, Bbar, A, B, ctx))
        return 1;

    fmpz_mod_mpoly_init(Anew, ctx);
    fmpz_mod_mpoly_init(Bnew, ctx);

    Ause = A;
    if (A->bits > FLINT_BITS)
    {
        if (!fmpz_mod_mpoly_repack_bits(Anew, A, FLINT_BITS, ctx))
            goto could_not_repack;
        Ause = Anew;
    }

    Buse = B;
    if (B->bits > FLINT_BITS)
    {
        if (!fmpz_mod_mpoly_repack_bits(Bnew, B, FLINT_BITS, ctx))
            goto could_not_repack;
        Buse = Bnew;
    }

    success = _fmpz_mod_mpoly_gcd_algo_small(G, Abar, Bbar, Ause, Buse, ctx, algo);
    goto cleanup;

could_not_repack:

    /*
        One of A or B could not be repacked into FLINT_BITS. See if
        they both fit into FLINT_BITS after deflation.
    */

    Ashift  = _fmpz_vec_init(ctx->minfo->nvars);
    Astride = _fmpz_vec_init(ctx->minfo->nvars);
    Bshift  = _fmpz_vec_init(ctx->minfo->nvars);
    Bstride = _fmpz_vec_init(ctx->minfo->nvars);
    Gshift  = _fmpz_vec_init(ctx->minfo->nvars);
    Gstride = _fmpz_vec_init(ctx->minfo->nvars);

    fmpz_mod_mpoly_deflation(Ashift, Astride, A, ctx);
    fmpz_mod_mpoly_deflation(Bshift, Bstride, B, ctx);
    _fmpz_vec_min(Gshift, Ashift, Bshift, ctx->minfo->nvars);
    for (k = 0; k < ctx->minfo->nvars; k++)
        fmpz_gcd(Gstride + k, Astride + k, Bstride + k);

    fmpz_mod_mpoly_deflate(Anew, A, Ashift, Gstride, ctx);
    if (Anew->bits > FLINT_BITS)
    {
        success = fmpz_mod_mpoly_repack_bits(Anew, Anew, FLINT_BITS, ctx);
        if (!success)
            goto deflate_cleanup;
    }

    fmpz_mod_mpoly_deflate(Bnew, B, Bshift, Gstride, ctx);
    if (Bnew->bits > FLINT_BITS)
    {
        success = fmpz_mod_mpoly_repack_bits(Bnew, Bnew, FLINT_BITS, ctx);
        if (!success)
            goto deflate_cleanup;
    }

    success = _fmpz_mod_mpoly_gcd_algo_small(G, Abar, Bbar, Anew, Bnew, ctx, algo);
    if (!success)
        goto deflate_cleanup;

    for (k = 0; k < ctx->minfo->nvars; k++)
    {
        fmpz_sub(Ashift + k, Ashift + k, Gshift + k);
        fmpz_sub(Bshift + k, Bshift + k, Gshift + k);
        FLINT_ASSERT(fmpz_sgn(Ashift + k) >= 0);
        FLINT_ASSERT(fmpz_sgn(Bshift + k) >= 0);
    }

    fmpz_mod_mpoly_inflate(G, G, Gshift, Gstride, ctx);
    if (Abar != NULL)
        fmpz_mod_mpoly_inflate(Abar, Abar, Ashift, Gstride, ctx);
    if (Bbar != NULL)
        fmpz_mod_mpoly_inflate(Bbar, Bbar, Bshift, Gstride, ctx);

    FLINT_ASSERT(G->length > 0);
    if (!fmpz_is_one(G->coeffs + 0))
    {
        if (Abar != NULL)
            _fmpz_mod_vec_scalar_mul_fmpz_mod(Abar->coeffs, Abar->coeffs,
                                     Abar->length, G->coeffs + 0, ctx->ffinfo);

        if (Bbar != NULL)
            _fmpz_mod_vec_scalar_mul_fmpz_mod(Bbar->coeffs, Bbar->coeffs,
                                     Bbar->length, G->coeffs + 0, ctx->ffinfo);

        _fmpz_mod_vec_scalar_div_fmpz_mod(G->coeffs, G->coeffs, G->length,
                                                   G->coeffs + 0, ctx->ffinfo);
    }

deflate_cleanup:

    _fmpz_vec_clear(Ashift, ctx->minfo->nvars);
    _fmpz_vec_clear(Astride, ctx->minfo->nvars);
    _fmpz_vec_clear(Bshift, ctx->minfo->nvars);
    _fmpz_vec_clear(Bstride, ctx->minfo->nvars);
    _fmpz_vec_clear(Gshift, ctx->minfo->nvars);
    _fmpz_vec_clear(Gstride, ctx->minfo->nvars);

cleanup:

    fmpz_mod_mpoly_clear(Anew, ctx);
    fmpz_mod_mpoly_clear(Bnew, ctx);

    return success;
}


int _fmpz_mod_mpoly_gcd_algo(
    fmpz_mod_mpoly_t G,
    fmpz_mod_mpoly_t Abar,
    fmpz_mod_mpoly_t Bbar,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    FLINT_ASSERT(!fmpz_mod_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!fmpz_mod_mpoly_is_zero(B, ctx));

    if (A->bits <= FLINT_BITS && B->bits <= FLINT_BITS)
        return _fmpz_mod_mpoly_gcd_algo_small(G, Abar, Bbar, A, B, ctx, algo);
    else
        return _fmpz_mod_mpoly_gcd_algo_large(G, Abar, Bbar, A, B, ctx, algo);
}

