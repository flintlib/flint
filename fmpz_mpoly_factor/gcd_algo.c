/*
    Copyright (C) 2018-2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"

/*
    For each j, set out[j] to the evaluation of A at x_i = alpha[i] (i != j)
    i.e. if nvars = 3
        out[0] = A(x, alpha[1], alpha[2])
        out[1] = A(alpha[0], x, alpha[2])
        out[2] = A(alpha[0], alpha[1], x)

    If ignore[j] is nonzero, then out[j] need not be calculated, probably
    because we shouldn't calculate it in dense form.
*/
void fmpz_mpoly_evals(
    nmod_poly_struct * out,
    const int * ignore,
    const fmpz_mpoly_t A,
    ulong * Amin_exp,
    ulong * Amax_exp,
    ulong * Astride,
    mp_limb_t * alpha,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    slong nvars = ctx->minfo->nvars;
    slong total_limit, total_length;
    int use_direct_LUT;
    ulong varexp;
    ulong mask;
    slong * offsets, * shifts;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    ulong * Aexp = A->exps;
    fmpz * Acoeff = A->coeffs;
    mp_limb_t meval;
    mp_limb_t t;

    FLINT_ASSERT(A->bits <= FLINT_BITS);

    mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    offsets = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    shifts = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));

    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        nmod_poly_zero(out + j);
        mpoly_gen_offset_shift_sp(offsets + j, shifts + j, j, A->bits, ctx->minfo);
    }

    /*
        two cases:
        (1) the Amax_exp[j] are small enough to calculate a direct LUT
        (2) use a LUT for exponents that are powers of two
    */

    total_limit = A->length/256;
    total_limit = FLINT_MAX(WORD(9999), total_limit);
    total_length = 0;
    use_direct_LUT = 1;
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        total_length += Amax_exp[j] + 1;
        if ((ulong) total_length > (ulong) total_limit)
            use_direct_LUT = 0;
    }

    if (use_direct_LUT)
    {
        slong off;
        mp_limb_t * LUT, ** LUTvalue, ** LUTvalueinv;

        /* value of powers of alpha[j] */
        LUT = (mp_limb_t *) flint_malloc(2*total_length*sizeof(mp_limb_t));

        /* pointers into LUT */
        LUTvalue    = (mp_limb_t **) flint_malloc(nvars*sizeof(mp_limb_t *));
        LUTvalueinv = (mp_limb_t **) flint_malloc(nvars*sizeof(mp_limb_t *));

        off = 0;
        for (j = 0; j < nvars; j++)
        {
            ulong k;
            mp_limb_t alphainvj = nmod_inv(alpha[j], (out + 0)->mod);

            LUTvalue[j] = LUT + off;
            LUTvalueinv[j] = LUT + total_length + off;
            LUTvalue[j][0] = 1;
            LUTvalueinv[j][0] = 1;
            for (k = 0; k < Amax_exp[j]; k++)
            {
                LUTvalue[j][k + 1] = nmod_mul(LUTvalue[j][k], alpha[j],
                                                               (out + 0)->mod);
                LUTvalueinv[j][k + 1] = nmod_mul(LUTvalueinv[j][k], alphainvj,
                                                               (out + 0)->mod);
            }

            off += Amax_exp[j] + 1;
        }
        FLINT_ASSERT(off == total_length);

        for (i = 0; i < A->length; i++)
        {
            meval = fmpz_fdiv_ui(Acoeff + i, out->mod.n);

            for (j = 0; j < nvars; j++)
            {
                varexp = ((Aexp + N*i)[offsets[j]]>>shifts[j])&mask;
                FLINT_ASSERT(varexp <= Amax_exp[j]);
                meval = nmod_mul(meval, LUTvalue[j][varexp], (out + 0)->mod);
            }

            for (j = 0; j < nvars; j++)
            {
                varexp = ((Aexp + N*i)[offsets[j]]>>shifts[j])&mask;

                if (ignore[j])
                    continue;

                t = nmod_mul(meval, LUTvalueinv[j][varexp], (out + j)->mod);

                FLINT_ASSERT((Astride[j] == 0 && varexp == Amin_exp[j])
                                  || (varexp - Amin_exp[j]) % Astride[j] == 0);

                varexp = Astride[j] < 2 ? varexp - Amin_exp[j] :
                                           (varexp - Amin_exp[j])/Astride[j];

                t = nmod_add(t, nmod_poly_get_coeff_ui(out + j, varexp),
                                                               (out + j)->mod);
                nmod_poly_set_coeff_ui(out + j, varexp, t);
            }
        }

        flint_free(LUT);
        flint_free(LUTvalue);
        flint_free(LUTvalueinv);
    }
    else
    {
        slong LUTlen;
        ulong * LUTmask;
        slong * LUToffset, * LUTvar;
        mp_limb_t * LUTvalue, * LUTvalueinv;
        mp_limb_t * vieval;
        mp_limb_t t, xpoweval, xinvpoweval;

        LUToffset   = (slong *) flint_malloc(N*FLINT_BITS*sizeof(slong));
        LUTmask     = (ulong *) flint_malloc(N*FLINT_BITS*sizeof(ulong));
        LUTvalue    = (mp_limb_t *) flint_malloc(N*FLINT_BITS*sizeof(mp_limb_t));
        LUTvar      = (slong *) flint_malloc(N*FLINT_BITS*sizeof(slong));
        LUTvalueinv = (mp_limb_t *) flint_malloc(N*FLINT_BITS*sizeof(mp_limb_t));

        vieval = (mp_limb_t *) flint_malloc(nvars*sizeof(mp_limb_t));

        LUTlen = 0;
        for (j = nvars - 1; j >= 0; j--)
        {
            flint_bitcnt_t bits = FLINT_BIT_COUNT(Amax_exp[j]);
            xpoweval = alpha[j]; /* xpoweval = alpha[j]^(2^i) */
            xinvpoweval = nmod_inv(xpoweval, (out + 0)->mod); /* alpha[j]^(-2^i) */
            for (i = 0; i < bits; i++)
            {
                LUToffset[LUTlen] = offsets[j];
                LUTmask[LUTlen] = (UWORD(1) << (shifts[j] + i));
                LUTvalue[LUTlen] = xpoweval;
                LUTvalueinv[LUTlen] = xinvpoweval;
                LUTvar[LUTlen] = j;
                LUTlen++;
                xpoweval = nmod_mul(xpoweval, xpoweval, (out + 0)->mod);
                xinvpoweval = nmod_mul(xinvpoweval, xinvpoweval, (out + 0)->mod);
            }

            vieval[j] = 1;
        }
        FLINT_ASSERT(LUTlen < N*FLINT_BITS);

        for (i = 0; i < A->length; i++)
        {
            meval = fmpz_fdiv_ui(Acoeff + i, (out + 0)->mod.n);

            for (j = 0; j < LUTlen; j++)
            {
                if (((Aexp + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
                {
                    meval = nmod_mul(meval, LUTvalue[j], (out + 0)->mod);
                    vieval[LUTvar[j]] = nmod_mul(vieval[LUTvar[j]],
                                               LUTvalueinv[j], (out + 0)->mod);
                }
            }

            for (j = 0; j < nvars; j++)
            {
                varexp = ((Aexp + N*i)[offsets[j]]>>shifts[j])&mask;

                FLINT_ASSERT((Astride[j] == 0 && varexp == Amin_exp[j])
                                  || (varexp - Amin_exp[j]) % Astride[j] == 0);

                varexp = Astride[j] < 2 ? varexp - Amin_exp[j] :
                                           (varexp - Amin_exp[j])/Astride[j];

                t = nmod_mul(meval, vieval[j], (out + j)->mod);
                t = nmod_add(t, nmod_poly_get_coeff_ui(out + j, varexp),
                                                               (out + j)->mod);
                nmod_poly_set_coeff_ui(out + j, varexp, t);
                vieval[j] = 1;
            }
        }

        flint_free(LUToffset);
        flint_free(LUTmask);
        flint_free(LUTvalue);
        flint_free(LUTvar);
        flint_free(LUTvalueinv);

        flint_free(vieval);
    }

    flint_free(offsets);
    flint_free(shifts);
}


void _set_estimates(
    mpoly_gcd_info_t I,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    int try_count = 0;
    slong i, j;
    nmod_poly_t Geval;
    nmod_poly_struct * Aevals, * Bevals;
    mp_limb_t p = UWORD(1) << (FLINT_BITS - 1);
    mp_limb_t * alpha;
    flint_rand_t randstate;
    slong ignore_limit;
    int * ignore;

    flint_randinit(randstate);

    ignore = (int *) flint_malloc(ctx->minfo->nvars*sizeof(int));
    alpha = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
    Aevals = (nmod_poly_struct *) flint_malloc(
                                   ctx->minfo->nvars*sizeof(nmod_poly_struct));
    Bevals = (nmod_poly_struct *) flint_malloc(
                                   ctx->minfo->nvars*sizeof(nmod_poly_struct));

    nmod_poly_init(Geval, p);
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        nmod_poly_init(Aevals + j, p);
        nmod_poly_init(Bevals + j, p);
    }

    ignore_limit = A->length/4096 + B->length/4096;
    ignore_limit = FLINT_MAX(WORD(9999), ignore_limit);
    I->Gdeflate_deg_bounds_are_nice = 1;
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        if (   I->Adeflate_deg[j] > ignore_limit
            || I->Bdeflate_deg[j] > ignore_limit)
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
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            I->Gdeflate_deg_bound[j] = FLINT_MIN(I->Adeflate_deg[j],
                                                 I->Bdeflate_deg[j]);
            I->Gterm_count_est[j] = (I->Gdeflate_deg_bound[j] + 1)/2;
        }

        goto cleanup;
    }

    p = n_nextprime(p, 1);
    nmod_init(&Geval->mod, p);
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        alpha[j] = n_urandint(randstate, p - 1) + 1;
        nmod_init(&(Aevals + j)->mod, p);
        nmod_init(&(Bevals + j)->mod, p);
    }

    fmpz_mpoly_evals(Aevals, ignore, A, I->Amin_exp, I->Amax_exp, I->Gstride, alpha, ctx);
    fmpz_mpoly_evals(Bevals, ignore, B, I->Bmin_exp, I->Bmax_exp, I->Gstride, alpha, ctx);

    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        if (ignore[j])
        {
            I->Gdeflate_deg_bound[j] = FLINT_MIN(I->Adeflate_deg[j],
                                                 I->Bdeflate_deg[j]);
            I->Gterm_count_est[j] = (I->Gdeflate_deg_bound[j] + 1)/2;
        }
        else
        {
            if (   I->Adeflate_deg[j] != nmod_poly_degree(Aevals + j)
                || I->Bdeflate_deg[j] != nmod_poly_degree(Bevals + j))
            {
                goto try_again;
            }

            nmod_poly_gcd(Geval, Aevals + j, Bevals + j);

            I->Gterm_count_est[j] = 0;
            I->Gdeflate_deg_bound[j] = nmod_poly_degree(Geval);
            for (i = I->Gdeflate_deg_bound[j]; i >= 0; i--)
            {
                I->Gterm_count_est[j] += (Geval->coeffs[i] != 0);
            }
        }
    }

cleanup:

    nmod_poly_clear(Geval);
    for (j = 0; j < ctx->minfo->nvars; j++)
    {
        nmod_poly_clear(Aevals + j);
        nmod_poly_clear(Bevals + j);
    }

    flint_free(ignore);
    flint_free(alpha);
    flint_free(Aevals);
    flint_free(Bevals);

    flint_randclear(randstate);

    return;
}

/* (Abar, Bbar) = (A, B) */
static void _parallel_set(
    fmpz_mpoly_t Abar, /* could be NULL */
    fmpz_mpoly_t Bbar, /* could be NULL */
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    if (Abar == B && Bbar == A)
    {
        FLINT_ASSERT(Abar != NULL && Bbar != NULL);
        fmpz_mpoly_set(Abar, B, ctx);
        fmpz_mpoly_set(Bbar, A, ctx);
        fmpz_mpoly_swap(Abar, Bbar, ctx);
    }
    else if (Abar == B && Bbar != A)
    {
        FLINT_ASSERT(Abar != NULL);
        if (Bbar != NULL)
            fmpz_mpoly_set(Bbar, B, ctx);
        fmpz_mpoly_set(Abar, A, ctx);
    }
    else
    {
        if (Abar != NULL)
            fmpz_mpoly_set(Abar, A, ctx);
        if (Bbar != NULL)
            fmpz_mpoly_set(Bbar, B, ctx);
    }
}

/*
    The variables in ess(A) and ess(B) are disjoint.
    gcd is trivial to compute.
*/
static void _do_trivial(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,  /* could be NULL */
    fmpz_mpoly_t Bbar,  /* could be NULL */
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const mpoly_gcd_info_t I,
    const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t cG;

    fmpz_init(cG);
    _fmpz_vec_content(cG, A->coeffs, A->length);
    _fmpz_vec_content_chained(cG, B->coeffs, B->length, cG);

    _parallel_set(Abar, Bbar, A, B, ctx);

    if (Abar != NULL)
    {
        _fmpz_vec_scalar_divexact_fmpz(Abar->coeffs, Abar->coeffs,
                                                             Abar->length, cG);
        mpoly_monomials_shift_right_ui(Abar->exps, Abar->bits, Abar->length,
                                                      I->Gmin_exp, ctx->minfo);
    }

    if (Bbar != NULL)
    {
        _fmpz_vec_scalar_divexact_fmpz(Bbar->coeffs, Bbar->coeffs,
                                                             Bbar->length, cG);
        mpoly_monomials_shift_right_ui(Bbar->exps, Bbar->bits, Bbar->length,
                                                      I->Gmin_exp, ctx->minfo);
    }

    fmpz_mpoly_fit_length_reset_bits(G, 1, I->Gbits, ctx);
    mpoly_set_monomial_ui(G->exps, I->Gmin_exp, I->Gbits, ctx->minfo);
    fmpz_swap(G->coeffs + 0, cG);
    _fmpz_mpoly_set_length(G, 1, ctx);

    fmpz_clear(cG);
}


/*********************** Easy when B is a monomial ***************************/
static int _do_monomial_gcd(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,  /* could be NULL */
    fmpz_mpoly_t Bbar,  /* could be NULL */
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_t g;
    flint_bitcnt_t Gbits = FLINT_MIN(A->bits, B->bits);
    fmpz * minAfields, * minAdegs, * minBdegs;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length == 1);

    fmpz_init(g);

    /* compute the coefficient of G */
    _fmpz_vec_content_chained(g, A->coeffs, A->length, B->coeffs + 0);
    fmpz_abs(g, g);

    if (B->exps[0] == 0 && mpoly_monomial_is_zero(B->exps,
                                     mpoly_words_per_exp(B->bits, ctx->minfo)))
    {
        _parallel_set(Abar, Bbar, A, B, ctx);

        if (Abar != NULL && !fmpz_is_one(g))
            _fmpz_vec_scalar_divexact_fmpz(Abar->coeffs, Abar->coeffs,
                                                              Abar->length, g);
        if (Bbar != NULL && !fmpz_is_one(g))
            _fmpz_vec_scalar_divexact_fmpz(Bbar->coeffs, Bbar->coeffs,
                                                              Bbar->length, g);
        fmpz_mpoly_fit_length(G, 1, ctx);
        mpoly_monomial_zero(G->exps, mpoly_words_per_exp(G->bits, ctx->minfo));

        goto done;
    }

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
    {
        _fmpz_vec_scalar_divexact_fmpz(Abar->coeffs, Abar->coeffs, Abar->length, g);
        mpoly_monomials_shift_right_ffmpz(Abar->exps, Abar->bits, Abar->length,
                                                         minBdegs, ctx->minfo);
    }

    if (Bbar != NULL)
    {
        _fmpz_vec_scalar_divexact_fmpz(Bbar->coeffs, Bbar->coeffs, Bbar->length, g);
        mpoly_monomials_shift_right_ffmpz(Bbar->exps, Bbar->bits, Bbar->length,
                                                         minBdegs, ctx->minfo);
    }

    /* write G */
    fmpz_mpoly_fit_length_reset_bits(G, 1, Gbits, ctx);
    mpoly_set_monomial_ffmpz(G->exps, minBdegs, Gbits, ctx->minfo);

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

done:

    fmpz_swap(G->coeffs + 0, g);
    _fmpz_mpoly_set_length(G, 1, ctx);

    fmpz_clear(g);

    return 1;
}


/********************** See if cofactors are monomials ***********************/
static int _try_monomial_cofactors(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,  /* could be NULL */
    fmpz_mpoly_t Bbar,  /* could be NULL */
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    slong NA, NG;
    slong nvars = ctx->minfo->nvars;
    fmpz * Abarexps, * Bbarexps, * Texps;
    fmpz_t t1, t2;
    fmpz_t gA, gB;
    fmpz_mpoly_t T;
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
    fmpz_init_set(gA, A->coeffs + 0);
    fmpz_init_set(gB, B->coeffs + 0);

    for (i = A->length - 1; i > 0; i--)
    {
        fmpz_mul(t1, A->coeffs + 0, B->coeffs + i);
        fmpz_mul(t2, B->coeffs + 0, A->coeffs + i);
        success = fmpz_equal(t1, t2);
        if (!success)
            goto cleanup;

        fmpz_gcd(gA, gA, A->coeffs + i);
        fmpz_gcd(gB, gB, B->coeffs + i);
    }

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

    /* put A's cofactor coefficient in t1 */
    fmpz_gcd(t2, gA, gB);
    fmpz_divexact(t1, gA, t2);
    if (fmpz_sgn(A->coeffs + 0) < 0)
        fmpz_neg(t1, t1);

    /* put B's cofactor coefficient in t2 */
    fmpz_divexact(gA, A->coeffs + 0, t1);
    fmpz_divexact(t2, B->coeffs + 0, gA);

    fmpz_mpoly_init3(T, A->length, Gbits, ctx);
    NG = mpoly_words_per_exp(Gbits, ctx->minfo);
    NA = mpoly_words_per_exp(A->bits, ctx->minfo);
    T->length = A->length;
    for (i = 0; i < A->length; i++)
    {
        mpoly_get_monomial_ffmpz(Texps, A->exps + NA*i, A->bits, ctx->minfo);
        _fmpz_vec_sub(Texps, Texps, Abarexps, nvars);
        mpoly_set_monomial_ffmpz(T->exps + NG*i, Texps, Gbits, ctx->minfo);
        fmpz_divexact(T->coeffs + i, A->coeffs + i, t1);
    }
    fmpz_mpoly_swap(G, T, ctx);
    fmpz_mpoly_clear(T, ctx);

    if (Abar != NULL)
    {
        fmpz_mpoly_fit_length_reset_bits(Abar, 1, Abarbits, ctx);
        mpoly_set_monomial_ffmpz(Abar->exps, Abarexps, Abarbits, ctx->minfo);
        fmpz_swap(Abar->coeffs + 0, t1);
        _fmpz_mpoly_set_length(Abar, 1, ctx);
    }

    if (Bbar != NULL)
    {
        fmpz_mpoly_fit_length_reset_bits(Bbar, 1, Bbarbits, ctx);
        mpoly_set_monomial_ffmpz(Bbar->exps, Bbarexps, Bbarbits, ctx->minfo);
        fmpz_swap(Bbar->coeffs + 0, t2);
        _fmpz_mpoly_set_length(Bbar, 1, ctx);
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

cleanup:

    fmpz_clear(t1);
    fmpz_clear(t2);
    fmpz_clear(gA);
    fmpz_clear(gB);

    return success;
}


/***  ess(A) and ess(B) depend on only one variable v_in_both ****************/
static int _do_univar(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    slong v_in_both,
    const mpoly_gcd_info_t I,
    const fmpz_mpoly_ctx_t ctx)
{
    fmpz_poly_t a, b, g, t, r;

    fmpz_poly_init(a);
    fmpz_poly_init(b);
    fmpz_poly_init(g);
    fmpz_poly_init(t);
    fmpz_poly_init(r);

    _fmpz_mpoly_to_fmpz_poly_deflate(a, A, v_in_both, I->Amin_exp, I->Gstride, ctx);
    _fmpz_mpoly_to_fmpz_poly_deflate(b, B, v_in_both, I->Bmin_exp, I->Gstride, ctx);

    fmpz_poly_gcd(g, a, b);

    _fmpz_mpoly_from_fmpz_poly_inflate(G, I->Gbits, g, v_in_both,
                                                 I->Gmin_exp, I->Gstride, ctx);

    if (Abar != NULL)
    {
        fmpz_poly_divrem(t, r, a, g);
        _fmpz_mpoly_from_fmpz_poly_inflate(Abar, I->Abarbits, t, v_in_both,
                                              I->Abarmin_exp, I->Gstride, ctx);
    }

    if (Bbar != NULL)
    {
        fmpz_poly_divrem(t, r, b, g);
        _fmpz_mpoly_from_fmpz_poly_inflate(Bbar, I->Bbarbits, t, v_in_both,
                                              I->Bbarmin_exp, I->Gstride, ctx);
    }

    fmpz_poly_clear(a);
    fmpz_poly_clear(b);
    fmpz_poly_clear(g);
    fmpz_poly_clear(t);
    fmpz_poly_clear(r);

    return 1;
}



/********* Assume B has length one when converted to univar format ***********/
static int _try_missing_var(
    fmpz_mpoly_t G, flint_bitcnt_t Gbits,
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    slong var,
    const fmpz_mpoly_t A, ulong Ashift,
    const fmpz_mpoly_t B, ulong Bshift,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    fmpz_mpoly_t tG;
    fmpz_mpoly_univar_t Au;

    fmpz_mpoly_init(tG, ctx);
    fmpz_mpoly_univar_init(Au, ctx);

#if FLINT_WANT_ASSERT
    fmpz_mpoly_to_univar(Au, B, var, ctx);
    FLINT_ASSERT(Au->length == 1);
#endif

    fmpz_mpoly_to_univar(Au, A, var, ctx);

    fmpz_mpoly_univar_fit_length(Au, Au->length + 1, ctx);
    fmpz_mpoly_set(Au->coeffs + Au->length, B, ctx);
    Au->length++;

    if (Abar == NULL && Bbar == NULL)
    {
        success = _fmpz_mpoly_vec_content_mpoly(G, Au->coeffs, Au->length, ctx);
        if (!success)
            goto cleanup;

        fmpz_mpoly_repack_bits_inplace(G, Gbits, ctx);
        _mpoly_gen_shift_left(G->exps, G->bits, G->length,
                                   var, FLINT_MIN(Ashift, Bshift), ctx->minfo);
    }
    else
    {
        success = _fmpz_mpoly_vec_content_mpoly(tG, Au->coeffs, Au->length, ctx);
        if (!success)
            goto cleanup;

        fmpz_mpoly_repack_bits_inplace(tG, Gbits, ctx);
        _mpoly_gen_shift_left(tG->exps, tG->bits, tG->length,
                                   var, FLINT_MIN(Ashift, Bshift), ctx->minfo);

        if (Abar != NULL)
        {
            success = fmpz_mpoly_divides(Abar, A, tG, ctx);
            FLINT_ASSERT(success);
        }

        if (Bbar != NULL)
        {
            success = fmpz_mpoly_divides(Bbar, Au->coeffs + Au->length - 1, tG, ctx);
            FLINT_ASSERT(success);
        }

        fmpz_mpoly_swap(G, tG, ctx);
    }

    success = 1;

cleanup:

    fmpz_mpoly_clear(tG, ctx);
    fmpz_mpoly_univar_clear(Au, ctx);

    return success;
}


/************************ See if B divides A ********************************/
static int _try_divides(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t BB,
    const fmpz_mpoly_ctx_t ctx)
{
    int success = 0;
    fmpz_mpoly_t Q, B, M;

    fmpz_mpoly_init(Q, ctx);
    fmpz_mpoly_init(B, ctx);
    fmpz_mpoly_init(M, ctx);

    /* BB = M*B */
    fmpz_mpoly_term_content(M, BB, ctx);
    fmpz_mpoly_divides(B, BB, M, ctx);

    success = fmpz_mpoly_divides(Q, A, B, ctx);
    if (success)
    {
        _do_monomial_gcd(G, Abar, Bbar, Q, M, ctx);
        fmpz_mpoly_mul(G, G, B, ctx);
    }

    fmpz_mpoly_clear(Q, ctx);
    fmpz_mpoly_clear(B, ctx);
    fmpz_mpoly_clear(M, ctx);

    return success;
}


/********************* Hit A and B with a prs ********************************/
static int _try_prs(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const mpoly_gcd_info_t I,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong j, var = -WORD(1);
    ulong score, best_score = WORD(1) << 24;
    fmpz_mpoly_t Ac, Bc, Gc, s, t;
    fmpz_mpoly_univar_t Ax, Bx, Gx;

    /*
        Pick a main variable. TODO:
            (1) consider the expected degree of G, and
            (2) add a size limit to univar_pseudo_gcd to abort the calculation
                if things are getting too big.
    */
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

        score = FLINT_BIT_COUNT(I->Alead_count[j] - 1);
        score *= FLINT_BIT_COUNT(I->Blead_count[j] - 1);
        score *= FLINT_BIT_COUNT(I->Atail_count[j] - 1);
        score *= FLINT_BIT_COUNT(I->Btail_count[j] - 1);
        score = FLINT_MAX(score, 1);
        if (n_mul_checked(&score, score, I->Amax_exp[j]))
            continue;
        if (n_mul_checked(&score, score, I->Bmax_exp[j]))
            continue;

        if (score < best_score)
        {
            best_score = score;
            var = j;
        }
    }

    if (var < 0)
        return 0;

    fmpz_mpoly_init(Ac, ctx);
    fmpz_mpoly_init(Bc, ctx);
    fmpz_mpoly_init(Gc, ctx);
    fmpz_mpoly_init(s, ctx);
    fmpz_mpoly_init(t, ctx);
    fmpz_mpoly_univar_init(Ax, ctx);
    fmpz_mpoly_univar_init(Bx, ctx);
    fmpz_mpoly_univar_init(Gx, ctx);

    fmpz_mpoly_to_univar(Ax, A, var, ctx);
    fmpz_mpoly_to_univar(Bx, B, var, ctx);

    success = _fmpz_mpoly_vec_content_mpoly(Ac, Ax->coeffs, Ax->length, ctx) &&
              _fmpz_mpoly_vec_content_mpoly(Bc, Bx->coeffs, Bx->length, ctx) &&
              fmpz_mpoly_gcd(Gc, Ac, Bc, ctx);
    if (!success)
        goto cleanup;

    _fmpz_mpoly_vec_divexact_mpoly(Ax->coeffs, Ax->length, Ac, ctx);
    _fmpz_mpoly_vec_divexact_mpoly(Bx->coeffs, Bx->length, Bc, ctx);

    success = fmpz_cmp(Ax->exps + 0, Bx->exps + 0) > 0 ?
                            fmpz_mpoly_univar_pseudo_gcd(Gx, Ax, Bx, ctx) :
                            fmpz_mpoly_univar_pseudo_gcd(Gx, Bx, Ax, ctx);
    if (!success)
        goto cleanup;

    if (fmpz_mpoly_gcd(t, Ax->coeffs + 0, Bx->coeffs + 0, ctx) &&
                                                                t->length == 1)
    {
        fmpz_mpoly_term_content(s, Gx->coeffs + 0, ctx);
        fmpz_mpoly_divexact(t, Gx->coeffs + 0, s, ctx);
        _fmpz_mpoly_vec_divexact_mpoly(Gx->coeffs, Gx->length, t, ctx);
    }
    else if (fmpz_mpoly_gcd(t, Ax->coeffs + Ax->length - 1,
                           Bx->coeffs + Bx->length - 1, ctx) && t->length == 1)
    {
        fmpz_mpoly_term_content(s, Gx->coeffs + Gx->length - 1, ctx);
        fmpz_mpoly_divexact(t, Gx->coeffs + Gx->length - 1, s, ctx);
        _fmpz_mpoly_vec_divexact_mpoly(Gx->coeffs, Gx->length, t, ctx);
    }

    success = _fmpz_mpoly_vec_content_mpoly(t, Gx->coeffs, Gx->length, ctx);
    if (!success)
        goto cleanup;

    _fmpz_mpoly_vec_divexact_mpoly(Gx->coeffs, Gx->length, t, ctx);
    _fmpz_mpoly_vec_mul_mpoly(Gx->coeffs, Gx->length, Gc, ctx);
    _fmpz_mpoly_from_univar(Gc, I->Gbits, Gx, var, ctx);

    if (Abar != NULL)
        fmpz_mpoly_divexact(s, A, Gc, ctx);

    if (Bbar != NULL)
        fmpz_mpoly_divexact(t, B, Gc, ctx);

    fmpz_mpoly_swap(G, Gc, ctx);

    if (Abar != NULL)
        fmpz_mpoly_swap(Abar, s, ctx);

    if (Bbar != NULL)
        fmpz_mpoly_swap(Bbar, t, ctx);

    success = 1;

cleanup:

    fmpz_mpoly_clear(Ac, ctx);
    fmpz_mpoly_clear(Bc, ctx);
    fmpz_mpoly_clear(Gc, ctx);
    fmpz_mpoly_clear(s, ctx);
    fmpz_mpoly_clear(t, ctx);
    fmpz_mpoly_univar_clear(Ax, ctx);
    fmpz_mpoly_univar_clear(Bx, ctx);
    fmpz_mpoly_univar_clear(Gx, ctx);

    return success;
}

/********************** Hit A and B with zippel ******************************/
static int _try_zippel(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const mpoly_gcd_info_t I,
    const fmpz_mpoly_ctx_t ctx)
{
    slong m = I->mvars;
    int success;
    flint_bitcnt_t wbits;
    flint_rand_t state;
    fmpz_mpoly_ctx_t lctx;
    fmpz_mpoly_t Al, Bl, Gl, Abarl, Bbarl;
    fmpz_mpoly_t Ac, Bc, Gc, Abarc, Bbarc;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);

    if (!(I->can_use & MPOLY_GCD_USE_ZIPPEL))
        return 0;

    FLINT_ASSERT(m > 1);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    flint_randinit(state);

    fmpz_mpoly_ctx_init(lctx, m, ORD_LEX);

    wbits = FLINT_MAX(A->bits, B->bits);

    fmpz_mpoly_init3(Al, 0, wbits, lctx);
    fmpz_mpoly_init3(Bl, 0, wbits, lctx);
    fmpz_mpoly_init3(Gl, 0, wbits, lctx);
    fmpz_mpoly_init3(Abarl, 0, wbits, lctx);
    fmpz_mpoly_init3(Bbarl, 0, wbits, lctx);
    fmpz_mpoly_init3(Ac, 0, wbits, lctx);
    fmpz_mpoly_init3(Bc, 0, wbits, lctx);
    fmpz_mpoly_init3(Gc, 0, wbits, lctx);
    fmpz_mpoly_init3(Abarc, 0, wbits, lctx);
    fmpz_mpoly_init3(Bbarc, 0, wbits, lctx);

    fmpz_mpoly_to_mpolyl_perm_deflate(Al, lctx, A, ctx,
                                      I->zippel_perm, I->Amin_exp, I->Gstride);
    fmpz_mpoly_to_mpolyl_perm_deflate(Bl, lctx, B, ctx,
                                      I->zippel_perm, I->Bmin_exp, I->Gstride);

    FLINT_ASSERT(Al->bits == wbits);
    FLINT_ASSERT(Bl->bits == wbits);
    FLINT_ASSERT(Al->length > 1);
    FLINT_ASSERT(Bl->length > 1);

    success = fmpz_mpolyl_content(Ac, Al, 1, lctx) &&
              fmpz_mpolyl_content(Bc, Bl, 1, lctx);
    if (!success)
        goto cleanup;

    fmpz_mpoly_repack_bits_inplace(Ac, wbits, lctx);
    fmpz_mpoly_repack_bits_inplace(Bc, wbits, lctx);

    success = _fmpz_mpoly_gcd_algo(Gc, Abar == NULL ? NULL : Abarc,
                                       Bbar == NULL ? NULL : Bbarc,
                                              Ac, Bc, lctx, MPOLY_GCD_USE_ALL);
    if (!success)
        goto cleanup;

    success = fmpz_mpoly_divides(Al, Al, Ac, lctx);
    FLINT_ASSERT(success);
    success = fmpz_mpoly_divides(Bl, Bl, Bc, lctx);
    FLINT_ASSERT(success);

    fmpz_mpoly_repack_bits_inplace(Al, wbits, lctx);
    fmpz_mpoly_repack_bits_inplace(Bl, wbits, lctx);

    success = fmpz_mpolyl_gcd_zippel(Gl, Abarl, Bbarl, Al, Bl, lctx, state);
    if (!success)
        goto cleanup;

    fmpz_mpoly_mul(Gl, Gl, Gc, lctx);
    fmpz_mpoly_repack_bits_inplace(Gl, wbits, lctx);
    fmpz_mpoly_from_mpolyl_perm_inflate(G, I->Gbits, ctx, Gl, lctx,
                                      I->zippel_perm, I->Gmin_exp, I->Gstride);

    if (Abar != NULL)
    {
        fmpz_mpoly_mul(Abarl, Abarl, Abarc, lctx);
        fmpz_mpoly_repack_bits_inplace(Abarl, wbits, lctx);
        fmpz_mpoly_from_mpolyl_perm_inflate(Abar, I->Abarbits, ctx, Abarl, lctx,
                                   I->zippel_perm, I->Abarmin_exp, I->Gstride);
    }

    if (Bbar != NULL)
    {
        fmpz_mpoly_mul(Bbarl, Bbarl, Bbarc, lctx);
        fmpz_mpoly_repack_bits_inplace(Bbarl, wbits, lctx);
        fmpz_mpoly_from_mpolyl_perm_inflate(Bbar, I->Bbarbits, ctx, Bbarl, lctx,
                                   I->zippel_perm, I->Bbarmin_exp, I->Gstride);
    }

    success = 1;

cleanup:

    fmpz_mpoly_clear(Al, lctx);
    fmpz_mpoly_clear(Bl, lctx);
    fmpz_mpoly_clear(Gl, lctx);
    fmpz_mpoly_clear(Abarl, lctx);
    fmpz_mpoly_clear(Bbarl, lctx);
    fmpz_mpoly_clear(Ac, lctx);
    fmpz_mpoly_clear(Bc, lctx);
    fmpz_mpoly_clear(Gc, lctx);
    fmpz_mpoly_clear(Abarc, lctx);
    fmpz_mpoly_clear(Bbarc, lctx);

    fmpz_mpoly_ctx_clear(lctx);

    flint_randclear(state);

    return success;
}


/************************ Hit A and B with bma *******************************/

static int _try_bma(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const mpoly_gcd_info_t I,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, k;
    slong m = I->mvars;
    int success;
    flint_bitcnt_t wbits;
    fmpz_mpoly_ctx_t lctx;
    fmpz_mpoly_t Al, Bl, Gl, Abarl, Bbarl;
    fmpz_mpoly_t Ac, Bc, Gc, Abarc, Bbarc, Gamma, lcAl, lcBl;
    slong max_deg;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (!(I->can_use & MPOLY_GCD_USE_ZIPPEL2))
        return 0;

    FLINT_ASSERT(m >= 3);

    fmpz_mpoly_ctx_init(lctx, m, ORD_LEX);

    max_deg = 0;
    for (i = 0; i < m; i++)
    {
        k = I->zippel2_perm[i];
        max_deg = FLINT_MAX(max_deg, I->Adeflate_deg[k]);
        max_deg = FLINT_MAX(max_deg, I->Bdeflate_deg[k]);
    }

    wbits = 1 + FLINT_BIT_COUNT(max_deg);
    wbits = FLINT_MAX(MPOLY_MIN_BITS, wbits);
    wbits = mpoly_fix_bits(wbits, lctx->minfo);
    FLINT_ASSERT(wbits <= FLINT_BITS);

    fmpz_mpoly_init3(Al, A->length, wbits, lctx);
    fmpz_mpoly_init3(Bl, B->length, wbits, lctx);
    fmpz_mpoly_init3(Gl, 0, wbits, lctx);
    fmpz_mpoly_init3(Abarl, 0, wbits, lctx);
    fmpz_mpoly_init3(Bbarl, 0, wbits, lctx);
    fmpz_mpoly_init3(Ac, 0, wbits, lctx);
    fmpz_mpoly_init3(Bc, 0, wbits, lctx);
    fmpz_mpoly_init3(Gc, 0, wbits, lctx);
    fmpz_mpoly_init3(Abarc, 0, wbits, lctx);
    fmpz_mpoly_init3(Bbarc, 0, wbits, lctx);
    fmpz_mpoly_init3(Gamma, 0, wbits, lctx);
    fmpz_mpoly_init3(lcAl, 0, wbits, lctx);
    fmpz_mpoly_init3(lcBl, 0, wbits, lctx);

    fmpz_mpoly_to_mpolyl_perm_deflate(Al, lctx, A, ctx,
                                     I->zippel2_perm, I->Amin_exp, I->Gstride);

    fmpz_mpoly_to_mpolyl_perm_deflate(Bl, lctx, B, ctx,
                                     I->zippel2_perm, I->Bmin_exp, I->Gstride);

    success = fmpz_mpolyl_content(Ac, Al, 2, lctx) &&
              fmpz_mpolyl_content(Bc, Bl, 2, lctx);
    if (!success)
        goto cleanup;

    success = _fmpz_mpoly_gcd_algo(Gc, Abar == NULL ? NULL : Abarc,
                                       Bbar == NULL ? NULL : Bbarc,
                                              Ac, Bc, lctx, MPOLY_GCD_USE_ALL);
    if (!success)
        goto cleanup;

    success = fmpz_mpoly_divides(Al, Al, Ac, lctx);
    FLINT_ASSERT(success);
    success = fmpz_mpoly_divides(Bl, Bl, Bc, lctx);
    FLINT_ASSERT(success);

    fmpz_mpoly_repack_bits_inplace(Al, wbits, lctx);
    fmpz_mpoly_repack_bits_inplace(Bl, wbits, lctx);

    fmpz_mpolyl_lead_coeff(lcAl, Al, 2, lctx);
    fmpz_mpolyl_lead_coeff(lcBl, Bl, 2, lctx);
    success = fmpz_mpoly_gcd(Gamma, lcAl, lcBl, lctx);
    if (!success)
        goto cleanup;

    fmpz_mpoly_repack_bits_inplace(Gamma, wbits, lctx);

    success = fmpz_mpolyl_gcd_zippel2(Gl, Abarl, Bbarl, Al, Bl, Gamma, lctx);
    if (!success)
        goto cleanup;

    fmpz_mpoly_mul(Gl, Gl, Gc, lctx);
    fmpz_mpoly_repack_bits_inplace(Gl, wbits, lctx);
    fmpz_mpoly_from_mpolyl_perm_inflate(G, I->Gbits, ctx, Gl, lctx,
                                     I->zippel2_perm, I->Gmin_exp, I->Gstride);

    if (Abar != NULL)
    {
        fmpz_mpoly_mul(Abarl, Abarl, Abarc, lctx);
        fmpz_mpoly_repack_bits_inplace(Abarl, wbits, lctx);
        fmpz_mpoly_from_mpolyl_perm_inflate(Abar, I->Abarbits, ctx, Abarl, lctx,
                                  I->zippel2_perm, I->Abarmin_exp, I->Gstride);
    }

    if (Bbar != NULL)
    {
        fmpz_mpoly_mul(Bbarl, Bbarl, Bbarc, lctx);
        fmpz_mpoly_repack_bits_inplace(Bbarl, wbits, lctx);
        fmpz_mpoly_from_mpolyl_perm_inflate(Bbar, I->Bbarbits, ctx, Bbarl, lctx,
                                  I->zippel2_perm, I->Bbarmin_exp, I->Gstride);
    }

    success = 1;

cleanup:

    fmpz_mpoly_clear(Al, lctx);
    fmpz_mpoly_clear(Bl, lctx);
    fmpz_mpoly_clear(Gl, lctx);
    fmpz_mpoly_clear(Abarl, lctx);
    fmpz_mpoly_clear(Bbarl, lctx);
    fmpz_mpoly_clear(Ac, lctx);
    fmpz_mpoly_clear(Bc, lctx);
    fmpz_mpoly_clear(Gc, lctx);
    fmpz_mpoly_clear(Abarc, lctx);
    fmpz_mpoly_clear(Bbarc, lctx);
    fmpz_mpoly_clear(Gamma, lctx);
    fmpz_mpoly_clear(lcAl, lctx);
    fmpz_mpoly_clear(lcBl, lctx);

    fmpz_mpoly_ctx_clear(lctx);

    return success;
}


/******************** Hit A and B with hensel lifting ************************/
static int _try_hensel(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const mpoly_gcd_info_t I,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, k;
    slong m = I->mvars;
    int success;
    flint_bitcnt_t wbits;
    fmpz_mpoly_ctx_t lctx;
    fmpz_mpoly_t Al, Bl, Gl, Abarl, Bbarl;
    fmpz_mpoly_t Ac, Bc, Gc, Abarc, Bbarc;
    slong max_deg;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (!(I->can_use & MPOLY_GCD_USE_HENSEL))
        return 0;

    FLINT_ASSERT(m >= 2);

    fmpz_mpoly_ctx_init(lctx, m, ORD_LEX);

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

    fmpz_mpoly_init3(Al, 0, wbits, lctx);
    fmpz_mpoly_init3(Bl, 0, wbits, lctx);
    fmpz_mpoly_init3(Gl, 0, wbits, lctx);
    fmpz_mpoly_init3(Abarl, 0, wbits, lctx);
    fmpz_mpoly_init3(Bbarl, 0, wbits, lctx);
    fmpz_mpoly_init3(Ac, 0, wbits, lctx);
    fmpz_mpoly_init3(Bc, 0, wbits, lctx);
    fmpz_mpoly_init3(Gc, 0, wbits, lctx);
    fmpz_mpoly_init3(Abarc, 0, wbits, lctx);
    fmpz_mpoly_init3(Bbarc, 0, wbits, lctx);

    fmpz_mpoly_to_mpolyl_perm_deflate(Al, lctx, A, ctx,
                                      I->hensel_perm, I->Amin_exp, I->Gstride);
    fmpz_mpoly_to_mpolyl_perm_deflate(Bl, lctx, B, ctx,
                                      I->hensel_perm, I->Bmin_exp, I->Gstride);

    success = fmpz_mpolyl_content(Ac, Al, 1, lctx) &&
              fmpz_mpolyl_content(Bc, Bl, 1, lctx);
    if (!success)
        goto cleanup;

    success = _fmpz_mpoly_gcd_algo(Gc, Abar == NULL ? NULL : Abarc,
                 Bbar == NULL ? NULL : Bbarc, Ac, Bc, lctx, MPOLY_GCD_USE_ALL);
    if (!success)
        goto cleanup;

    success = fmpz_mpoly_divides(Al, Al, Ac, lctx);
    FLINT_ASSERT(success);
    success = fmpz_mpoly_divides(Bl, Bl, Bc, lctx);
    FLINT_ASSERT(success);

    fmpz_mpoly_repack_bits_inplace(Al, wbits, lctx);
    fmpz_mpoly_repack_bits_inplace(Bl, wbits, lctx);

    max_deg = I->Gdeflate_deg_bound[I->hensel_perm[0]];
    success = fmpz_mpolyl_gcd_hensel(Gl, max_deg, Abarl, Bbarl, Al, Bl, lctx);
    if (!success)
        goto cleanup;

    fmpz_mpoly_mul(Gl, Gl, Gc, lctx);
    fmpz_mpoly_from_mpolyl_perm_inflate(G, I->Gbits, ctx, Gl, lctx,
                                      I->hensel_perm, I->Gmin_exp, I->Gstride);
    if (Abar != NULL)
    {
        fmpz_mpoly_mul(Abarl, Abarl, Abarc, lctx);
        fmpz_mpoly_from_mpolyl_perm_inflate(Abar, I->Abarbits, ctx, Abarl, lctx,
                                   I->hensel_perm, I->Abarmin_exp, I->Gstride);
    }

    if (Bbar != NULL)
    {
        fmpz_mpoly_mul(Bbarl, Bbarl, Bbarc, lctx);
        fmpz_mpoly_from_mpolyl_perm_inflate(Bbar, I->Bbarbits, ctx, Bbarl, lctx,
                                   I->hensel_perm, I->Bbarmin_exp, I->Gstride);
    }

    success = 1;

cleanup:

    fmpz_mpoly_clear(Al, lctx);
    fmpz_mpoly_clear(Bl, lctx);
    fmpz_mpoly_clear(Gl, lctx);
    fmpz_mpoly_clear(Abarl, lctx);
    fmpz_mpoly_clear(Bbarl, lctx);
    fmpz_mpoly_clear(Ac, lctx);
    fmpz_mpoly_clear(Bc, lctx);
    fmpz_mpoly_clear(Gc, lctx);
    fmpz_mpoly_clear(Abarc, lctx);
    fmpz_mpoly_clear(Bbarc, lctx);

    fmpz_mpoly_ctx_clear(lctx);

    return success;
}

/*********************** Hit A and B with brown ******************************/
static int _try_brown(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    mpoly_gcd_info_t I,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong m = I->mvars;
    flint_bitcnt_t wbits;
    fmpz_mpoly_ctx_t lctx;
    fmpz_mpoly_t Al, Bl, Gl, Abarl, Bbarl;

    if (!(I->can_use & MPOLY_GCD_USE_BROWN))
        return 0;

    FLINT_ASSERT(m >= 2);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    wbits = FLINT_MAX(A->bits, B->bits);

    fmpz_mpoly_ctx_init(lctx, m, ORD_LEX);
    fmpz_mpoly_init3(Al, 0, wbits, lctx);
    fmpz_mpoly_init3(Bl, 0, wbits, lctx);
    fmpz_mpoly_init3(Gl, 0, wbits, lctx);
    fmpz_mpoly_init3(Abarl, 0, wbits, lctx);
    fmpz_mpoly_init3(Bbarl, 0, wbits, lctx);

    fmpz_mpoly_to_mpolyl_perm_deflate(Al, lctx, A, ctx,
                                       I->brown_perm, I->Amin_exp, I->Gstride);
    fmpz_mpoly_to_mpolyl_perm_deflate(Bl, lctx, B, ctx,
                                       I->brown_perm, I->Bmin_exp, I->Gstride);
    FLINT_ASSERT(Al->bits == wbits);
    FLINT_ASSERT(Bl->bits == wbits);
    FLINT_ASSERT(Al->length > 1);
    FLINT_ASSERT(Bl->length > 1);

    success = fmpz_mpolyl_gcd_brown(Gl, Abarl, Bbarl, Al, Bl, lctx, I);
    if (!success)
        goto cleanup;

    fmpz_mpoly_from_mpolyl_perm_inflate(G, I->Gbits, ctx, Gl, lctx,
                                       I->brown_perm, I->Gmin_exp, I->Gstride);
    if (Abar != NULL)
        fmpz_mpoly_from_mpolyl_perm_inflate(Abar, I->Abarbits, ctx, Abarl, lctx,
                                    I->brown_perm, I->Abarmin_exp, I->Gstride);

    if (Bbar != NULL)
        fmpz_mpoly_from_mpolyl_perm_inflate(Bbar, I->Bbarbits, ctx, Bbarl, lctx,
                                    I->brown_perm, I->Bbarmin_exp, I->Gstride);

    success = 1;

cleanup:

    fmpz_mpoly_clear(Al, lctx);
    fmpz_mpoly_clear(Bl, lctx);
    fmpz_mpoly_clear(Gl, lctx);
    fmpz_mpoly_clear(Abarl, lctx);
    fmpz_mpoly_clear(Bbarl, lctx);
    fmpz_mpoly_ctx_clear(lctx);

    return success;
}

/*
    Both A and B have to be packed into bits <= FLINT_BITS.
    return is 1 for success, 0 for failure.
*/
static int _fmpz_mpoly_gcd_algo_small(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar, /* could be NULL */
    fmpz_mpoly_t Bbar, /* could be NULL */
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx,
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
    fmpz_mpoly_t T, Asave, Bsave;
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
    fmpz_mpoly_init(T, ctx);
    fmpz_mpoly_init(Asave, ctx);
    fmpz_mpoly_init(Bsave, ctx);
    fmpz_mpoly_set(Asave, A, ctx);
    fmpz_mpoly_set(Bsave, B, ctx);
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

    mpoly_gcd_info_stride(I->Gstride,
            A->exps, A->bits, A->length, I->Amax_exp, I->Amin_exp,
            B->exps, B->bits, B->length, I->Bmax_exp, I->Bmin_exp, ctx->minfo);

    for (j = 0; j < nvars; j++)
    {
        ulong t = I->Gstride[j];

        if (t == 0)
        {
            FLINT_ASSERT(  I->Amax_exp[j] == I->Amin_exp[j]
                        || I->Bmax_exp[j] == I->Bmin_exp[j]);
        }
        else
        {
            FLINT_ASSERT((I->Amax_exp[j] - I->Amin_exp[j]) % t == 0);
            FLINT_ASSERT((I->Bmax_exp[j] - I->Bmin_exp[j]) % t == 0);
        }

        I->Adeflate_deg[j] = t == 0 ? 0 : (I->Amax_exp[j] - I->Amin_exp[j])/t;
        I->Bdeflate_deg[j] = t == 0 ? 0 : (I->Bmax_exp[j] - I->Bmin_exp[j])/t;
        I->Gmin_exp[j] = FLINT_MIN(I->Amin_exp[j], I->Bmin_exp[j]);

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

    /*
        TODO: (1) the single-bit algo == MPOLY_GCD_USE_xxx is supposed to only
                  use that one algorithm xxx, and
              (2) multi-bit algo is supposed to choose among the present bits.
        (2) is not implemented completely.
    */

    if (algo == MPOLY_GCD_USE_HENSEL)
    {
        mpoly_gcd_info_measure_hensel(I, A->length, B->length, ctx->minfo);
        success = _try_hensel(G, Abar, Bbar, A, B, I, ctx);
        goto cleanup;
    }
    else if (algo == MPOLY_GCD_USE_PRS)
    {
        success = _try_prs(G, Abar, Bbar, A, B, I, ctx);
        goto cleanup;
    }
    else if (algo == MPOLY_GCD_USE_ZIPPEL)
    {
        mpoly_gcd_info_measure_zippel(I, A->length, B->length, ctx->minfo);
        success = _try_zippel(G, Abar, Bbar, A, B, I, ctx);
        goto cleanup;
    }

    mpoly_gcd_info_measure_brown(I, A->length, B->length, ctx->minfo);
    mpoly_gcd_info_measure_bma(I, A->length, B->length, ctx->minfo);

    if (I->mvars < 3)
    {
        /* TODO: bivariate heuristic here */

        if (_try_brown(G, Abar, Bbar, A, B, I, ctx))
            goto successful;
    }
    else if (algo == MPOLY_GCD_USE_BROWN)
    {
        success = _try_brown(G, Abar, Bbar, A, B, I, ctx);
        goto cleanup;
    }
    else if (algo == MPOLY_GCD_USE_ZIPPEL2)
    {
        success = _try_bma(G, Abar, Bbar, A, B, I, ctx);
        goto cleanup;
    }
    else if ((I->can_use & MPOLY_GCD_USE_BROWN) &&
             (I->can_use & MPOLY_GCD_USE_ZIPPEL2) &&
             I->zippel2_time < I->brown_time &&
             (I->mvars*(I->Adensity + I->Bdensity) < 1 ||
              I->zippel2_time < 0.01*I->brown_time))
    {
        if (_try_bma(G, Abar, Bbar, A, B, I, ctx))
            goto successful;

        if (_try_brown(G, Abar, Bbar, A, B, I, ctx))
            goto successful;
    }
    else
    {
        if (_try_brown(G, Abar, Bbar, A, B, I, ctx))
            goto successful;

        if (_try_bma(G, Abar, Bbar, A, B, I, ctx))
            goto successful;
    }

    mpoly_gcd_info_measure_zippel(I, A->length, B->length, ctx->minfo);

    if (_try_zippel(G, Abar, Bbar, A, B, I, ctx))
        goto successful;

    success = 0;
    goto cleanup;

successful:

    success = 1;

cleanup:

    mpoly_gcd_info_clear(I);

    if (success)
    {
        FLINT_ASSERT(G->length > 0);
        if (fmpz_sgn(G->coeffs + 0) < 0)
        {
            _fmpz_vec_neg(G->coeffs, G->coeffs, G->length);

            if (Abar != NULL)
                _fmpz_vec_neg(Abar->coeffs, Abar->coeffs, Abar->length);

            if (Bbar != NULL)
                _fmpz_vec_neg(Bbar->coeffs, Bbar->coeffs, Bbar->length);
        }

        FLINT_ASSERT(fmpz_mpoly_divides(T, Asave, G, ctx));
        FLINT_ASSERT(Abar == NULL || fmpz_mpoly_equal(T, Abar, ctx));
            
        FLINT_ASSERT(fmpz_mpoly_divides(T, Bsave, G, ctx));
        FLINT_ASSERT(Bbar == NULL || fmpz_mpoly_equal(T, Bbar, ctx));
    }

#if FLINT_WANT_ASSERT
    fmpz_mpoly_clear(T, ctx);
    fmpz_mpoly_clear(Asave, ctx);
    fmpz_mpoly_clear(Bsave, ctx);
#endif

    return success;
}


static int _fmpz_mpoly_gcd_algo_large(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong k;
    fmpz * Ashift, * Astride;
    fmpz * Bshift, * Bstride;
    fmpz * Gshift, * Gstride;
    fmpz_mpoly_t Anew, Bnew;
    const fmpz_mpoly_struct * Ause, * Buse;

    if (A->length == 1)
        return _do_monomial_gcd(G, Bbar, Abar, B, A, ctx);

    if (B->length == 1)
        return _do_monomial_gcd(G, Abar, Bbar, A, B, ctx);

    if (_try_monomial_cofactors(G, Abar, Bbar, A, B, ctx))
        return 1;

    fmpz_mpoly_init(Anew, ctx);
    fmpz_mpoly_init(Bnew, ctx);

    Ause = A;
    if (A->bits > FLINT_BITS)
    {
        if (!fmpz_mpoly_repack_bits(Anew, A, FLINT_BITS, ctx))
            goto could_not_repack;
        Ause = Anew;
    }

    Buse = B;
    if (B->bits > FLINT_BITS)
    {
        if (!fmpz_mpoly_repack_bits(Bnew, B, FLINT_BITS, ctx))
            goto could_not_repack;
        Buse = Bnew;
    }

    success = _fmpz_mpoly_gcd_algo_small(G, Abar, Bbar, Ause, Buse, ctx, algo);
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

    fmpz_mpoly_deflation(Ashift, Astride, A, ctx);
    fmpz_mpoly_deflation(Bshift, Bstride, B, ctx);
    _fmpz_vec_min(Gshift, Ashift, Bshift, ctx->minfo->nvars);
    for (k = 0; k < ctx->minfo->nvars; k++)
        fmpz_gcd(Gstride + k, Astride + k, Bstride + k);

    fmpz_mpoly_deflate(Anew, A, Ashift, Gstride, ctx);
    if (Anew->bits > FLINT_BITS)
    {
        success = fmpz_mpoly_repack_bits(Anew, Anew, FLINT_BITS, ctx);
        if (!success)
            goto deflate_cleanup;
    }

    fmpz_mpoly_deflate(Bnew, B, Bshift, Gstride, ctx);
    if (Bnew->bits > FLINT_BITS)
    {
        success = fmpz_mpoly_repack_bits(Bnew, Bnew, FLINT_BITS, ctx);
        if (!success)
            goto deflate_cleanup;
    }

    success = _fmpz_mpoly_gcd_algo_small(G, Abar, Bbar, Anew, Bnew, ctx, algo);
    if (!success)
        goto deflate_cleanup;

    for (k = 0; k < ctx->minfo->nvars; k++)
    {
        fmpz_sub(Ashift + k, Ashift + k, Gshift + k);
        fmpz_sub(Bshift + k, Bshift + k, Gshift + k);
        FLINT_ASSERT(fmpz_sgn(Ashift + k) >= 0);
        FLINT_ASSERT(fmpz_sgn(Bshift + k) >= 0);
    }

    fmpz_mpoly_inflate(G, G, Gshift, Gstride, ctx);

    if (Abar != NULL)
        fmpz_mpoly_inflate(Abar, Abar, Ashift, Gstride, ctx);

    if (Bbar != NULL)
        fmpz_mpoly_inflate(Bbar, Bbar, Bshift, Gstride, ctx);

    FLINT_ASSERT(G->length > 0);
    if (fmpz_sgn(G->coeffs + 0) < 0)
    {
        _fmpz_vec_neg(G->coeffs, G->coeffs, G->length);

        if (Abar != NULL)
            _fmpz_vec_neg(Abar->coeffs, Abar->coeffs, Abar->length);

        if (Bbar != NULL)
            _fmpz_vec_neg(Bbar->coeffs, Bbar->coeffs, Bbar->length);
    }

deflate_cleanup:

    _fmpz_vec_clear(Ashift, ctx->minfo->nvars);
    _fmpz_vec_clear(Astride, ctx->minfo->nvars);
    _fmpz_vec_clear(Bshift, ctx->minfo->nvars);
    _fmpz_vec_clear(Bstride, ctx->minfo->nvars);
    _fmpz_vec_clear(Gshift, ctx->minfo->nvars);
    _fmpz_vec_clear(Gstride, ctx->minfo->nvars);

cleanup:

    fmpz_mpoly_clear(Anew, ctx);
    fmpz_mpoly_clear(Bnew, ctx);

    return success;
}


int _fmpz_mpoly_gcd_algo(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,  /* could be NULL */
    fmpz_mpoly_t Bbar,  /* could be NULL */
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx,
    unsigned int algo)
{
    FLINT_ASSERT(!fmpz_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!fmpz_mpoly_is_zero(B, ctx));

    if (A->bits <= FLINT_BITS && B->bits <= FLINT_BITS)
        return _fmpz_mpoly_gcd_algo_small(G, Abar, Bbar, A, B, ctx, algo);
    else
        return _fmpz_mpoly_gcd_algo_large(G, Abar, Bbar, A, B, ctx, algo);
}

