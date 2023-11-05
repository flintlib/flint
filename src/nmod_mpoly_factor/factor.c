/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "long_extras.h"
#include "nmod_poly_factor.h"
#include "nmod_mpoly_factor.h"

static slong _deflate(
    nmod_mpoly_t A,
    slong tot_deg,
    const ulong * strides,
    const slong * perm,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong nvars = ctx->minfo->nvars;
    flint_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    ulong * texps, * sexps;
    TMP_INIT;

    for (j = 0; j < nvars; j++)
    {
        if (strides[j] != 1 || perm[j] != j)
            goto do_it;
    }

    return tot_deg;

do_it:

    TMP_START;

    texps = (ulong *) TMP_ALLOC(2*nvars*sizeof(ulong));
    sexps = texps + nvars;

    tot_deg = 1;
    for (i = 0; i < A->length; i++)
    {
        slong this_deg = 0;

        mpoly_get_monomial_ui(texps, A->exps + N*i, bits, ctx->minfo);

        for (j = 0; j < nvars; j++)
        {
            FLINT_ASSERT(0 == texps[j] % strides[j]);
            texps[j] = texps[j]/strides[j];
        }

        for (j = 0; j < nvars; j++)
        {
            sexps[j] = texps[perm[j]];
            this_deg += sexps[j];
        }

        tot_deg = FLINT_MAX(tot_deg, this_deg);

        mpoly_set_monomial_ui(A->exps + N*i, sexps, bits, ctx->minfo);
    }

    TMP_END;

    nmod_mpoly_sort_terms(A, ctx);
    nmod_mpoly_make_monic(A, A, ctx);

    return tot_deg;
}


static void _inflate(
    nmod_mpoly_t A,
    flint_bitcnt_t bits,
    const ulong * strides,
    const slong * perm,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong nvars = ctx->minfo->nvars;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    ulong * texps, * sexps;
    TMP_INIT;

    for (j = 0; j < nvars; j++)
    {
        if (strides[j] != 1 || perm[j] != j)
            goto do_it;
    }

    return;

do_it:

    nmod_mpoly_repack_bits_inplace(A, bits, ctx);

    TMP_START;

    texps = (ulong *) TMP_ALLOC(2*nvars*sizeof(ulong));
    sexps = texps + nvars;

    for (i = 0; i < A->length; i++)
    {
        mpoly_get_monomial_ui(sexps, A->exps + N*i, bits, ctx->minfo);

        for (j = 0; j < nvars; j++)
            texps[perm[j]] = sexps[j];

        for (j = 0; j < nvars; j++)
            texps[j] = texps[j]*strides[j];

        mpoly_set_monomial_ui(A->exps + N*i, texps, bits, ctx->minfo);
    }

    TMP_END;

    nmod_mpoly_sort_terms(A, ctx);
    nmod_mpoly_make_monic(A, A, ctx);

    return;
}


/* A has degree 2 wrt gen(0) */
static int _apply_quadratic(
    nmod_mpolyv_t Af,
    nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, shift, off, N;
    flint_bitcnt_t bits = A->bits;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    nmod_mpoly_t a_mock, b_mock, c_mock;
    nmod_mpoly_t t0, t1, t2, t3;

    FLINT_ASSERT(A->length > 1 || A->coeffs[0] != 0);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    nmod_mpoly_init(t0, ctx);
    nmod_mpoly_init(t1, ctx);
    nmod_mpoly_init(t2, ctx);
    nmod_mpoly_init(t3, ctx);

    mpoly_gen_offset_shift_sp(&off, &shift, 0, bits, ctx->minfo);
    N = mpoly_words_per_exp_sp(bits, ctx->minfo);

    i = 0;
    a_mock->exps = A->exps + N*i;
    a_mock->coeffs = A->coeffs + i;
    a_mock->bits = bits;
    while (i < A->length && (mask & (A->exps[N*i + off] >> shift)) == 2)
        i++;
    a_mock->length = i;
    a_mock->coeffs_alloc = a_mock->length;
    a_mock->exps_alloc = N*a_mock->length;

    b_mock->exps = A->exps + N*i;
    b_mock->coeffs = A->coeffs + i;
    b_mock->bits = bits;
    while (i < A->length && (mask & (A->exps[N*i + off] >> shift)) == 1)
        i++;
    b_mock->length = i - a_mock->length;
    b_mock->coeffs_alloc = b_mock->length;
    b_mock->exps_alloc = N*b_mock->length;

    c_mock->exps = A->exps + N*i;
    c_mock->coeffs = A->coeffs + i;
    c_mock->bits = bits;
    c_mock->length = A->length - i;
    c_mock->coeffs_alloc = c_mock->length;
    c_mock->exps_alloc = N*c_mock->length;

    FLINT_ASSERT(a_mock->length > 0);
    FLINT_ASSERT(c_mock->length > 0);

    nmod_mpoly_mul(t1, a_mock, c_mock, ctx);
    nmod_mpoly_neg(t1, t1, ctx);
    if (!nmod_mpoly_quadratic_root(t2, b_mock, t1, ctx))
    {
        nmod_mpolyv_fit_length(Af, 1, ctx);
        Af->length = 1;
        nmod_mpoly_swap(Af->coeffs + 0, A, ctx);
        success = 1;
        goto cleanup;
    }
    nmod_mpoly_neg(t2, t2, ctx);

    success = nmod_mpoly_gcd_cofactors(t0, t1, t2, a_mock, t2, ctx);
    if (!success)
        goto cleanup;

    nmod_mpoly_divides(t3, c_mock, t2, ctx);

    nmod_mpolyv_fit_length(Af, 2, ctx);
    Af->length = 2;
    nmod_mpoly_add(Af->coeffs + 0, t1, t2, ctx);
    nmod_mpoly_add(Af->coeffs + 1, t0, t3, ctx);

    success = 1;

cleanup:

    nmod_mpoly_clear(t0, ctx);
    nmod_mpoly_clear(t1, ctx);
    nmod_mpoly_clear(t2, ctx);
    nmod_mpoly_clear(t3, ctx);

    return success;
}

/*
    The property "sep" used here is that of the returned factors of
    _nmod_mpoly_factor_separable with sep = 1, namely:
        (1) monic
        (2) primitive wrt each variable
        (3) for all i, derivative(A, gen(i)) = 0, or
                       gcd(A, derivative(A, gen(i))) = 1
        (4) there is at least one i for which derivative(A, gen(i)) != 0

    Input A is sep and compressed.

    return 1 for success, 0 for failure
*/
static int _factor_irred_compressed(
    nmod_mpolyv_t Af,
    nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong i, j, tot_deg;
    slong nvars = ctx->minfo->nvars;
    slong * perm;
    ulong * strides, * texps;
    flint_bitcnt_t Abits;
    flint_rand_t state;
#ifdef FLINT_WANT_ASSERT
    nmod_mpoly_t Aorg;

    nmod_mpoly_init(Aorg, ctx);
    nmod_mpoly_set(Aorg, A, ctx);
#endif

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->coeffs[0] == 1);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    if (A->length < 2)
    {
        FLINT_ASSERT(A->length == 1);
        FLINT_ASSERT(!nmod_mpoly_is_ui(A, ctx));
        nmod_mpolyv_fit_length(Af, 1, ctx);
        nmod_mpoly_swap(Af->coeffs + 0, A, ctx);
        Af->length = 1;
        return 1;
    }

    if (A->bits > FLINT_BITS &&
        !nmod_mpoly_repack_bits_inplace(A, FLINT_BITS, ctx))
    {
        return 0;
    }

    Abits = A->bits;

    flint_randinit(state);

    strides = FLINT_ARRAY_ALLOC(2*nvars, ulong);
    texps = strides + nvars;

    perm = FLINT_ARRAY_ALLOC(nvars, slong);

    /* fill perm with id, and fill in strides */
    {
        ulong ppowt, ppow = ctx->mod.n;
        slong N = mpoly_words_per_exp_sp(Abits, ctx->minfo);

        while (!n_mul_checked(&ppowt, ppow, ctx->mod.n))
            ppow = ppowt;

        for (j = 0; j < nvars; j++)
        {
            strides[j] = ppow;
            perm[j] = j;
        }

        tot_deg = 1;
        for (i = 0; i < A->length; i++)
        {
            slong this_deg = 0;
	        mpoly_get_monomial_ui(texps, A->exps + N*i, Abits, ctx->minfo);
            for (j = 0; j < nvars; j++)
            {
                if (z_add_checked(&this_deg, this_deg, texps[j]))
                {
                    success = 0;
                    goto cleanup;
                }
                strides[j] = n_gcd(strides[j], texps[j]);
            }

            tot_deg = FLINT_MAX(tot_deg, this_deg);
        }
    }

    /* find permutation with gcd(A, derivative(A, gen(perm[0]))) = 1 */
    for (i = 0; i < nvars; i++)
    {
        if (strides[i] == 1)
        {
            FLINT_SWAP(slong, perm[0], perm[i]);
            break;
        }
    }

    if (nvars < 2)
    {
        nmod_poly_t Au;
        nmod_poly_factor_t Auf;

        FLINT_ASSERT(nvars == 1);

        nmod_poly_init_mod(Au, ctx->mod);
        nmod_poly_factor_init(Auf);

        FLINT_ASSERT(nmod_mpoly_is_nmod_poly(A, perm[0], ctx));
        success = nmod_mpoly_get_nmod_poly(Au, A, perm[0], ctx);
        FLINT_ASSERT(success);
        nmod_poly_factor(Auf, Au);

        nmod_mpolyv_fit_length(Af, Auf->num, ctx);
        Af->length = Auf->num;
        for (i = 0; i < Auf->num; i++)
        {
            FLINT_ASSERT(Auf->exp[i] == 1);
            _nmod_mpoly_set_nmod_poly(Af->coeffs + i, Abits,
                             Auf->p[i].coeffs, Auf->p[i].length, perm[0], ctx);
        }

        nmod_poly_clear(Au);
        nmod_poly_factor_clear(Auf);

        success = 1;
    }
    else if (nvars == 2)
    {
        n_poly_t c;
        n_bpoly_t Ab;
        n_tpoly_t Abf;

        n_poly_init(c);
        n_bpoly_init(Ab);
        n_tpoly_init(Abf);

        nmod_mpoly_get_bpoly(Ab, A, perm[0], perm[1], ctx);
        success = n_bpoly_mod_factor_smprime(c, Abf, Ab, 1, ctx->mod);
        if (!success)
        {
            nmod_mpoly_get_bpoly(Ab, A, perm[0], perm[1], ctx);
            n_bpoly_mod_factor_lgprime(c, Abf, Ab, ctx->mod);
        }

        FLINT_ASSERT(n_poly_degree(c) == 0);

        nmod_mpolyv_fit_length(Af, Abf->length, ctx);
        Af->length = Abf->length;
        for (i = 0; i < Abf->length; i++)
        {
            nmod_mpoly_set_bpoly(Af->coeffs + i, Abits, Abf->coeffs + i,
                                                        perm[0], perm[1], ctx);
            nmod_mpoly_make_monic(Af->coeffs + i, Af->coeffs + i, ctx);
        }

        n_poly_clear(c);
        n_bpoly_clear(Ab);
        n_tpoly_clear(Abf);

        success = 1;
    }
    else
    {
        slong Adeg0;
        nmod_mpoly_t lcA;
        nmod_mpoly_factor_t lcAf;

        nmod_mpoly_init(lcA, ctx);
        nmod_mpoly_factor_init(lcAf, ctx);

        tot_deg = _deflate(A, tot_deg, strides, perm, ctx);

        #ifdef FLINT_WANT_ASSERT
        {
            nmod_mpoly_t g;
            nmod_mpoly_init(g, ctx);
            nmod_mpoly_derivative(g, A, 0, ctx);
            FLINT_ASSERT(nmod_mpoly_gcd(g, g, A, ctx));
            FLINT_ASSERT(nmod_mpoly_is_one(g, ctx));
            nmod_mpoly_clear(g, ctx);
        }
        #endif

        Adeg0 = nmod_mpoly_degree_si(A, 0, ctx);
        if (Adeg0 == 1)
        {
            nmod_mpolyv_fit_length(Af, 1, ctx);
            Af->length = 1;
            nmod_mpoly_swap(Af->coeffs + 0, A, ctx);
            success = 1;
            goto cleanup_inflate;
        }
        else if (Adeg0 == 2)
        {
            success = _apply_quadratic(Af, A, ctx);
            goto cleanup_inflate;
        }

        success = 0;

        if (!(algo & (MPOLY_FACTOR_USE_WANG | MPOLY_FACTOR_USE_ZIP)))
            goto try_zassenhaus;

        /* TODO lcc_kaltofen */
        _nmod_mpoly_get_lead0(lcA, A, ctx);
        if (!nmod_mpoly_factor(lcAf, lcA, ctx))
            goto try_zassenhaus;

        if (!(algo & MPOLY_FACTOR_USE_ZIP))
        {
            if (success == 0)
                success = nmod_mpoly_factor_irred_smprime_wang(
                                                Af, A, lcAf, lcA, ctx, state);
            if (success == 0)
                success = nmod_mpoly_factor_irred_medprime_wang(
                                                 Af, A, lcAf, lcA, ctx, state);
            if (success == 0)
                success = nmod_mpoly_factor_irred_lgprime_wang(
                                                 Af, A, lcAf, lcA, ctx, state);
        }
        else if (!(algo & MPOLY_FACTOR_USE_WANG))
        {
            if (success == 0)
                success = nmod_mpoly_factor_irred_smprime_zippel(
                                                 Af, A, lcAf, lcA, ctx, state);
            if (success == 0)
                success = nmod_mpoly_factor_irred_medprime_zippel(
                                                 Af, A, lcAf, lcA, ctx, state);
            if (success == 0)
                success = nmod_mpoly_factor_irred_lgprime_zippel(
                                                 Af, A, lcAf, lcA, ctx, state);
        }
        else
        {
            double tdensity;
            fmpz_t x;
            fmpz_init(x);
            fmpz_bin_uiui(x, (ulong)tot_deg + nvars, nvars);
            tdensity = A->length/fmpz_get_d(x);
            fmpz_clear(x);

            if (tdensity > 0.005)
            {
                if (success == 0)
                    success = nmod_mpoly_factor_irred_smprime_wang(
                                                 Af, A, lcAf, lcA, ctx, state);
                if (success == 0)
                    success = nmod_mpoly_factor_irred_medprime_wang(
                                                 Af, A, lcAf, lcA, ctx, state);
                if (success == 0)
                    success = nmod_mpoly_factor_irred_smprime_zippel(
                                                 Af, A, lcAf, lcA, ctx, state);
                if (success == 0)
                    success = nmod_mpoly_factor_irred_medprime_zippel(
                                                 Af, A, lcAf, lcA, ctx, state);
            }
            else
            {
                if (success == 0)
                    success = nmod_mpoly_factor_irred_smprime_zippel(
                                                 Af, A, lcAf, lcA, ctx, state);
                if (success == 0)
                    success = nmod_mpoly_factor_irred_medprime_zippel(
                                                 Af, A, lcAf, lcA, ctx, state);
                if (success == 0)
                    success = nmod_mpoly_factor_irred_smprime_wang(
                                                 Af, A, lcAf, lcA, ctx, state);
                if (success == 0)
                    success = nmod_mpoly_factor_irred_medprime_wang(
                                                 Af, A, lcAf, lcA, ctx, state);
            }

            if (tdensity > 0.001)
            {
                if (success == 0)
                    success = nmod_mpoly_factor_irred_lgprime_wang(
                                                 Af, A, lcAf, lcA, ctx, state);
                if (success == 0)
                    success = nmod_mpoly_factor_irred_lgprime_zippel(
                                                 Af, A, lcAf, lcA, ctx, state);
            }
            else
            {
                if (success == 0)
                    success = nmod_mpoly_factor_irred_lgprime_zippel(
                                                 Af, A, lcAf, lcA, ctx, state);
                if (success == 0)
                    success = nmod_mpoly_factor_irred_lgprime_wang(
                                                 Af, A, lcAf, lcA, ctx, state);
            }
        }

    try_zassenhaus:

        if (algo & MPOLY_FACTOR_USE_ZAS)
        {
            if (success == 0)
    		    success = nmod_mpoly_factor_irred_smprime_zassenhaus(
                                                            Af, A, ctx, state);
            if (success == 0)
        		success = nmod_mpoly_factor_irred_medprime_zassenhaus(
                                                            Af, A, ctx, state);
            if (success == 0)
        		success = nmod_mpoly_factor_irred_lgprime_zassenhaus(
                                                            Af, A, ctx, state);
        }

    cleanup_inflate:

        success = (success > 0);
		if (success)
        {
		    for (i = 0; i < Af->length; i++)
                _inflate(Af->coeffs + i, Abits, strides, perm, ctx);
        }

	    nmod_mpoly_clear(lcA, ctx);
        nmod_mpoly_factor_clear(lcAf, ctx);
    }

cleanup:

    flint_randclear(state);
    flint_free(strides);
    flint_free(perm);

#ifdef FLINT_WANT_ASSERT
    if (success)
    {
        nmod_mpoly_t prod;
        nmod_mpoly_init(prod, ctx);
        nmod_mpoly_one(prod, ctx);
        for (i = 0; i < Af->length; i++)
            nmod_mpoly_mul(prod, prod, Af->coeffs + i, ctx);

        FLINT_ASSERT(nmod_mpoly_equal(prod, Aorg, ctx));
        nmod_mpoly_clear(prod, ctx);
        nmod_mpoly_clear(Aorg, ctx);
    }
#endif

    FLINT_ASSERT(success == 0 || success == 1);
    return success;
}


/*
    f is already squarefree
    make the factors in f have the sep property
*/
static int _refine_sep(
    nmod_mpolyv_t f,
    const nmod_mpoly_ctx_t ctx,
    nmod_mpolyv_t g)    /* temp */
{
    int success;
    slong v, i;
    nmod_mpoly_struct * t;
    nmod_mpoly_univar_t u;

    nmod_mpoly_univar_init(u, ctx);

    /* first make primitive */
    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        g->length = 0;
        for (i = 0; i < f->length; i++)
        {
            nmod_mpoly_to_univar(u, f->coeffs + i, v, ctx);
            FLINT_ASSERT(u->length > 0);
            FLINT_ASSERT(fmpz_is_zero(u->exps + u->length - 1));

            nmod_mpolyv_fit_length(g, g->length + 2, ctx);
            success = _nmod_mpoly_vec_content_mpoly(g->coeffs + g->length,
                                                    u->coeffs, u->length, ctx);
            if (!success)
                goto cleanup;

            if (nmod_mpoly_is_ui(g->coeffs + g->length, ctx))
            {
                nmod_mpoly_swap(g->coeffs + g->length, f->coeffs + i, ctx);
                g->length++;
            }
            else
            {
                success = nmod_mpoly_divides(g->coeffs + g->length + 1,
                                    f->coeffs + i, g->coeffs + g->length, ctx);
                FLINT_ASSERT(success);

                if (nmod_mpoly_is_ui(g->coeffs + g->length + 1, ctx))
                    g->length += 1;
                else
                    g->length += 2;
            }
        }

        nmod_mpolyv_swap(f, g, ctx);
    }

    /* now make separable/derivative zero wrt each variable */
    nmod_mpolyv_fit_length(g, 1, ctx);
    t = g->coeffs + 0;
    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        i = 0;
        while (i < f->length)
        {
            nmod_mpoly_derivative(t, f->coeffs + i, v, ctx);
            if (nmod_mpoly_is_zero(t, ctx))
            {
                /* f[i] has zero derivative */
                i++;
                continue;
            }

            nmod_mpolyv_fit_length(f, f->length + 1, ctx);

            success = nmod_mpoly_gcd_cofactors(f->coeffs + f->length,
                                      f->coeffs + i, t, f->coeffs + i, t, ctx);
            if (!success)
                goto cleanup;

            if (nmod_mpoly_is_ui(f->coeffs + f->length, ctx))
            {
                /* f[i] is comprime with its derivative */
                i++;
            }
            else
            {
                /* f[i] and f[end] at least got smaller */
                f->length++;
            }
        }
    }

    success = 1;

cleanup:

    nmod_mpoly_univar_clear(u, ctx);

    return 1;
}


/*
    A is sep.

    return 1 for success, 0 for failure
*/
static int _factor_irred(
    nmod_mpolyv_t Af,
    nmod_mpoly_t A,
    const nmod_mpoly_ctx_t Actx,
    unsigned int algo)
{
    int success;
    slong i, j;
    flint_bitcnt_t Abits;
    mpoly_compression_t M;
#ifdef FLINT_WANT_ASSERT
    nmod_mpoly_t Aorg;

    nmod_mpoly_init(Aorg, Actx);
    nmod_mpoly_set(Aorg, A, Actx);
#endif

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->coeffs[0] == 1);

    if (A->length < 2)
    {
        FLINT_ASSERT(A->length == 1);
        FLINT_ASSERT(!nmod_mpoly_is_ui(A, Actx));
        nmod_mpolyv_fit_length(Af, 1, Actx);
        Af->length = 1;
        nmod_mpoly_swap(Af->coeffs + 0, A, Actx);
        success = 1;
        goto cleanup_less;
    }

    if (A->bits > FLINT_BITS &&
        !nmod_mpoly_repack_bits_inplace(A, FLINT_BITS, Actx))
    {
        success = 0;
        goto cleanup_less;
    }

    Abits = A->bits;

    mpoly_compression_init(M);
    mpoly_compression_set(M, A->exps, A->bits, A->length, Actx->minfo);

    if (M->is_irred)
    {
        nmod_mpolyv_fit_length(Af, 1, Actx);
        Af->length = 1;
        nmod_mpoly_swap(Af->coeffs + 0, A, Actx);
        success = 1;
    }
    else if (M->is_trivial)
    {
        success = _factor_irred_compressed(Af, A, Actx, algo);
    }
    else
    {
        nmod_mpoly_ctx_t Lctx;
        nmod_mpolyv_t Lf, Lft, Lfs;

        nmod_mpoly_ctx_init(Lctx, M->mvars, ORD_LEX, Actx->mod.n);
        nmod_mpolyv_init(Lf, Lctx);
        nmod_mpolyv_init(Lft, Lctx);
        nmod_mpolyv_init(Lfs, Lctx);

        nmod_mpolyv_fit_length(Lft, 1, Lctx);
        Lft->length = 1;
        nmod_mpoly_compression_do(Lft->coeffs + 0, Lctx, A->coeffs, A->length, M);

        _refine_sep(Lft, Lctx, Lf);

        if (Lft->length == 1)
        {
            success = _factor_irred_compressed(Lf, Lft->coeffs + 0, Lctx, algo);
        }
        else
        {
            success = 1;
            Lf->length = 0;
            for (i = 0; i < Lft->length; i++)
            {
                success = _factor_irred(Lfs, Lft->coeffs + i, Lctx, algo);
                if (!success)
                    break;

                nmod_mpolyv_fit_length(Lf, Lf->length + Lfs->length, Lctx);
                for (j = 0; j < Lfs->length; j++)
                {
                    nmod_mpoly_swap(Lf->coeffs + Lf->length, Lfs->coeffs + j, Lctx);
                    Lf->length++;
                }
            }
        }

        if (success)
        {
            nmod_mpolyv_fit_length(Af, Lf->length, Actx);
            Af->length = Lf->length;
            for (i = 0; i < Lf->length; i++)
            {
                nmod_mpoly_compression_undo(Af->coeffs + i, Abits, Actx,
                                                      Lf->coeffs + i, Lctx, M);
            }
        }

        nmod_mpolyv_clear(Lf, Lctx);
        nmod_mpolyv_clear(Lft, Lctx);
        nmod_mpolyv_clear(Lfs, Lctx);
        nmod_mpoly_ctx_clear(Lctx);
    }

    mpoly_compression_clear(M);

cleanup_less:

#ifdef FLINT_WANT_ASSERT
    if (success)
    {
        nmod_mpoly_t prod;
        nmod_mpoly_init(prod, Actx);
        nmod_mpoly_one(prod, Actx);
        for (i = 0; i < Af->length; i++)
            nmod_mpoly_mul(prod, prod, Af->coeffs + i, Actx);
        FLINT_ASSERT(nmod_mpoly_equal(prod, Aorg, Actx));
        nmod_mpoly_clear(prod, Actx);
        nmod_mpoly_clear(Aorg, Actx);
    }
#endif

    FLINT_ASSERT(success == 0 || success == 1);
    return success;
}


/*
    Assume each factor in f is sep.
    Replace f by an irreducible factorization.
*/
int nmod_mpoly_factor_irred(
    nmod_mpoly_factor_t f,
    const nmod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong i, j;
    nmod_mpolyv_t t;
    nmod_mpoly_factor_t g;

    nmod_mpolyv_init(t, ctx);
    nmod_mpoly_factor_init(g, ctx);

    g->constant = f->constant;
    g->num = 0;
    for (j = 0; j < f->num; j++)
    {
        success = _factor_irred(t, f->poly + j, ctx, algo);
        if (!success)
            goto cleanup;

        nmod_mpoly_factor_fit_length(g, g->num + t->length, ctx);
        for (i = 0; i < t->length; i++)
        {
            fmpz_set(g->exp + g->num, f->exp + j);
            nmod_mpoly_swap(g->poly + g->num, t->coeffs + i, ctx);
            g->num++;
        }
    }
    nmod_mpoly_factor_swap(f, g, ctx);

    success = 1;

cleanup:

    nmod_mpolyv_clear(t, ctx);
    nmod_mpoly_factor_clear(g, ctx);

    return success;
}


/*
    append factor(f)^e to g
    assuming f is compressed and content free
*/
static int _compressed_content_to_irred(
    nmod_mpoly_factor_t g,
    nmod_mpoly_t f,
    const fmpz_t e,
    const nmod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong j, k;
    nmod_mpoly_factor_t h;
    nmod_mpolyv_t v;

    nmod_mpoly_factor_init(h, ctx);
    nmod_mpolyv_init(v, ctx);

    success = _nmod_mpoly_factor_separable(h, f, ctx, 1);
    if (!success)
        goto cleanup;

    for (j = 0; j < h->num; j++)
    {
        success = h->num > 1 ? _factor_irred(v, h->poly + j, ctx, algo) :
                           _factor_irred_compressed(v, h->poly + j, ctx, algo);
        if (!success)
            goto cleanup;

        nmod_mpoly_factor_fit_length(g, g->num + v->length, ctx);
        for (k = 0; k < v->length; k++)
        {
            fmpz_mul(g->exp + g->num, h->exp + j, e);
            nmod_mpoly_swap(g->poly + g->num, v->coeffs + k, ctx);
            g->num++;
        }
    }

cleanup:

    nmod_mpoly_factor_clear(h, ctx);
    nmod_mpolyv_clear(v, ctx);

    return success;
}


int nmod_mpoly_factor_algo(
    nmod_mpoly_factor_t f,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong i, j;
    flint_bitcnt_t bits;
    nmod_mpoly_factor_t g;
    mpoly_compression_t M;

    if (!nmod_mpoly_factor_content(f, A, ctx))
        return 0;

    nmod_mpoly_factor_init(g, ctx);
    mpoly_compression_init(M);

    /* write into g */
    g->constant = f->constant;
    g->num = 0;
    for (i = 0; i < f->num; i++)
    {
        if (f->poly[i].length < 2)
        {
            nmod_mpoly_factor_fit_length(g, g->num + 1, ctx);
            nmod_mpoly_swap(g->poly + g->num, f->poly + i, ctx);
            fmpz_swap(g->exp + g->num, f->exp + i);
            g->num++;
            continue;
        }

        if (f->poly[i].bits > FLINT_BITS &&
            !nmod_mpoly_repack_bits_inplace(f->poly + i, FLINT_BITS, ctx))
        {
            success = 0;
            goto cleanup;
        }

        bits = f->poly[i].bits;

        mpoly_compression_set(M, f->poly[i].exps, bits, f->poly[i].length, ctx->minfo);
        if (M->is_irred)
        {
            nmod_mpoly_factor_fit_length(g, g->num + 1, ctx);
            nmod_mpoly_swap(g->poly + g->num, f->poly + i, ctx);
            fmpz_swap(g->exp + g->num, f->exp + i);
            g->num++;
        }
        else if (M->is_trivial)
        {
            success = _compressed_content_to_irred(g, f->poly + i, f->exp + i, ctx, algo);
            if (!success)
                goto cleanup;
        }
        else
        {
            nmod_mpoly_ctx_t Lctx;
            nmod_mpoly_t L;
            nmod_mpoly_factor_t h;

            /* compression may have messed up the content factorization */
            nmod_mpoly_ctx_init(Lctx, M->mvars, ORD_LEX, ctx->mod.n);
            nmod_mpoly_init(L, Lctx);
            nmod_mpoly_factor_init(h, Lctx);

            nmod_mpoly_compression_do(L, Lctx, f->poly[i].coeffs,
                                               f->poly[i].length, M);
            if (M->is_perm)
            {
                success = _compressed_content_to_irred(h, L, f->exp + i, Lctx, algo);
                fmpz_one(f->exp + i);
            }
            else
            {
                success = nmod_mpoly_factor_separable(h, L, Lctx, 1) &&
                          nmod_mpoly_factor_irred(h, Lctx, algo);
            }

            if (success)
            {
                FLINT_ASSERT(h->constant == 1);
                nmod_mpoly_factor_fit_length(g, g->num + h->num, ctx);
                for (j = 0; j < h->num; j++)
                {
                    fmpz_mul(g->exp + g->num, f->exp + i, h->exp + j);
                    nmod_mpoly_compression_undo(g->poly + g->num, bits, ctx,
                                                         h->poly + j, Lctx, M);
                    g->num++;
                }
            }

            nmod_mpoly_factor_clear(h, Lctx);
            nmod_mpoly_clear(L, Lctx);
            nmod_mpoly_ctx_clear(Lctx);

            if (!success)
                goto cleanup;
        }
    }

    nmod_mpoly_factor_swap(f, g, ctx);

    success = 1;

cleanup:

    nmod_mpoly_factor_clear(g, ctx);
    mpoly_compression_clear(M);

    FLINT_ASSERT(!success || nmod_mpoly_factor_matches(A, f, ctx));
    return success;
}


int nmod_mpoly_factor(
    nmod_mpoly_factor_t f,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    return nmod_mpoly_factor_algo(f, A, ctx, MPOLY_FACTOR_USE_ALL);
}


int nmod_mpoly_factor_zassenhaus(
    nmod_mpoly_factor_t f,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    return nmod_mpoly_factor_algo(f, A, ctx, MPOLY_FACTOR_USE_ZAS);
}


int nmod_mpoly_factor_wang(
    nmod_mpoly_factor_t f,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    return nmod_mpoly_factor_algo(f, A, ctx, MPOLY_FACTOR_USE_WANG);
}


int nmod_mpoly_factor_zippel(
    nmod_mpoly_factor_t f,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    return nmod_mpoly_factor_algo(f, A, ctx, MPOLY_FACTOR_USE_ZIP);
}

