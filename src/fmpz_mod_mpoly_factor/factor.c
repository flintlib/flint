/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_factor.h"
#include "long_extras.h"

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
    fmpz_mod_mpolyv_t Af,
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    slong i;
    int success;
    slong nvars = ctx->minfo->nvars;
    flint_bitcnt_t Abits;
    flint_rand_t state;
#if FLINT_WANT_ASSERT
    fmpz_mod_mpoly_t Aorg;
#endif

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(fmpz_is_one(A->coeffs + 0));
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(!fmpz_abs_fits_ui(fmpz_mod_ctx_modulus(ctx->ffinfo)));

    if (A->length < 2)
    {
        FLINT_ASSERT(A->length == 1);
        FLINT_ASSERT(!fmpz_mod_mpoly_is_fmpz(A, ctx));
        fmpz_mod_mpolyv_fit_length(Af, 1, ctx);
        fmpz_mod_mpoly_swap(Af->coeffs + 0, A, ctx);
        Af->length = 1;
        return 1;
    }

    if (A->bits > FLINT_BITS &&
        !fmpz_mod_mpoly_repack_bits_inplace(A, FLINT_BITS, ctx))
    {
        return 0;
    }

#if FLINT_WANT_ASSERT
    fmpz_mod_mpoly_init(Aorg, ctx);
    fmpz_mod_mpoly_set(Aorg, A, ctx);
#endif

    Abits = A->bits;

    flint_randinit(state);

    if (nvars < 2)
    {
        fmpz_mod_poly_t Au;
        fmpz_mod_poly_factor_t Auf;

        FLINT_ASSERT(nvars == 1);

        fmpz_mod_poly_init(Au, ctx->ffinfo);
        fmpz_mod_poly_factor_init(Auf, ctx->ffinfo);

        FLINT_ASSERT(fmpz_mod_mpoly_is_fmpz_mod_poly(A, 0, ctx));
        success = fmpz_mod_mpoly_get_fmpz_mod_poly(Au, A, 0, ctx);
        FLINT_ASSERT(success);
        fmpz_mod_poly_factor(Auf, Au, ctx->ffinfo);

        fmpz_mod_mpolyv_fit_length(Af, Auf->num, ctx);
        Af->length = Auf->num; 
        for (i = 0; i < Auf->num; i++)
        {
            FLINT_ASSERT(Auf->exp[i] == 1);
            _fmpz_mod_mpoly_set_fmpz_mod_poly(Af->coeffs + i, Abits,
                             Auf->poly[i].coeffs, Auf->poly[i].length, 0, ctx);
        }

        fmpz_mod_poly_clear(Au, ctx->ffinfo);
        fmpz_mod_poly_factor_clear(Auf, ctx->ffinfo);

        success = 1;
    }
    else if (nvars == 2)
    {
        fmpz_mod_poly_t c;
        fmpz_mod_bpoly_t Ab;
        fmpz_mod_tpoly_t Abf;

        fmpz_mod_poly_init(c, ctx->ffinfo);
        fmpz_mod_bpoly_init(Ab, ctx->ffinfo);
        fmpz_mod_tpoly_init(Abf, ctx->ffinfo);

        fmpz_mod_mpoly_get_fmpz_mod_bpoly(Ab, A, 0, 1, ctx);
        success = fmpz_mod_bpoly_factor_smprime(c, Abf, Ab, 1, ctx->ffinfo);

        FLINT_ASSERT(!success || fmpz_mod_poly_degree(c, ctx->ffinfo) == 0);

        fmpz_mod_mpolyv_fit_length(Af, Abf->length, ctx);
        Af->length = Abf->length;
        for (i = 0; i < Abf->length; i++)
        {
            fmpz_mod_mpoly_set_fmpz_mod_bpoly(Af->coeffs + i, Abits,
                                                   Abf->coeffs + i, 0, 1, ctx);
            fmpz_mod_mpoly_make_monic(Af->coeffs + i, Af->coeffs + i, ctx);
        }

        fmpz_mod_poly_clear(c, ctx->ffinfo);
        fmpz_mod_bpoly_clear(Ab, ctx->ffinfo);
        fmpz_mod_tpoly_clear(Abf, ctx->ffinfo);
    }
    else
    {
		fmpz_mod_mpoly_t lcA;
        fmpz_mod_mpoly_factor_t lcAf;

        fmpz_mod_mpoly_init(lcA, ctx);
        fmpz_mod_mpoly_factor_init(lcAf, ctx);

        #if FLINT_WANT_ASSERT
        {
            fmpz_mod_mpoly_t g;
            fmpz_mod_mpoly_init(g, ctx);
            fmpz_mod_mpoly_derivative(g, A, 0, ctx);
            FLINT_ASSERT(fmpz_mod_mpoly_gcd(g, g, A, ctx));
            FLINT_ASSERT(fmpz_mod_mpoly_is_one(g, ctx));
            fmpz_mod_mpoly_clear(g, ctx);
        }
        #endif

        success = 0;

        if (!(algo & (MPOLY_FACTOR_USE_WANG | MPOLY_FACTOR_USE_ZIP)))
            goto try_zassenhaus;

        /* TODO lcc_kaltofen */
        _fmpz_mod_mpoly_get_lead0(lcA, A, ctx);
        if (!fmpz_mod_mpoly_factor(lcAf, lcA, ctx))
            goto try_zassenhaus;

        if (!(algo & MPOLY_FACTOR_USE_ZIP))
        {
            if (success == 0)
                success = fmpz_mod_mpoly_factor_irred_smprime_wang(
                                                Af, A, lcAf, lcA, ctx, state);
        }
        else if (!(algo & MPOLY_FACTOR_USE_WANG))
        {
            if (success == 0)
                success = fmpz_mod_mpoly_factor_irred_smprime_zippel(
                                                 Af, A, lcAf, lcA, ctx, state);
        }
        else
        {
            double tdensity = 0;
            fmpz_t x;
            fmpz_init(x);
            fmpz_mod_mpoly_total_degree_fmpz(x, A, ctx);
            if (fmpz_fits_si(x))
            {
                fmpz_bin_uiui(x, fmpz_get_si(x) + nvars, nvars);
                tdensity = A->length/fmpz_get_d(x);
            }
            fmpz_clear(x);

            if (tdensity > 0.005)
            {
                if (success == 0)
                    success = fmpz_mod_mpoly_factor_irred_smprime_wang(
                                                 Af, A, lcAf, lcA, ctx, state);
                if (success == 0)
                    success = fmpz_mod_mpoly_factor_irred_smprime_zippel(
                                                 Af, A, lcAf, lcA, ctx, state);
            }
            else
            {
                if (success == 0)
                    success = fmpz_mod_mpoly_factor_irred_smprime_zippel(
                                                 Af, A, lcAf, lcA, ctx, state);
                if (success == 0)
                    success = fmpz_mod_mpoly_factor_irred_smprime_wang(
                                                 Af, A, lcAf, lcA, ctx, state);
            }
        }

    try_zassenhaus:

        if (algo & MPOLY_FACTOR_USE_ZAS)
        {
            if (success == 0)
    		    success = fmpz_mod_mpoly_factor_irred_smprime_zassenhaus(
                                                            Af, A, ctx, state);
        }

        success = (success > 0);

	    fmpz_mod_mpoly_clear(lcA, ctx);
        fmpz_mod_mpoly_factor_clear(lcAf, ctx);
    }

    flint_randclear(state);

#if FLINT_WANT_ASSERT
    if (success)
    {
        fmpz_mod_mpoly_t prod;
        fmpz_mod_mpoly_init(prod, ctx);
        fmpz_mod_mpoly_one(prod, ctx);
        for (i = 0; i < Af->length; i++)
            fmpz_mod_mpoly_mul(prod, prod, Af->coeffs + i, ctx);

        FLINT_ASSERT(fmpz_mod_mpoly_equal(prod, Aorg, ctx));
        fmpz_mod_mpoly_clear(prod, ctx);
    }
    fmpz_mod_mpoly_clear(Aorg, ctx);
#endif

    FLINT_ASSERT(success == 0 || success == 1);
    return success;
}


/*
    f is already squarefree
    make the factors in f have the sep property
*/
static int _refine_sep(
    fmpz_mod_mpolyv_t f,
    const fmpz_mod_mpoly_ctx_t ctx,
    fmpz_mod_mpolyv_t g)    /* temp */
{
    int success;
    slong v, i;
    fmpz_mod_mpoly_struct * t;
    fmpz_mod_mpoly_univar_t u;

    fmpz_mod_mpoly_univar_init(u, ctx);

    /* first make primitive */
    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        g->length = 0;
        for (i = 0; i < f->length; i++)
        {
            fmpz_mod_mpoly_to_univar(u, f->coeffs + i, v, ctx);
            FLINT_ASSERT(u->length > 0);
            FLINT_ASSERT(fmpz_is_zero(u->exps + u->length - 1));

            fmpz_mod_mpolyv_fit_length(g, g->length + 2, ctx);
            success = _fmpz_mod_mpoly_vec_content_mpoly(g->coeffs + g->length,
                                                    u->coeffs, u->length, ctx);
            if (!success)
                goto cleanup;

            if (fmpz_mod_mpoly_is_fmpz(g->coeffs + g->length, ctx))
            {
                fmpz_mod_mpoly_swap(g->coeffs + g->length, f->coeffs + i, ctx);
                g->length++;
            }
            else
            {
                success = fmpz_mod_mpoly_divides(g->coeffs + g->length + 1,
                                    f->coeffs + i, g->coeffs + g->length, ctx);
                FLINT_ASSERT(success);

                if (fmpz_mod_mpoly_is_fmpz(g->coeffs + g->length + 1, ctx))
                    g->length += 1;
                else
                    g->length += 2;
            }
        }

        fmpz_mod_mpolyv_swap(f, g, ctx);
    }

    /* now make separable/derivative zero wrt each variable */
    fmpz_mod_mpolyv_fit_length(g, 1, ctx);
    t = g->coeffs + 0;
    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        i = 0;
        while (i < f->length)
        {
            fmpz_mod_mpoly_derivative(t, f->coeffs + i, v, ctx);
            if (fmpz_mod_mpoly_is_zero(t, ctx))
            {
                /* f[i] has zero derivative */
                FLINT_ASSERT(fmpz_mod_mpoly_degree_si(f->coeffs + i, v, ctx) == 0);
                i++;
                continue;
            }

            fmpz_mod_mpolyv_fit_length(f, f->length + 1, ctx);

            success = fmpz_mod_mpoly_gcd_cofactors(f->coeffs + f->length,
                                      f->coeffs + i, t, f->coeffs + i, t, ctx);
            if (!success)
                goto cleanup;

            if (fmpz_mod_mpoly_is_fmpz(f->coeffs + f->length, ctx))
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

    fmpz_mod_mpoly_univar_clear(u, ctx);

    return 1;
}


/*
    A is sep.

    return 1 for success, 0 for failure
*/
static int _factor_irred(
    fmpz_mod_mpolyv_t Af,
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t Actx,
    unsigned int algo)
{
    int success;
    slong i, j;
    flint_bitcnt_t Abits;
    mpoly_compression_t M;
#if FLINT_WANT_ASSERT
    fmpz_mod_mpoly_t Aorg;

    fmpz_mod_mpoly_init(Aorg, Actx);
    fmpz_mod_mpoly_set(Aorg, A, Actx);
#endif

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->coeffs[0] == 1);
    FLINT_ASSERT(!fmpz_abs_fits_ui(fmpz_mod_ctx_modulus(Actx->ffinfo)));

    if (A->length < 2)
    {
        FLINT_ASSERT(A->length == 1);
        FLINT_ASSERT(!fmpz_mod_mpoly_is_fmpz(A, Actx));
        fmpz_mod_mpolyv_fit_length(Af, 1, Actx);
        Af->length = 1;
        fmpz_mod_mpoly_swap(Af->coeffs + 0, A, Actx);
        success = 1;
        goto cleanup_less;
    }

    if (A->bits > FLINT_BITS &&
        !fmpz_mod_mpoly_repack_bits_inplace(A, FLINT_BITS, Actx))
    {
        success = 0;
        goto cleanup_less;
    }

    Abits = A->bits;

    mpoly_compression_init(M);
    mpoly_compression_set(M, A->exps, A->bits, A->length, Actx->minfo);

    if (M->is_irred)
    {
        fmpz_mod_mpolyv_fit_length(Af, 1, Actx);
        Af->length = 1;
        fmpz_mod_mpoly_swap(Af->coeffs + 0, A, Actx);
        success = 1;
    }
    else if (M->is_trivial)
    {
        success = _factor_irred_compressed(Af, A, Actx, algo);
    }
    else
    {
        fmpz_mod_mpoly_ctx_t Lctx;
        fmpz_mod_mpolyv_t Lf, Lft, Lfs;

        fmpz_mod_mpoly_ctx_init(Lctx, M->mvars, ORD_LEX, fmpz_mod_ctx_modulus(Actx->ffinfo));
        fmpz_mod_mpolyv_init(Lf, Lctx);
        fmpz_mod_mpolyv_init(Lft, Lctx);
        fmpz_mod_mpolyv_init(Lfs, Lctx);

        fmpz_mod_mpolyv_fit_length(Lft, 1, Lctx);
        Lft->length = 1;
        fmpz_mod_mpoly_compression_do(Lft->coeffs + 0, Lctx, A->coeffs, A->length, M);

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

                fmpz_mod_mpolyv_fit_length(Lf, Lf->length + Lfs->length, Lctx);
                for (j = 0; j < Lfs->length; j++)
                {
                    fmpz_mod_mpoly_swap(Lf->coeffs + Lf->length, Lfs->coeffs + j, Lctx);
                    Lf->length++;
                }
            }
        }

        if (success)
        {
            fmpz_mod_mpolyv_fit_length(Af, Lf->length, Actx);
            Af->length = Lf->length;
            for (i = 0; i < Lf->length; i++)
            {
                fmpz_mod_mpoly_compression_undo(Af->coeffs + i, Abits, Actx,
                                                      Lf->coeffs + i, Lctx, M);
            }
        }

        fmpz_mod_mpolyv_clear(Lf, Lctx);
        fmpz_mod_mpolyv_clear(Lft, Lctx);
        fmpz_mod_mpolyv_clear(Lfs, Lctx);
        fmpz_mod_mpoly_ctx_clear(Lctx);
    }

    mpoly_compression_clear(M);

cleanup_less:

#if FLINT_WANT_ASSERT
    if (success)
    {
        fmpz_mod_mpoly_t prod;
        fmpz_mod_mpoly_init(prod, Actx);
        fmpz_mod_mpoly_one(prod, Actx);
        for (i = 0; i < Af->length; i++)
            fmpz_mod_mpoly_mul(prod, prod, Af->coeffs + i, Actx);
        FLINT_ASSERT(fmpz_mod_mpoly_equal(prod, Aorg, Actx));
        fmpz_mod_mpoly_clear(prod, Actx);
        fmpz_mod_mpoly_clear(Aorg, Actx);
    }
#endif

    FLINT_ASSERT(success == 0 || success == 1);
    return success;
}


/*
    Assume each factor in f is sep.
    Replace f by an irreducible factorization.
*/
int fmpz_mod_mpoly_factor_irred(
    fmpz_mod_mpoly_factor_t f,
    const fmpz_mod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong i, j;
    fmpz_mod_mpolyv_t t;
    fmpz_mod_mpoly_factor_t g;

    fmpz_mod_mpolyv_init(t, ctx);
    fmpz_mod_mpoly_factor_init(g, ctx);

    fmpz_swap(g->constant, f->constant);
    g->num = 0;
    for (j = 0; j < f->num; j++)
    {
        success = _factor_irred(t, f->poly + j, ctx, algo);
        if (!success)
            goto cleanup;

        fmpz_mod_mpoly_factor_fit_length(g, g->num + t->length, ctx);
        for (i = 0; i < t->length; i++)
        {
            fmpz_set(g->exp + g->num, f->exp + j);
            fmpz_mod_mpoly_swap(g->poly + g->num, t->coeffs + i, ctx);
            g->num++;
        }
    }
    fmpz_mod_mpoly_factor_swap(f, g, ctx);

    success = 1;

cleanup:

    fmpz_mod_mpolyv_clear(t, ctx);
    fmpz_mod_mpoly_factor_clear(g, ctx);

    return success;
}


/*
    append factor(f)^e to g
    assuming f is compressed and content free
*/
static int _compressed_content_to_irred(
    fmpz_mod_mpoly_factor_t g,
    fmpz_mod_mpoly_t f,
    const fmpz_t e,
    const fmpz_mod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong j, k;
    fmpz_mod_mpoly_factor_t h;
    fmpz_mod_mpolyv_t v;

    FLINT_ASSERT(!fmpz_abs_fits_ui(fmpz_mod_ctx_modulus(ctx->ffinfo)));

    fmpz_mod_mpoly_factor_init(h, ctx);
    fmpz_mod_mpolyv_init(v, ctx);

    success = _fmpz_mod_mpoly_factor_separable(h, f, ctx, 1);
    if (!success)
        goto cleanup;

    for (j = 0; j < h->num; j++)
    {
        success = h->num > 1 ? _factor_irred(v, h->poly + j, ctx, algo) :
                           _factor_irred_compressed(v, h->poly + j, ctx, algo);
        if (!success)
            goto cleanup;

        fmpz_mod_mpoly_factor_fit_length(g, g->num + v->length, ctx);
        for (k = 0; k < v->length; k++)
        {
            fmpz_mul(g->exp + g->num, h->exp + j, e);
            fmpz_mod_mpoly_swap(g->poly + g->num, v->coeffs + k, ctx);
            g->num++;
        }
    }

cleanup:

    fmpz_mod_mpoly_factor_clear(h, ctx);
    fmpz_mod_mpolyv_clear(v, ctx);

    return success;
}


int fmpz_mod_mpoly_factor_algo(
    fmpz_mod_mpoly_factor_t f,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong i, j;
    flint_bitcnt_t bits;
    fmpz_mod_mpoly_factor_t g;
    mpoly_compression_t M;

    if (fmpz_abs_fits_ui(fmpz_mod_ctx_modulus(ctx->ffinfo)))
    {
        nmod_mpoly_ctx_t nctx;
        nmod_mpoly_t nA;
        nmod_mpoly_factor_t nf;

        *nctx->minfo = *ctx->minfo;
        nmod_init(&nctx->mod, fmpz_get_ui(fmpz_mod_ctx_modulus(ctx->ffinfo)));
        nmod_mpoly_init(nA, nctx);
        nmod_mpoly_factor_init(nf, nctx);

        _fmpz_mod_mpoly_get_nmod_mpoly(nA, nctx, A, ctx);
        success = nmod_mpoly_factor_algo(nf, nA, nctx, algo);
        _fmpz_mod_mpoly_factor_set_nmod_mpoly_factor(f, ctx, nf, nctx);

        nmod_mpoly_factor_clear(nf, nctx);
        nmod_mpoly_clear(nA, nctx);

        return success;
    }

    if (!fmpz_mod_mpoly_factor_content(f, A, ctx))
        return 0;

    fmpz_mod_mpoly_factor_init(g, ctx);
    mpoly_compression_init(M);

    /* write into g */
    fmpz_swap(g->constant, f->constant);
    g->num = 0;
    for (i = 0; i < f->num; i++)
    {
        if (f->poly[i].length < 2)
        {
            fmpz_mod_mpoly_factor_fit_length(g, g->num + 1, ctx);
            fmpz_mod_mpoly_swap(g->poly + g->num, f->poly + i, ctx);
            fmpz_swap(g->exp + g->num, f->exp + i);
            g->num++;
            continue;
        }

        if (f->poly[i].bits > FLINT_BITS &&
            !fmpz_mod_mpoly_repack_bits_inplace(f->poly + i, FLINT_BITS, ctx))
        {
            success = 0;
            goto cleanup;
        }

        bits = f->poly[i].bits;

        mpoly_compression_set(M, f->poly[i].exps, bits, f->poly[i].length, ctx->minfo);
        if (M->is_irred)
        {
            fmpz_mod_mpoly_factor_fit_length(g, g->num + 1, ctx);
            fmpz_mod_mpoly_swap(g->poly + g->num, f->poly + i, ctx);
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
            fmpz_mod_mpoly_ctx_t Lctx;
            fmpz_mod_mpoly_t L;
            fmpz_mod_mpoly_factor_t h;

            /* compression may have messed up the content factorization */
            fmpz_mod_mpoly_ctx_init(Lctx, M->mvars, ORD_LEX, fmpz_mod_ctx_modulus(ctx->ffinfo));
            fmpz_mod_mpoly_init(L, Lctx);
            fmpz_mod_mpoly_factor_init(h, Lctx);

            fmpz_mod_mpoly_compression_do(L, Lctx, f->poly[i].coeffs,
                                               f->poly[i].length, M);
            if (M->is_perm)
            {
                success = _compressed_content_to_irred(h, L, f->exp + i, Lctx, algo);
                fmpz_one(f->exp + i);
            }
            else
            {
                success = fmpz_mod_mpoly_factor_separable(h, L, Lctx, 1) &&
                          fmpz_mod_mpoly_factor_irred(h, Lctx, algo);
            }

            if (success)
            {
                FLINT_ASSERT(fmpz_is_one(h->constant));
                fmpz_mod_mpoly_factor_fit_length(g, g->num + h->num, ctx);
                for (j = 0; j < h->num; j++)
                {
                    fmpz_mul(g->exp + g->num, f->exp + i, h->exp + j);
                    fmpz_mod_mpoly_compression_undo(g->poly + g->num, bits, ctx,
                                                         h->poly + j, Lctx, M);
                    g->num++;
                }
            }

            fmpz_mod_mpoly_factor_clear(h, Lctx);
            fmpz_mod_mpoly_clear(L, Lctx);
            fmpz_mod_mpoly_ctx_clear(Lctx);

            if (!success)
                goto cleanup;
        }
    }

    fmpz_mod_mpoly_factor_swap(f, g, ctx);

    success = 1;

cleanup:

    fmpz_mod_mpoly_factor_clear(g, ctx);
    mpoly_compression_clear(M);

    FLINT_ASSERT(!success || fmpz_mod_mpoly_factor_matches(A, f, ctx));
    return success;
}


int fmpz_mod_mpoly_factor(
    fmpz_mod_mpoly_factor_t f,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    return fmpz_mod_mpoly_factor_algo(f, A, ctx, MPOLY_FACTOR_USE_ALL);
}


int fmpz_mod_mpoly_factor_zassenhaus(
    fmpz_mod_mpoly_factor_t f,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    return fmpz_mod_mpoly_factor_algo(f, A, ctx, MPOLY_FACTOR_USE_ZAS);
}


int fmpz_mod_mpoly_factor_wang(
    fmpz_mod_mpoly_factor_t f,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    return fmpz_mod_mpoly_factor_algo(f, A, ctx, MPOLY_FACTOR_USE_WANG);
}


int fmpz_mod_mpoly_factor_zippel(
    fmpz_mod_mpoly_factor_t f,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    return fmpz_mod_mpoly_factor_algo(f, A, ctx, MPOLY_FACTOR_USE_ZIP);
}
