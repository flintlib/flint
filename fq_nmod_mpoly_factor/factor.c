/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"

static void _deflate(
    fq_nmod_mpoly_t A,
    const ulong * strides,
    const slong * perm,
    const fq_nmod_mpoly_ctx_t ctx)
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

    return;

do_it:

    TMP_START;

    texps = (ulong *) TMP_ALLOC(2*nvars*sizeof(ulong));
    sexps = texps + nvars;

    for (i = 0; i < A->length; i++)
    {
        mpoly_get_monomial_ui(texps, A->exps + N*i, bits, ctx->minfo);

        for (j = 0; j < nvars; j++)
        {
            FLINT_ASSERT(0 == texps[j] % strides[j]);
            texps[j] = texps[j]/strides[j];
        }

        for (j = 0; j < nvars; j++)
            sexps[j] = texps[perm[j]];

        mpoly_set_monomial_ui(A->exps + N*i, sexps, bits, ctx->minfo);
    }

    TMP_END;

    fq_nmod_mpoly_sort_terms(A, ctx);
    fq_nmod_mpoly_make_monic(A, A, ctx);

    return;
}


static void _inflate(
    fq_nmod_mpoly_t A,
    flint_bitcnt_t bits,
    const ulong * strides,
    const slong * perm,
    const fq_nmod_mpoly_ctx_t ctx)
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

    fq_nmod_mpoly_repack_bits_inplace(A, bits, ctx);

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

    fq_nmod_mpoly_sort_terms(A, ctx);
    fq_nmod_mpoly_make_monic(A, A, ctx);

    return;
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
    fq_nmod_mpolyv_t Af,
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong i, j;
    slong nvars = ctx->minfo->nvars;
    slong * perm;
    ulong * strides, * texps;
    flint_bitcnt_t Abits;
    flint_rand_t state;
#if FLINT_WANT_ASSERT
    fq_nmod_mpoly_t Aorg;

    fq_nmod_mpoly_init(Aorg, ctx);
    fq_nmod_mpoly_set(Aorg, A, ctx);
#endif

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(n_fq_is_one(A->coeffs + 0, ctx->fqctx));
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    if (A->length < 2)
    {
        FLINT_ASSERT(A->length == 1);
        FLINT_ASSERT(!fq_nmod_mpoly_is_fq_nmod(A, ctx));
        fq_nmod_mpolyv_fit_length(Af, 1, ctx);
        fq_nmod_mpoly_swap(Af->coeffs + 0, A, ctx);
        Af->length = 1;
        return 1;
    }

    if (A->bits > FLINT_BITS &&
        !fq_nmod_mpoly_repack_bits_inplace(A, FLINT_BITS, ctx))
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
        ulong ppowt, ppow = fq_nmod_ctx_mod(ctx->fqctx).n;
        slong N = mpoly_words_per_exp_sp(Abits, ctx->minfo);

        while (!n_mul_checked(&ppowt, ppow, fq_nmod_ctx_mod(ctx->fqctx).n))
            ppow = ppowt;

        for (j = 0; j < nvars; j++)
        {
            strides[j] = ppow;
            perm[j] = j;
        }

        for (i = 0; i < A->length; i++)
        {
	        mpoly_get_monomial_ui(texps, A->exps + N*i, Abits, ctx->minfo);
            for (j = 0; j < nvars; j++)
                strides[j] = n_gcd(strides[j], texps[j]);
        }
    }

    /* find permutation with gcd(A, derivative(A, gen(perm[i]))) = 1 */
    for (i = 0; i < nvars; i++)
    {
        if (strides[i] == 1)
        {
            SLONG_SWAP(perm[0], perm[i]);
            break;
        }
    }

    if (nvars < 2)
    {
        fq_nmod_t c;
        fq_nmod_poly_t Au;
        fq_nmod_poly_factor_t Auf;

        FLINT_ASSERT(nvars == 1);

        fq_nmod_init(c, ctx->fqctx);
        fq_nmod_poly_init(Au, ctx->fqctx);
        fq_nmod_poly_factor_init(Auf, ctx->fqctx);

        FLINT_ASSERT(fq_nmod_mpoly_is_fq_nmod_poly(A, perm[0], ctx));
        success = fq_nmod_mpoly_get_fq_nmod_poly(Au, A, perm[0], ctx);
        FLINT_ASSERT(success);
        fq_nmod_poly_factor(Auf, c, Au, ctx->fqctx);

        FLINT_ASSERT(fq_nmod_is_one(c, ctx->fqctx));

        fq_nmod_mpolyv_fit_length(Af, Auf->num, ctx);
        Af->length = Auf->num; 
        for (i = 0; i < Auf->num; i++)
        {
            FLINT_ASSERT(Auf->exp[i] == 1);
            _fq_nmod_mpoly_set_fq_nmod_poly(Af->coeffs + i, Abits,
                       Auf->poly[i].coeffs, Auf->poly[i].length, perm[0], ctx);
        }

        fq_nmod_clear(c, ctx->fqctx);
        fq_nmod_poly_clear(Au, ctx->fqctx);
        fq_nmod_poly_factor_clear(Auf, ctx->fqctx);

        success = 1;
    }
    else if (nvars == 2)
    {
        n_fq_poly_t c;
        n_fq_bpoly_t Ab;
        n_fq_tpoly_t Abf;

        n_poly_init(c);
        n_bpoly_init(Ab);
        n_tpoly_init(Abf);

        fq_nmod_mpoly_get_n_fq_bpoly(Ab, A, perm[0], perm[1], ctx);
        success = n_fq_bpoly_factor_smprime(c, Abf, Ab, 1, ctx->fqctx);
        if (!success)
        {
            fq_nmod_mpoly_get_n_fq_bpoly(Ab, A, perm[0], perm[1], ctx);
            n_fq_bpoly_factor_lgprime(c, Abf, Ab, ctx->fqctx, state);
        }
        FLINT_ASSERT(n_poly_degree(c) == 0);

        fq_nmod_mpolyv_fit_length(Af, Abf->length, ctx);
        Af->length = Abf->length;
        for (i = 0; i < Abf->length; i++)
        {
            fq_nmod_mpoly_set_n_fq_bpoly(Af->coeffs + i, Abits,
                                       Abf->coeffs + i, perm[0], perm[1], ctx);
            fq_nmod_mpoly_make_monic(Af->coeffs + i, Af->coeffs + i, ctx);
        }

        n_poly_clear(c);
        n_bpoly_clear(Ab);
        n_tpoly_clear(Abf);

        success = 1;
    }
    else
    {
		fq_nmod_mpoly_t lcA;
        fq_nmod_mpoly_factor_t lcAf;

        fq_nmod_mpoly_init(lcA, ctx);
        fq_nmod_mpoly_factor_init(lcAf, ctx);

        _deflate(A, strides, perm, ctx);

        #if FLINT_WANT_ASSERT
        {
            fq_nmod_mpoly_t g;
            fq_nmod_mpoly_init(g, ctx);
            fq_nmod_mpoly_derivative(g, A, 0, ctx);
            FLINT_ASSERT(fq_nmod_mpoly_gcd(g, g, A, ctx));
            FLINT_ASSERT(fq_nmod_mpoly_is_one(g, ctx));
            fq_nmod_mpoly_clear(g, ctx);
        }
        #endif

        success = 0;

        if (algo & (MPOLY_FACTOR_USE_WANG | MPOLY_FACTOR_USE_ZIP))
        {
            _fq_nmod_mpoly_get_lead0(lcA, A, ctx);
            if (fq_nmod_mpoly_factor(lcAf, lcA, ctx))
            {
                if (algo & MPOLY_FACTOR_USE_ZIP)
                {
                    success = fq_nmod_mpoly_factor_irred_smprime_zippel(
                                                 Af, A, lcAf, lcA, ctx, state);
                    if (success == 0)
                        success = fq_nmod_mpoly_factor_irred_lgprime_zippel(
                                                 Af, A, lcAf, lcA, ctx, state);
                }

                if (algo & MPOLY_FACTOR_USE_WANG)
                {
                    if (success == 0)
                        success = fq_nmod_mpoly_factor_irred_smprime_wang(
                                                 Af, A, lcAf, lcA, ctx, state);
                    if (success == 0)
                        success = fq_nmod_mpoly_factor_irred_lgprime_wang(
                                                 Af, A, lcAf, lcA, ctx, state);
                }
            }
        }

        if (algo & MPOLY_FACTOR_USE_ZAS)
        {
            if (success == 0)
    		    success = fq_nmod_mpoly_factor_irred_smprime_zassenhaus(
                                                            Af, A, ctx, state);
            if (success == 0)
        		success = fq_nmod_mpoly_factor_irred_lgprime_zassenhaus(
                                                            Af, A, ctx, state);
        }

        success = (success > 0);
		if (success)
        {
		    for (i = 0; i < Af->length; i++)
                _inflate(Af->coeffs + i, Abits, strides, perm, ctx);
        }

		fq_nmod_mpoly_clear(lcA, ctx);
		fq_nmod_mpoly_factor_clear(lcAf, ctx);
    }

    flint_randclear(state);
    flint_free(strides);
    flint_free(perm);

#if FLINT_WANT_ASSERT
    if (success)
    {
        fq_nmod_mpoly_t prod;
        fq_nmod_mpoly_init(prod, ctx);
        fq_nmod_mpoly_one(prod, ctx);
        for (i = 0; i < Af->length; i++)
            fq_nmod_mpoly_mul(prod, prod, Af->coeffs + i, ctx);
        FLINT_ASSERT(fq_nmod_mpoly_equal(prod, Aorg, ctx));
        fq_nmod_mpoly_clear(prod, ctx);
        fq_nmod_mpoly_clear(Aorg, ctx);
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
    fq_nmod_mpolyv_t f,
    const fq_nmod_mpoly_ctx_t ctx,
    fq_nmod_mpolyv_t g)    /* temp */
{
    int success;
    slong v, i;
    fq_nmod_mpoly_struct * t;
    fq_nmod_mpoly_univar_t u;

    fq_nmod_mpoly_univar_init(u, ctx);

    /* first make primitive */
    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        g->length = 0;
        for (i = 0; i < f->length; i++)
        {
            fq_nmod_mpoly_to_univar(u, f->coeffs + i, v, ctx);
            FLINT_ASSERT(u->length > 0);
            FLINT_ASSERT(fmpz_is_zero(u->exps + u->length - 1));

            fq_nmod_mpolyv_fit_length(g, g->length + 2, ctx);
            success = _fq_nmod_mpoly_vec_content_mpoly(g->coeffs + g->length,
                                                    u->coeffs, u->length, ctx);
            if (!success)
                goto cleanup;

            if (fq_nmod_mpoly_is_fq_nmod(g->coeffs + g->length, ctx))
            {
                fq_nmod_mpoly_swap(g->coeffs + g->length, f->coeffs + i, ctx);
                g->length++;
            }
            else
            {
                success = fq_nmod_mpoly_divides(g->coeffs + g->length + 1,
                                    f->coeffs + i, g->coeffs + g->length, ctx);
                FLINT_ASSERT(success);

                if (fq_nmod_mpoly_is_fq_nmod(g->coeffs + g->length + 1, ctx))
                    g->length += 1;
                else
                    g->length += 2;
            }
        }

        fq_nmod_mpolyv_swap(f, g, ctx);
    }

    /* now make separable/derivative zero wrt each variable */
    fq_nmod_mpolyv_fit_length(g, 1, ctx);
    t = g->coeffs + 0;
    for (v = 0; v < ctx->minfo->nvars; v++)
    {
        i = 0;
        while (i < f->length)
        {
            fq_nmod_mpoly_derivative(t, f->coeffs + i, v, ctx);
            if (fq_nmod_mpoly_is_zero(t, ctx))
            {
                /* f[i] has zero derivative */
                i++;
                continue;
            }

            fq_nmod_mpolyv_fit_length(f, f->length + 1, ctx);

            success = fq_nmod_mpoly_gcd_cofactors(f->coeffs + f->length,
                                      f->coeffs + i, t, f->coeffs + i, t, ctx);
            if (!success)
                goto cleanup;

            if (fq_nmod_mpoly_is_fq_nmod(f->coeffs + f->length, ctx))
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

    fq_nmod_mpoly_univar_clear(u, ctx);

    return 1;
}

/*
    A is sep.

    return 1 for success, 0 for failure
*/
static int _factor_irred(
    fq_nmod_mpolyv_t Af,
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t Actx,
    unsigned int algo)
{
    int success;
    slong i, j;
    flint_bitcnt_t Abits;
    mpoly_compression_t M;
#if FLINT_WANT_ASSERT
    fq_nmod_mpoly_t Aorg;

    fq_nmod_mpoly_init(Aorg, Actx);
    fq_nmod_mpoly_set(Aorg, A, Actx);
#endif

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->coeffs[0] == 1);

    if (A->length < 2)
    {
        FLINT_ASSERT(A->length == 1);
        FLINT_ASSERT(!fq_nmod_mpoly_is_fq_nmod(A, Actx));
        fq_nmod_mpolyv_fit_length(Af, 1, Actx);
        Af->length = 1;
        fq_nmod_mpoly_swap(Af->coeffs + 0, A, Actx);
        success = 1;
        goto cleanup_less;
    }

    if (A->bits > FLINT_BITS &&
        !fq_nmod_mpoly_repack_bits_inplace(A, FLINT_BITS, Actx))
    {
        success = 0;
        goto cleanup_less;
    }

    Abits = A->bits;

    mpoly_compression_init(M);
    mpoly_compression_set(M, A->exps, A->bits, A->length, Actx->minfo);

    if (M->is_irred)
    {
        fq_nmod_mpolyv_fit_length(Af, 1, Actx);
        Af->length = 1;
        fq_nmod_mpoly_swap(Af->coeffs + 0, A, Actx);
        success = 1;
    }
    else if (M->is_trivial)
    {
        success = _factor_irred_compressed(Af, A, Actx, algo);
    }
    else
    {
        fq_nmod_mpoly_ctx_t Lctx;
        fq_nmod_mpolyv_t Lf, Lft, Lfs;

        fq_nmod_mpoly_ctx_init(Lctx, M->mvars, ORD_LEX, Actx->fqctx);
        fq_nmod_mpolyv_init(Lf, Lctx);
        fq_nmod_mpolyv_init(Lft, Lctx);
        fq_nmod_mpolyv_init(Lfs, Lctx);

        fq_nmod_mpolyv_fit_length(Lft, 1, Lctx);
        Lft->length = 1;
        fq_nmod_mpoly_compression_do(Lft->coeffs + 0, Lctx, A->coeffs, A->length, M);

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

                fq_nmod_mpolyv_fit_length(Lf, Lf->length + Lfs->length, Lctx);
                for (j = 0; j < Lfs->length; j++)
                {
                    fq_nmod_mpoly_swap(Lf->coeffs + Lf->length, Lfs->coeffs + j, Lctx);
                    Lf->length++;
                }
            }
        }

        if (success)
        {
            fq_nmod_mpolyv_fit_length(Af, Lf->length, Actx);
            Af->length = Lf->length;
            for (i = 0; i < Lf->length; i++)
            {
                fq_nmod_mpoly_compression_undo(Af->coeffs + i, Abits, Actx,
                                                      Lf->coeffs + i, Lctx, M);
            }
        }

        fq_nmod_mpolyv_clear(Lf, Lctx);
        fq_nmod_mpolyv_clear(Lft, Lctx);
        fq_nmod_mpolyv_clear(Lfs, Lctx);
        fq_nmod_mpoly_ctx_clear(Lctx);
    }

    mpoly_compression_clear(M);

cleanup_less:

#if FLINT_WANT_ASSERT
    if (success)
    {
        fq_nmod_mpoly_t prod;
        fq_nmod_mpoly_init(prod, Actx);
        fq_nmod_mpoly_one(prod, Actx);
        for (i = 0; i < Af->length; i++)
            fq_nmod_mpoly_mul(prod, prod, Af->coeffs + i, Actx);
        FLINT_ASSERT(fq_nmod_mpoly_equal(prod, Aorg, Actx));
        fq_nmod_mpoly_clear(prod, Actx);
        fq_nmod_mpoly_clear(Aorg, Actx);
    }
#endif

    FLINT_ASSERT(success == 0 || success == 1);
    return success;
}


/*
    Assum each factor in f is sep.
    Replace f by an irreducible factorization.
*/
int fq_nmod_mpoly_factor_irred(
    fq_nmod_mpoly_factor_t f,
    const fq_nmod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong i, j;
    fq_nmod_mpolyv_t t;
    fq_nmod_mpoly_factor_t g;

    fq_nmod_mpolyv_init(t, ctx);
    fq_nmod_mpoly_factor_init(g, ctx);

    fq_nmod_set(g->constant, f->constant, ctx->fqctx);
    g->num = 0;
    for (j = 0; j < f->num; j++)
    {
        success = _factor_irred(t, f->poly + j, ctx, algo);
        if (!success)
            goto cleanup;

        fq_nmod_mpoly_factor_fit_length(g, g->num + t->length, ctx);
        for (i = 0; i < t->length; i++)
        {
            fmpz_set(g->exp + g->num, f->exp + j);
            fq_nmod_mpoly_swap(g->poly + g->num, t->coeffs + i, ctx);
            g->num++;
        }
    }
    fq_nmod_mpoly_factor_swap(f, g, ctx);

    success = 1;

cleanup:

    fq_nmod_mpolyv_clear(t, ctx);
    fq_nmod_mpoly_factor_clear(g, ctx);

    return success;
}


/*
    append factor(f)^e to g
    assuming f is compressed and content free
*/
static int _compressed_content_to_irred(
    fq_nmod_mpoly_factor_t g,
    fq_nmod_mpoly_t f,
    const fmpz_t e,
    const fq_nmod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong j, k;
    fq_nmod_mpoly_factor_t h;
    fq_nmod_mpolyv_t v;

    fq_nmod_mpoly_factor_init(h, ctx);
    fq_nmod_mpolyv_init(v, ctx);

    success = _fq_nmod_mpoly_factor_separable(h, f, ctx, 1);
    if (!success)
        goto cleanup;

    for (j = 0; j < h->num; j++)
    {
        success = h->num > 1 ? _factor_irred(v, h->poly + j, ctx, algo) :
                           _factor_irred_compressed(v, h->poly + j, ctx, algo);
        if (!success)
            goto cleanup;

        fq_nmod_mpoly_factor_fit_length(g, g->num + v->length, ctx);
        for (k = 0; k < v->length; k++)
        {
            fmpz_mul(g->exp + g->num, h->exp + j, e);
            fq_nmod_mpoly_swap(g->poly + g->num, v->coeffs + k, ctx);
            g->num++;
        }
    }

cleanup:

    fq_nmod_mpoly_factor_clear(h, ctx);
    fq_nmod_mpolyv_clear(v, ctx);

    return success;
}

int fq_nmod_mpoly_factor_algo(
    fq_nmod_mpoly_factor_t f,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong i, j;
    flint_bitcnt_t bits;
    fq_nmod_mpoly_factor_t g;
    mpoly_compression_t M;

    if (!fq_nmod_mpoly_factor_content(f, A, ctx))
        return 0;

    fq_nmod_mpoly_factor_init(g, ctx);
    mpoly_compression_init(M);

    /* write into g */
    fq_nmod_set(g->constant, f->constant, ctx->fqctx);
    g->num = 0;
    for (i = 0; i < f->num; i++)
    {
        if (f->poly[i].length < 2)
        {
            fq_nmod_mpoly_factor_fit_length(g, g->num + 1, ctx);
            fq_nmod_mpoly_swap(g->poly + g->num, f->poly + i, ctx);
            fmpz_swap(g->exp + g->num, f->exp + i);
            g->num++;
            continue;
        }

        if (f->poly[i].bits > FLINT_BITS &&
            !fq_nmod_mpoly_repack_bits_inplace(f->poly + i, FLINT_BITS, ctx))
        {
            success = 0;
            goto cleanup;
        }

        bits = f->poly[i].bits;

        mpoly_compression_set(M, f->poly[i].exps, bits, f->poly[i].length, ctx->minfo);
        if (M->is_irred)
        {
            fq_nmod_mpoly_factor_fit_length(g, g->num + 1, ctx);
            fq_nmod_mpoly_swap(g->poly + g->num, f->poly + i, ctx);
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
            fq_nmod_mpoly_ctx_t Lctx;
            fq_nmod_mpoly_t L;
            fq_nmod_mpoly_factor_t h;

            /* compression may have messed up the content factorization */
            fq_nmod_mpoly_ctx_init(Lctx, M->mvars, ORD_LEX, ctx->fqctx);
            fq_nmod_mpoly_init(L, Lctx);
            fq_nmod_mpoly_factor_init(h, Lctx);

            fq_nmod_mpoly_compression_do(L, Lctx, f->poly[i].coeffs,
                                                  f->poly[i].length, M);
            if (M->is_perm)
            {
                success = _compressed_content_to_irred(h, L, f->exp + i, Lctx, algo);
                fmpz_one(f->exp + i);
            }
            else
            {
                success = fq_nmod_mpoly_factor_separable(h, L, Lctx, 1) &&
                          fq_nmod_mpoly_factor_irred(h, Lctx, algo);
            }

            if (success)
            {
                FLINT_ASSERT(fq_nmod_is_one(h->constant, ctx->fqctx));
                fq_nmod_mpoly_factor_fit_length(g, g->num + h->num, ctx);
                for (j = 0; j < h->num; j++)
                {
                    fmpz_mul(g->exp + g->num, f->exp + i, h->exp + j);
                    fq_nmod_mpoly_compression_undo(g->poly + g->num, bits, ctx,
                                                         h->poly + j, Lctx, M);
                    g->num++;
                }
            }

            fq_nmod_mpoly_factor_clear(h, Lctx);
            fq_nmod_mpoly_clear(L, Lctx);
            fq_nmod_mpoly_ctx_clear(Lctx);

            if (!success)
                goto cleanup;
        }
    }

    fq_nmod_mpoly_factor_swap(f, g, ctx);

    success = 1;

cleanup:

    fq_nmod_mpoly_factor_clear(g, ctx);
    mpoly_compression_clear(M);

    FLINT_ASSERT(!success || fq_nmod_mpoly_factor_matches(A, f, ctx));
    return success;
}


int fq_nmod_mpoly_factor(
    fq_nmod_mpoly_factor_t f,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    return fq_nmod_mpoly_factor_algo(f, A, ctx, MPOLY_FACTOR_USE_ALL);
}


int fq_nmod_mpoly_factor_zassenhaus(
    fq_nmod_mpoly_factor_t f,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    return fq_nmod_mpoly_factor_algo(f, A, ctx, MPOLY_FACTOR_USE_ZAS);
}


int fq_nmod_mpoly_factor_wang(
    fq_nmod_mpoly_factor_t f,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    return fq_nmod_mpoly_factor_algo(f, A, ctx, MPOLY_FACTOR_USE_WANG);
}


int fq_nmod_mpoly_factor_zippel(
    fq_nmod_mpoly_factor_t f,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    return fq_nmod_mpoly_factor_algo(f, A, ctx, MPOLY_FACTOR_USE_ZIP);
}

