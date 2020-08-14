/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"
#include "fq_nmod_mpoly_factor.h"


void nmod_mpoly_convert_perm(
    nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const nmod_mpoly_ctx_t Actx,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t Bctx,
    const slong * perm)
{
    slong n = Bctx->minfo->nvars;
    slong m = Actx->minfo->nvars;
    slong i, k, l;
    slong NA, NB;
    ulong * Aexps;
    ulong * Bexps;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    TMP_START;

    Aexps = (ulong *) TMP_ALLOC(m*sizeof(ulong));
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, Actx->minfo);
    NB = mpoly_words_per_exp(B->bits, Bctx->minfo);

    nmod_mpoly_fit_length_set_bits(A, B->length, Abits, Actx);
    A->length = B->length;
    for (i = 0; i < B->length; i++)
    {        
	    A->coeffs[i] = B->coeffs[i];
	    mpoly_get_monomial_ui(Bexps, B->exps + NB*i, B->bits, Bctx->minfo);
	    for (k = 0; k < m; k++)
	    {
	        l = perm[k];
	        Aexps[k] = l < 0 ? 0 : Bexps[l];
	    }
	    mpoly_set_monomial_ui(A->exps + NA*i, Aexps, Abits, Actx->minfo);
     }  
    TMP_END;
    nmod_mpoly_sort_terms(A, Actx);
}

#define USE_ZAS 1
#define USE_WANG 2
#define USE_ZIP 4

/*
    A is primitive w.r.t to any variable appearing in A.
    A is squarefree and monic.

    return 1 for success, 0 for failure
*/
static int _irreducible_factors(
    nmod_mpolyv_t Af,
    nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong i, t, nzdvar, mvars;
    slong nvars = ctx->minfo->nvars;
    slong * Adegs, * perm, * iperm;
    nmod_mpoly_t G, Abar, Bbar, nzdpoly;
    flint_bitcnt_t Lbits, Abits;
    int perm_is_id;
    flint_rand_t state;
#if WANT_ASSERT
    nmod_mpoly_t Aorg;

    nmod_mpoly_init(Aorg, ctx);
    nmod_mpoly_set(Aorg, A, ctx);
#endif

    nmod_mpoly_init(G, ctx);
    nmod_mpoly_init(Abar, ctx);
    nmod_mpoly_init(Bbar, ctx);
    nmod_mpoly_init(nzdpoly, ctx);
    flint_randinit(state);
    Adegs = (slong *) flint_malloc(3*nvars*sizeof(slong));
    perm = Adegs + nvars;
    iperm = perm + nvars;

    if (A->length < 2)
    {
        FLINT_ASSERT(A->length == 1);
        FLINT_ASSERT(!nmod_mpoly_is_ui(A, ctx));
        nmod_mpolyv_fit_length(Af, 1, ctx);
        Af->length = 1;
        nmod_mpoly_swap(Af->coeffs + 0, A, ctx);
        success = 1;
        goto cleanup;
    }

    if (!nmod_mpoly_degrees_fit_si(A, ctx))
    {
        success = 0;
        goto cleanup;
    }

    if (A->bits > FLINT_BITS &&
        !nmod_mpoly_repack_bits_inplace(A, FLINT_BITS, ctx))
    {
        success = 0;
        goto cleanup;
    }

    nmod_mpoly_degrees_si(Adegs, A, ctx);

    Abits = A->bits;

    mvars = 0;
    Lbits = 0;
    nzdvar = -1;
    for (i = 0; i < nvars; i++)
    {
		iperm[i] = -1;
        if (Adegs[i] > 0)
        {
            flint_bitcnt_t this_bits = FLINT_BIT_COUNT(Adegs[i]);
            Lbits = FLINT_MAX(Lbits, this_bits);
            perm[mvars] = i;
            mvars++;

            if (nzdvar < 0)
            {
                nmod_mpoly_derivative(nzdpoly, A, i, ctx);
                if (!nmod_mpoly_is_zero(nzdpoly, ctx))
                {
                    nzdvar = i;
                    t = perm[mvars - 1];
                    perm[mvars - 1] = perm[0];
                    perm[0] = t;
                }
            }
        }
    }

    FLINT_ASSERT(nzdvar >= 0);
    FLINT_ASSERT(perm[0] == nzdvar);
    FLINT_ASSERT(!nmod_mpoly_is_zero(nzdpoly, ctx));

    /*
        Check for annoying things like (x^2+y)(x^p+y). The squarefree
        requirement should already have ruled out (x^2+y)(x^p+y^p).
    */
    if (!nmod_mpoly_gcd_cofactors(G, Abar, Bbar, A, nzdpoly, ctx))
    {
        success = 0;
        goto cleanup;
    }
    if (!nmod_mpoly_is_one(G, ctx))
    {
        nmod_mpolyv_t Gf;
        nmod_mpolyv_init(Gf, ctx);
        success = _irreducible_factors(Af, Abar, ctx, algo);
        success = success && _irreducible_factors(Gf, G, ctx, algo);
        if (success)
        {
            nmod_mpolyv_fit_length(Af, Af->length + Gf->length, ctx);
            for (i = 0; i < Gf->length; i++)
            {
                nmod_mpoly_swap(Af->coeffs + Af->length, Gf->coeffs + i, ctx);
                Af->length++;
            }
        }
        nmod_mpolyv_clear(Gf, ctx);
        goto cleanup;
    }

	/* A is separable wrt gen(perm[0]) now */

    if (Lbits > FLINT_BITS - 10)
    {
        success = 0;
        goto cleanup;
    }

    /* TODO nice permutation */

    /* invert perm */
    perm_is_id = (mvars == nvars);
    for (i = 0; i < mvars; i++)
    {
        perm_is_id = perm_is_id && (perm[i] == i);
        iperm[perm[i]] = i;
    }

    if (mvars < 2)
    {
        nmod_poly_t Au;
        nmod_poly_factor_t Auf;

        FLINT_ASSERT(mvars == 1);

        nmod_poly_init_mod(Au, ctx->ffinfo->mod);
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
    else if (mvars == 2)
    {
        n_poly_t c;
        n_bpoly_t Ab;
        n_tpoly_t Abf;

        n_poly_init(c);
        n_bpoly_init(Ab);
        n_tpoly_init(Abf);

        nmod_mpoly_get_bpoly(Ab, A, perm[0], perm[1], ctx);
        success = n_bpoly_mod_factor_smprime(c, Abf, Ab, 1, ctx->ffinfo->mod);
        if (!success)
        {
            nmod_mpoly_get_bpoly(Ab, A, perm[0], perm[1], ctx);
            n_bpoly_mod_factor_lgprime(c, Abf, Ab, ctx->ffinfo->mod);
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
		nmod_mpoly_ctx_t Lctx;
		nmod_mpoly_t L, lcL;
		nmod_mpolyv_t Lf;
        nmod_mpoly_factor_t lcLf;

		nmod_mpoly_ctx_init(Lctx, mvars, ORD_LEX, ctx->ffinfo->mod.n);
		nmod_mpoly_init(L, Lctx);
        nmod_mpoly_init(lcL, Lctx);
		nmod_mpolyv_init(Lf, Lctx);
        nmod_mpoly_factor_init(lcLf, Lctx);

        Lbits = mpoly_fix_bits(Lbits + 1, Lctx->minfo);

		nmod_mpoly_convert_perm(L, Lbits, Lctx, A, ctx, perm);
        nmod_mpoly_make_monic(L, L, ctx);

        success = 0;

        if (algo & (USE_WANG | USE_ZIP))
        {
            _nmod_mpoly_get_lead0(lcL, L, Lctx);
            if (nmod_mpoly_factor(lcLf, lcL, Lctx))
            {
                if (algo & USE_ZIP)
                {
                    success = nmod_mpoly_factor_irred_smprime_zippel(
                                                Lf, L, lcLf, lcL, Lctx, state);
                    if (success == 0)
                        success = nmod_mpoly_factor_irred_medprime_zippel(
                                                Lf, L, lcLf, lcL, Lctx, state);
                    if (success == 0)
                        success = nmod_mpoly_factor_irred_lgprime_zippel(
                                                Lf, L, lcLf, lcL, Lctx, state);
                }

                if (algo & USE_WANG)
                {
                    if (success == 0)
                        success = nmod_mpoly_factor_irred_smprime_wang(
                                                Lf, L, lcLf, lcL, Lctx, state);
                    if (success == 0)
                        success = nmod_mpoly_factor_irred_medprime_wang(
                                                Lf, L, lcLf, lcL, Lctx, state);
                    if (success == 0)
                        success = nmod_mpoly_factor_irred_lgprime_wang(
                                                Lf, L, lcLf, lcL, Lctx, state);
                }
            }
        }

        if (algo & USE_ZAS)
        {
            if (success == 0)
    		    success = nmod_mpoly_factor_irred_smprime_zassenhaus(
                                                           Lf, L, Lctx, state);
            if (success == 0)
        		success = nmod_mpoly_factor_irred_medprime_zassenhaus(
                                                           Lf, L, Lctx, state);
            if (success == 0)
        		success = nmod_mpoly_factor_irred_lgprime_zassenhaus(
                                                           Lf, L, Lctx, state);
        }

        success = (success > 0);
		if (success)
        {
            nmod_mpolyv_fit_length(Af, Lf->length, ctx);
            Af->length = Lf->length;
		    for (i = 0; i < Lf->length; i++)
		    {
                nmod_mpoly_convert_perm(Af->coeffs + i, Abits, ctx,
                                                  Lf->coeffs + i, Lctx, iperm);
                nmod_mpoly_make_monic(Af->coeffs + i, Af->coeffs + i, ctx);
		    }
        }

	    nmod_mpoly_clear(L, Lctx);
	    nmod_mpoly_clear(lcL, Lctx);
	    nmod_mpolyv_clear(Lf, Lctx);
        nmod_mpoly_factor_clear(lcLf, Lctx);
	    nmod_mpoly_ctx_clear(Lctx);
    }

cleanup:

    nmod_mpoly_clear(G, ctx);
    nmod_mpoly_clear(Abar, ctx);
    nmod_mpoly_clear(Bbar, ctx);
    nmod_mpoly_clear(nzdpoly, ctx);
    flint_randclear(state);
    flint_free(Adegs);

#if WANT_ASSERT
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


static int nmod_mpoly_factor_algo(
    nmod_mpoly_factor_t f,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong i, j;
    nmod_mpolyv_t t;
    nmod_mpoly_factor_t g;

    nmod_mpolyv_init(t, ctx);
    nmod_mpoly_factor_init(g, ctx);

    success = nmod_mpoly_factor_squarefree(f, A, ctx);
    if (!success)
        goto cleanup;

    g->constant = f->constant;
    g->num = 0;
    for (j = 0; j < f->num; j++)
    {
        success = _irreducible_factors(t, f->poly + j, ctx, algo);
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

    FLINT_ASSERT(!success || nmod_mpoly_factor_matches(A, f, ctx));

    return success;
}



int nmod_mpoly_factor(
    nmod_mpoly_factor_t f,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    return nmod_mpoly_factor_algo(f, A, ctx, USE_ZAS | USE_WANG | USE_ZIP);
}


int nmod_mpoly_factor_zassenhaus(
    nmod_mpoly_factor_t f,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    return nmod_mpoly_factor_algo(f, A, ctx, USE_ZAS);
}


int nmod_mpoly_factor_wang(
    nmod_mpoly_factor_t f,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    return nmod_mpoly_factor_algo(f, A, ctx, USE_WANG);
}


int nmod_mpoly_factor_zippel(
    nmod_mpoly_factor_t f,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    return nmod_mpoly_factor_algo(f, A, ctx, USE_ZIP);
}

