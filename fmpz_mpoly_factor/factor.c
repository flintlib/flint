/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"
#include "fmpq_poly.h"
#include "fmpz_mod_mpoly.h"
#include "nmod_mpoly_factor.h"

#define USE_ZAS 1
#define USE_WANG 2
#define USE_ZIP 4



void fmpz_mpoly_convert_perm(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t Actx,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t Bctx,
    const slong * perm)
{
    slong n = Bctx->minfo->nvars;
    slong m = Actx->minfo->nvars;
    slong i, k, l;
    slong NA, NB;
    ulong * Aexps;
    ulong * Bexps;
    TMP_INIT;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(B->length > 0);
    TMP_START;

    Aexps = (ulong *) TMP_ALLOC(m*sizeof(ulong));
    Bexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, Actx->minfo);
    NB = mpoly_words_per_exp(B->bits, Bctx->minfo);

    fmpz_mpoly_fit_length_set_bits(A, B->length, Abits, Actx);
    A->length = B->length;
    for (i = 0; i < B->length; i++)
    {        
        fmpz_set(A->coeffs + i, B->coeffs + i);
        mpoly_get_monomial_ui(Bexps, B->exps + NB*i, B->bits, Bctx->minfo);
        for (k = 0; k < m; k++)
        {
            l = perm[k];
            Aexps[k] = l < 0 ? 0 : Bexps[l];
        }
        mpoly_set_monomial_ui(A->exps + NA*i, Aexps, Abits, Actx->minfo);
     }  
    TMP_END;
    fmpz_mpoly_sort_terms(A, Actx);
}


/*
    A is primitive w.r.t to any variable appearing in A.
    A is square free with positive lead coeff.

    return 1 for success, 0 for failure
*/
static int _irreducible_factors(
    fmpz_mpolyv_t Af,
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong i, j, mvars;
    slong * Adegs, * perm, * iperm;
    slong nvars = ctx->minfo->nvars;
    flint_bitcnt_t Lbits, Abits;
    double density;
    int perm_is_id;
    fmpz_poly_t Au;
    fmpz_poly_factor_t Auf;
    flint_rand_t state;
#if WANT_ASSERT
    fmpz_mpoly_t Aorg;

    fmpz_mpoly_init(Aorg, ctx);
    fmpz_mpoly_set(Aorg, A, ctx);
#endif

    FLINT_ASSERT(A->length == 0 || fmpz_sgn(A->coeffs + 0) > 0);

    flint_randinit(state);
    Adegs = (slong *) flint_malloc(3*nvars*sizeof(slong));
    perm = Adegs + nvars;
    iperm = perm + nvars;
    fmpz_poly_init(Au);
    fmpz_poly_factor_init(Auf);

    if (A->length < 2)
    {
        FLINT_ASSERT(A->length == 1);
        FLINT_ASSERT(!fmpz_mpoly_is_fmpz(A, ctx));
        fmpz_mpolyv_fit_length(Af, 1, ctx);
        Af->length = 1;
        fmpz_mpoly_swap(Af->coeffs + 0, A, ctx);
        success = 1;
        goto cleanup;
    }

    if (!fmpz_mpoly_degrees_fit_si(A, ctx))
    {
        success = 0;
        goto cleanup;
    }

    if (A->bits > FLINT_BITS &&
        !fmpz_mpoly_repack_bits_inplace(A, FLINT_BITS, ctx))
    {
        success = 0;
        goto cleanup;
    }

    fmpz_mpoly_degrees_si(Adegs, A, ctx);

    Abits = A->bits;

    density = A->length;

    mvars = 0;
    Lbits = 0;
    for (i = 0; i < nvars; i++)
    {
        iperm[i] = -1;
        if (Adegs[i] > 0)
        {
            flint_bitcnt_t this_bits = FLINT_BIT_COUNT(Adegs[i]);
            density /= Adegs[i] + 1;
            Lbits = FLINT_MAX(Lbits, this_bits);
            perm[mvars] = i;
            mvars++;
        }
    }

    if (Lbits > FLINT_BITS - 10)
    {
        success = 0;
        goto cleanup;
    }

    /* sort by degree */
    for (i = 1; i < mvars; i++)
    {
        for (j = i; j > 0 && Adegs[perm[j]] < Adegs[perm[j - 1]]; j--)
        {
            SLONG_SWAP(perm[j], perm[j - 1]);
        }
    }

    /* invert perm */
    perm_is_id = (mvars == nvars);
    for (i = 0; i < mvars; i++)
    {
        perm_is_id = perm_is_id && (perm[i] == i);
        iperm[perm[i]] = i;
    }

    if (mvars < 2)
    {
        FLINT_ASSERT(mvars == 1);

        FLINT_ASSERT(fmpz_mpoly_is_fmpz_poly(A, perm[0], ctx));
        success = fmpz_mpoly_get_fmpz_poly(Au, A, perm[0], ctx);
        FLINT_ASSERT(success);
        fmpz_poly_factor(Auf, Au);

        /* c = +-1 */
        FLINT_ASSERT(fmpz_is_pm1(&Auf->c));

        fmpz_mpolyv_fit_length(Af, Auf->num, ctx);
        Af->length = Auf->num;
        for (i = 0; i < Auf->num; i++)
        {
            FLINT_ASSERT(Auf->exp[i] == 1);
            _fmpz_mpoly_set_fmpz_poly(Af->coeffs + i, Abits,
                             Auf->p[i].coeffs, Auf->p[i].length, perm[0], ctx);
        }

        success = 1;
    }
    else if (mvars == 2)
    {
        fmpz_poly_struct * c = Au;
        fmpz_bpoly_t Ab;
        fmpz_tpoly_t Abf;

        fmpz_bpoly_init(Ab);
        fmpz_tpoly_init(Abf);

        fmpz_mpoly_get_bpoly(Ab, A, perm[0], perm[1], ctx);
        fmpz_bpoly_factor(c, Abf, Ab);

        FLINT_ASSERT(c->length == 1 && fmpz_is_pm1(c->coeffs + 0));

        fmpz_mpolyv_fit_length(Af, Abf->length, ctx);
        Af->length = Abf->length;
        for (i = 0; i < Abf->length; i++)
        {
            fmpz_mpoly_set_fmpz_bpoly(Af->coeffs + i, Abits, Abf->coeffs + i,
                                                        perm[0], perm[1], ctx);
            fmpz_mpoly_unit_normalize(Af->coeffs + i, ctx);
        }

        fmpz_bpoly_clear(Ab);
        fmpz_tpoly_clear(Abf);

        success = 1;
    }
    else
    {
        fmpz_mpoly_ctx_t Lctx;
        fmpz_mpoly_t Lcopy, lcL;
        fmpz_mpolyv_t Lf;
        fmpz_mpoly_factor_t lcLf;
        fmpz_mpoly_struct * L;
        zassenhaus_prune_t Z;
        fmpz * alpha;
        int zero_ok, trying_zero, image_count, sqrfree;
        ulong alpha_modulus;
        slong min_uni_r, max_uni_deg, max_deg;
        slong main_degree = Adegs[perm[0]];

        fmpz_mpoly_ctx_init(Lctx, mvars, ORD_LEX);
        fmpz_mpoly_init(Lcopy, Lctx);
        fmpz_mpoly_init(lcL, Lctx);
        fmpz_mpolyv_init(Lf, Lctx);
        fmpz_mpoly_factor_init(lcLf, Lctx);
        zassenhaus_prune_init(Z);

        zassenhaus_prune_set_degree(Z, main_degree);

        Lbits = mpoly_fix_bits(Lbits + 1, Lctx->minfo);

        perm_is_id = perm_is_id && Lbits == Abits && ctx->minfo->ord == ORD_LEX;
        if (perm_is_id)
        {
            L = A;
        }
        else
        {
            fmpz_mpoly_convert_perm(Lcopy, Lbits, Lctx, A, ctx, perm);
            fmpz_mpoly_unit_normalize(Lcopy, ctx);
            L = Lcopy;
        }

        /* some simple checks */

        alpha = _fmpz_vec_init(mvars - 1);
        zero_ok = 0;
        trying_zero = 1;
        alpha_modulus = 5;
        min_uni_r = WORD_MAX;
        max_uni_deg = 0;
        image_count = 0;

        goto got_alpha;

next_alpha:

        trying_zero = 0;

        if (++alpha_modulus > 10)
            goto done_alpha;

        for (i = 0; i < mvars - 1; i++)
        {
            slong a = n_urandint(state, alpha_modulus);
            a -= alpha_modulus/2;
            fmpz_set_si(alpha + i, a);
        }

got_alpha:

        _fmpz_mpoly_eval_rest_to_poly(Au, L, alpha, Lctx);

        if (fmpz_poly_degree(Au) != main_degree)
            goto next_alpha;

        fmpz_poly_factor(Auf, Au);

        zassenhaus_prune_start_add_factors(Z);                
        sqrfree = 1;
        max_deg = 0;
        for (i = 0; i < Auf->num; i++)
        {
            if (Auf->exp[i] != 1)
                sqrfree = 0;
            max_deg = FLINT_MAX(max_deg, fmpz_poly_degree(Auf->p + i));
            zassenhaus_prune_add_factor(Z, fmpz_poly_degree(Auf->p + i), Auf->exp[i]);
        }
        zassenhaus_prune_end_add_factors(Z);

        if (!sqrfree)
            goto next_alpha;

        if (Auf->num < min_uni_r)
        {
            min_uni_r = Auf->num;
            max_uni_deg = max_deg;
        }

        zero_ok = zero_ok || trying_zero;
        if (++image_count < 3)
            goto next_alpha;

done_alpha:

        _fmpz_vec_clear(alpha, mvars - 1);

        /* simple check done */

        if (zassenhaus_prune_must_be_irreducible(Z))
        {
            fmpz_mpolyv_fit_length(Af, 1, ctx);
            Af->length = 1;
            fmpz_mpoly_swap(Af->coeffs + 0, A, ctx);
            success = 1;
            goto cleanup2;
        }

        success = 0;

        if (algo & (USE_WANG | USE_ZIP))
        {
            _fmpz_mpoly_get_lead0(lcL, L, Lctx);

            if (fmpz_mpoly_factor(lcLf, lcL, Lctx))
            {
                if (!(algo & USE_ZIP))
                {
                    success = fmpz_mpoly_factor_irred_wang(Lf, L,
                                                 lcLf, lcL, Lctx, state, Z, 1);
                }
                else if (!(algo & USE_WANG))
                {
                    success = fmpz_mpoly_factor_irred_zippel(Lf, L,
                                                    lcLf, lcL, Lctx, state, Z);
                }
                else
                {
                    if (density > 0.002 && zero_ok && max_uni_deg < 40)
                    {
                        success = fmpz_mpoly_factor_irred_wang(Lf, L,
                                                 lcLf, lcL, Lctx, state, Z, 0);
                    }

                    if (success == 0 && density > 0.04 && max_uni_deg < 20)
                    {
                        success = fmpz_mpoly_factor_irred_wang(Lf, L,
                                                 lcLf, lcL, Lctx, state, Z, 1);
                    }

                    if (success == 0)
                    {
                        success = fmpz_mpoly_factor_irred_zippel(Lf, L,
                                                    lcLf, lcL, Lctx, state, Z);
                    }

                    if (success == 0)
                    {
                        success = fmpz_mpoly_factor_irred_wang(Lf, L,
                                                 lcLf, lcL, Lctx, state, Z, 1);
                    }
                }
            }
        }

        if (algo & USE_ZAS)
        {
            if (success == 0)
                success = fmpz_mpoly_factor_irred_zassenhaus(Lf, L, Lctx, Z);
        }

        success = (success > 0);
        if (success)
        {
            fmpz_mpolyv_fit_length(Af, Lf->length, ctx);
            Af->length = Lf->length;
            for (i = 0; i < Lf->length; i++)
            {
                if (perm_is_id)
                {
                    fmpz_mpoly_swap(Af->coeffs + i, Lf->coeffs + i, Lctx);
                }
                else
                {
                    fmpz_mpoly_convert_perm(Af->coeffs + i, Abits, ctx,
                                                  Lf->coeffs + i, Lctx, iperm);
                }
                fmpz_mpoly_unit_normalize(Af->coeffs + i, ctx);
            }
        }

cleanup2:

        fmpz_mpoly_clear(Lcopy, Lctx);
        fmpz_mpoly_clear(lcL, Lctx);
        fmpz_mpolyv_clear(Lf, Lctx);
        fmpz_mpoly_factor_clear(lcLf, Lctx);
        fmpz_mpoly_ctx_clear(Lctx);

        zassenhaus_prune_clear(Z);
    }

cleanup:

    fmpz_poly_clear(Au);
    fmpz_poly_factor_clear(Auf);

    flint_randclear(state);
    flint_free(Adegs);

#if WANT_ASSERT
    if (success)
    {
        fmpz_mpoly_t prod;
        fmpz_mpoly_init(prod, ctx);
        fmpz_mpoly_one(prod, ctx);
        for (i = 0; i < Af->length; i++)
            fmpz_mpoly_mul(prod, prod, Af->coeffs + i, ctx);
        FLINT_ASSERT(fmpz_mpoly_equal(prod, Aorg, ctx));
        fmpz_mpoly_clear(prod, ctx);
        fmpz_mpoly_clear(Aorg, ctx);
    }
#endif

    FLINT_ASSERT(success == 0 || success == 1);
    return success;
}


static int fmpz_mpoly_factor_algo(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong i, j;
    fmpz_mpolyv_t t;
    fmpz_mpoly_factor_t g;

    fmpz_mpolyv_init(t, ctx);
    fmpz_mpoly_factor_init(g, ctx);

    success = fmpz_mpoly_factor_squarefree(f, A, ctx);
    if (!success)
        goto cleanup;

    FLINT_ASSERT(fmpz_mpoly_factor_matches(A, f, ctx));

    fmpz_swap(g->constant, f->constant);
    g->num = 0;
    for (j = 0; j < f->num; j++)
    {
        success = _irreducible_factors(t, f->poly + j, ctx, algo);
        if (!success)
            goto cleanup;

        fmpz_mpoly_factor_fit_length(g, g->num + t->length, ctx);
        for (i = 0; i < t->length; i++)
        {
            fmpz_set(g->exp + g->num, f->exp + j);
            fmpz_mpoly_swap(g->poly + g->num, t->coeffs + i, ctx);
            g->num++;
        }
    }
    fmpz_mpoly_factor_swap(f, g, ctx);

    success = 1;

cleanup:

    fmpz_mpolyv_clear(t, ctx);
    fmpz_mpoly_factor_clear(g, ctx);

    FLINT_ASSERT(!success || fmpz_mpoly_factor_matches(A, f, ctx));

    return success;
}


int fmpz_mpoly_factor(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    return fmpz_mpoly_factor_algo(f, A, ctx, USE_ZAS | USE_WANG | USE_ZIP);
}


int fmpz_mpoly_factor_zassenhaus(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    return fmpz_mpoly_factor_algo(f, A, ctx, USE_ZAS);
}


int fmpz_mpoly_factor_wang(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    return fmpz_mpoly_factor_algo(f, A, ctx, USE_WANG);
}


int fmpz_mpoly_factor_zippel(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    return fmpz_mpoly_factor_algo(f, A, ctx, USE_ZIP);
}

