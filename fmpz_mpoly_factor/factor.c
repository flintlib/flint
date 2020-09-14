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


typedef struct {
    fmpz_mpoly_ctx_t ctx;
    slong nvars;
    slong exps_alloc;
    slong * exps;
    slong * umat;
    slong * deltas;
    slong * degs;
    int is_id;
} fmpz_mpoly_perm_struct;

typedef fmpz_mpoly_perm_struct fmpz_mpoly_perm_t[1];


static void fmpz_mpoly_perm_clear(
    fmpz_mpoly_perm_t P,
    fmpz_mpoly_t L)
{
    flint_free(P->umat);
    flint_free(P->deltas);
    flint_free(P->degs);
    flint_free(P->exps);

    fmpz_mpoly_clear(L, P->ctx);    
    fmpz_mpoly_ctx_clear(P->ctx);
}


static void fmpz_mpoly_perm_init(
    fmpz_mpoly_perm_t P,
    fmpz_mpoly_t L,
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t Actx)
{
    int same;
    slong i, j, max_deg;
    flint_bitcnt_t Lbits, Abits = A->bits;
    slong LN, AN = mpoly_words_per_exp_sp(Abits, Actx->minfo);
    slong Alen = A->length;
    slong mvars, nvars = Actx->minfo->nvars;

    P->nvars = nvars;
    P->umat = FLINT_ARRAY_ALLOC(nvars*nvars, slong);
    P->deltas = FLINT_ARRAY_ALLOC(nvars, slong);
    P->degs = FLINT_ARRAY_ALLOC(nvars, slong);
    P->exps = FLINT_ARRAY_ALLOC(Alen*nvars, slong);
    P->exps_alloc = Alen*nvars;
    for (i = 0; i < Alen; i++)
        mpoly_get_monomial_ui_sp((ulong *)P->exps + nvars*i,
                                           A->exps + AN*i, Abits, Actx->minfo);
    mvars = _mpoly_compress_exps(P->umat, P->deltas, P->degs, P->exps, nvars, Alen);

    FLINT_ASSERT(mvars > 0);

    fmpz_mpoly_ctx_init(P->ctx, mvars, ORD_LEX);
    fmpz_mpoly_init(L, P->ctx);

    same = (mvars == nvars) && (Actx->minfo->ord == ORD_LEX);
    if (same)
    {
        for (i = 0; i < nvars; i++)
        {
            same = same && (P->deltas[i] == 0);
            for (j = 0; j < nvars; j++)
                same = same && (P->umat[i*nvars + j] == (i == j));
        }
    }

    P->is_id = same;
    if (same)
    {
        fmpz_mpoly_swap(L, A, Actx);
        return;
    }

    max_deg = P->degs[0];
    for (i = 1; i < mvars; i++)
        max_deg = FLINT_MAX(max_deg, P->degs[i]);
    Lbits = mpoly_fix_bits(1 + FLINT_BIT_COUNT(max_deg), P->ctx->minfo);

    fmpz_mpoly_fit_length_set_bits(L, Alen, Lbits, P->ctx);

    LN = mpoly_words_per_exp(Lbits, P->ctx->minfo);

    L->length = Alen;
    for (i = 0; i < Alen; i++)
    {
        fmpz_swap(L->coeffs + i, A->coeffs + i);
        mpoly_set_monomial_ui(L->exps + LN*i,
                             (ulong *)P->exps + nvars*i, Lbits, P->ctx->minfo);
    }

    fmpz_mpoly_sort_terms(L, P->ctx);
    fmpz_mpoly_unit_normalize(L, P->ctx);
}


static void _slong_array_fit_length(slong ** array, slong * alloc, slong len)
{
    if (len <= *alloc)
        return;
    len = FLINT_MAX(len, *alloc + *alloc/2 + 1);
    *array = flint_realloc(*array, len*sizeof(slong));
    *alloc = len;
}


static void fmpz_mpoly_perm_expand_fmpz_mpoly(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t Actx,
    fmpz_mpoly_t B,
    fmpz_mpoly_perm_t P)
{
    slong i, k, l;
    slong nvars = Actx->minfo->nvars;
    slong NA = mpoly_words_per_exp(Abits, Actx->minfo);
    slong mvars = P->ctx->minfo->nvars;
    flint_bitcnt_t Bbits = B->bits;
    slong NB = mpoly_words_per_exp(Bbits, P->ctx->minfo);
    slong * mins, * texps;
    TMP_INIT;

    FLINT_ASSERT(fmpz_mpoly_degrees_fit_si(B, P->ctx));

    if (P->is_id)
    {
        fmpz_mpoly_swap(A, B, Actx);
        return;
    }

    TMP_START;
    texps = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    mins = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    for (k = 0; k < nvars; k++)
        mins[k] = WORD_MAX;

    _slong_array_fit_length(&P->exps, &P->exps_alloc, B->length*nvars);
    fmpz_mpoly_fit_length_set_bits(A, B->length, Abits, Actx);
    _fmpz_mpoly_set_length(A, B->length, Actx);
    for (i = 0; i < B->length; i++)
    {
        fmpz_swap(A->coeffs + i, B->coeffs + i);
        mpoly_get_monomial_ui((ulong *)texps, B->exps + NB*i, Bbits, P->ctx->minfo);
        for (k = 0; k < nvars; k++)
        {
            slong tot = P->deltas[k];
            for (l = 0; l < mvars; l++)
                tot += P->umat[k*nvars + l]*texps[l];
            P->exps[i*nvars + k] = tot;
            mins[k] = FLINT_MIN(mins[k], tot);
        }
    }

    for (i = 0; i < B->length; i++)
    {
        for (k = 0; k < nvars; k++)
            P->exps[i*nvars + k] -= mins[k];
        mpoly_set_monomial_ui(A->exps + NA*i, (ulong *)P->exps + i*nvars,
                                                           Abits, Actx->minfo);
    }

    TMP_END;

    fmpz_mpoly_sort_terms(A, Actx);
    fmpz_mpoly_unit_normalize(A, Actx);
}

/* A has degree 2 wrt gen(0) */
static void _apply_quadratic(
    fmpz_mpolyv_t Af,
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong shift, off, N;
    flint_bitcnt_t bits = A->bits;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    fmpz_mpoly_t a_mock, b_mock, c_mock;
    fmpz_mpoly_t t0, t1, t2, t3;

    FLINT_ASSERT(A->length > 1 || fmpz_sgn(A->coeffs + 0) > 0);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);

    fmpz_mpoly_init(t0, ctx);
    fmpz_mpoly_init(t1, ctx);
    fmpz_mpoly_init(t2, ctx);
    fmpz_mpoly_init(t3, ctx);

    mpoly_gen_offset_shift_sp(&off, &shift, 0, bits, ctx->minfo);
    N = mpoly_words_per_exp_sp(bits, ctx->minfo);

    i = 0;
    a_mock->exps = A->exps + N*i;
    a_mock->coeffs = A->coeffs + i;
    a_mock->bits = bits;
    while (i < A->length && (mask & (A->exps[N*i + off] >> shift)) == 2)
        i++;
    a_mock->length = i;
    a_mock->alloc = a_mock->length;

    b_mock->exps = A->exps + N*i;
    b_mock->coeffs = A->coeffs + i;
    b_mock->bits = bits;
    while (i < A->length && (mask & (A->exps[N*i + off] >> shift)) == 1)
        i++;
    b_mock->length = i - a_mock->length;
    b_mock->alloc = b_mock->length;

    c_mock->exps = A->exps + N*i;
    c_mock->coeffs = A->coeffs + i;
    c_mock->bits = bits;
    c_mock->length = A->length - i;
    c_mock->alloc = c_mock->length;

    FLINT_ASSERT(a_mock->length > 0);
    FLINT_ASSERT(c_mock->length > 0);

    fmpz_mpoly_mul(t0, b_mock, b_mock, ctx);
    fmpz_mpoly_mul(t1, a_mock, c_mock, ctx);
    fmpz_mpoly_scalar_mul_si(t1, t1, 4, ctx);
    fmpz_mpoly_sub(t2, t0, t1, ctx);
    if (!fmpz_mpoly_sqrt(t0, t2, ctx))
    {
        fmpz_mpolyv_fit_length(Af, 1, ctx);
        Af->length = 1;
        fmpz_mpoly_swap(Af->coeffs + 0, A, ctx);
        goto cleanup;
    }

    fmpz_mpoly_add(t2, t0, b_mock, ctx);
    fmpz_mpoly_scalar_divides_si(t2, t2, 2, ctx);
    fmpz_mpoly_gcd_cofactors(t0, t1, t2, a_mock, t2, ctx);
    fmpz_mpoly_divides(t3, c_mock, t2, ctx);

    fmpz_mpolyv_fit_length(Af, 2, ctx);
    Af->length = 2;
    fmpz_mpoly_add(Af->coeffs + 0, t1, t2, ctx);
    fmpz_mpoly_add(Af->coeffs + 1, t0, t3, ctx);

cleanup:

    fmpz_mpoly_clear(t0, ctx);
    fmpz_mpoly_clear(t1, ctx);
    fmpz_mpoly_clear(t2, ctx);
    fmpz_mpoly_clear(t3, ctx);
}


/* A is squarefree and primitive wrt gen(0) */
static int _apply_algo(
    fmpz_mpolyv_t Af,
    fmpz_mpoly_t A,
    const slong * Adegs,
    const fmpz_mpoly_ctx_t ctx,
    flint_rand_t state,
    unsigned int algo)
{
    int success;
    slong i;
    flint_bitcnt_t Abits;
    fmpz_poly_t u;
    fmpz_poly_factor_t uf;
    fmpz_mpoly_t lcA;
    fmpz_mpoly_factor_t lcAf;
    zassenhaus_prune_t Z;

    FLINT_ASSERT(A->length > 1 || fmpz_sgn(A->coeffs + 0) > 0);
    FLINT_ASSERT(fmpz_mpoly_degrees_fit_si(A, ctx));
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(fmpz_mpoly_is_canonical(A, ctx));

    if (A->bits > FLINT_BITS && !fmpz_mpoly_repack_bits_inplace(A, FLINT_BITS, ctx))
        return 0;

    Abits = A->bits;

    flint_randinit(state);
    fmpz_poly_init(u);
    fmpz_poly_factor_init(uf);
    fmpz_mpoly_init(lcA, ctx);
    fmpz_mpoly_factor_init(lcAf, ctx);
    zassenhaus_prune_init(Z);

    /* init done */

    if (Adegs[0] == 1)
    {
        fmpz_mpolyv_fit_length(Af, 1, ctx);
        Af->length = 1;
        fmpz_mpoly_swap(Af->coeffs + 0, A, ctx);
        success = 1;
        goto cleanup;
    }
    else if (Adegs[0] == 2)
    {
        _apply_quadratic(Af, A, ctx);
        success = 1;
        goto cleanup;
    }

    if (ctx->minfo->nvars < 2)
    {
        fmpz_mpoly_get_fmpz_poly(u, A, 0, ctx);
        fmpz_poly_factor(uf, u);

        FLINT_ASSERT(fmpz_is_pm1(&uf->c));

        fmpz_mpolyv_fit_length(Af, uf->num, ctx);
        Af->length = uf->num;
        for (i = 0; i < uf->num; i++)
        {
            FLINT_ASSERT(uf->exp[i] == 1);
            _fmpz_mpoly_set_fmpz_poly(Af->coeffs + i, Abits,
                                     uf->p[i].coeffs, uf->p[i].length, 0, ctx);
        }

        success = 1;
    }
    else if (ctx->minfo->nvars == 2)
    {
        fmpz_bpoly_t b;
        fmpz_tpoly_t bf;

        fmpz_bpoly_init(b);
        fmpz_tpoly_init(bf);

        fmpz_mpoly_get_bpoly(b, A, 0, 1, ctx);
        fmpz_bpoly_factor(u, bf, b);

        FLINT_ASSERT(u->length == 1 && fmpz_is_pm1(u->coeffs + 0));

        fmpz_mpolyv_fit_length(Af, bf->length, ctx);
        Af->length = bf->length;
        for (i = 0; i < bf->length; i++)
        {
            fmpz_mpoly_set_fmpz_bpoly(Af->coeffs + i, Abits,
                                                    bf->coeffs + i, 0, 1, ctx);
        }

        fmpz_bpoly_clear(b);
        fmpz_tpoly_clear(bf);

        success = 1;
    }
    else
    {
        int zero_ok, trying_zero, image_count, sqrfree;
        ulong alpha_modulus;
        double density;
        fmpz * alpha;

        zassenhaus_prune_set_degree(Z, Adegs[0]);

        /* some simple checks */

        alpha = _fmpz_vec_init(ctx->minfo->nvars - 1);
        zero_ok = 0;
        trying_zero = 1;
        alpha_modulus = 5;
        image_count = 0;

        goto got_alpha;

next_alpha:

        trying_zero = 0;

        if (++alpha_modulus > 10)
            goto done_alpha;

        for (i = 0; i < ctx->minfo->nvars - 1; i++)
        {
            slong a = n_urandint(state, alpha_modulus);
            a -= alpha_modulus/2;
            fmpz_set_si(alpha + i, a);
        }

got_alpha:

        _fmpz_mpoly_eval_rest_to_poly(u, A, alpha, ctx);

        if (fmpz_poly_degree(u) != Adegs[0])
            goto next_alpha;

        fmpz_poly_factor(uf, u);

        zassenhaus_prune_start_add_factors(Z);                
        sqrfree = 1;
        for (i = 0; i < uf->num; i++)
        {
            if (uf->exp[i] != 1)
                sqrfree = 0;
            zassenhaus_prune_add_factor(Z, fmpz_poly_degree(uf->p + i), uf->exp[i]);
        }
        zassenhaus_prune_end_add_factors(Z);

        if (!sqrfree)
            goto next_alpha;

        zero_ok = zero_ok || trying_zero;
        if (++image_count < 3)
            goto next_alpha;

done_alpha:

        _fmpz_vec_clear(alpha, ctx->minfo->nvars - 1);

        /* simple check done */

        if (zassenhaus_prune_must_be_irreducible(Z))
        {
            fmpz_mpolyv_fit_length(Af, 1, ctx);
            Af->length = 1;
            fmpz_mpoly_swap(Af->coeffs + 0, A, ctx);
            success = 1;
            goto cleanup;
        }

        density = A->length;
        for (i = 0; i < ctx->minfo->nvars; i++)
            density /= Adegs[i] + 1;

        success = 0;

        if (algo & (USE_WANG | USE_ZIP))
        {
            _fmpz_mpoly_get_lead0(lcA, A, ctx);

            if (fmpz_mpoly_factor_squarefree(lcAf, lcA, ctx))
            {
                int irr_fac = 1;
                for (i = 0; i < lcAf->num; i++)
                    irr_fac = irr_fac && lcAf->poly[i].length < 4;

                irr_fac = irr_fac && fmpz_mpoly_factor_irred(lcAf, ctx, algo);

                if (!(algo & USE_ZIP))
                {
                    success = fmpz_mpoly_factor_irred_wang(Af, A,
                                         lcAf, irr_fac, lcA, ctx, state, Z, 1);
                }
                else if (!(algo & USE_WANG))
                {
                    success = fmpz_mpoly_factor_irred_zippel(Af, A,
                                            lcAf, irr_fac, lcA, ctx, state, Z);
                    FLINT_ASSERT(success);
                }
                else
                {
                    if (density > 0.002 && zero_ok)
                    {
                        success = fmpz_mpoly_factor_irred_wang(Af, A,
                                         lcAf, irr_fac, lcA, ctx, state, Z, 0);
                    }

                    if (success == 0 && density > 0.04)
                    {
                        success = fmpz_mpoly_factor_irred_wang(Af, A,
                                         lcAf, irr_fac, lcA, ctx, state, Z, 1);
                    }

                    if (success == 0)
                    {
                        success = fmpz_mpoly_factor_irred_zippel(Af, A,
                                            lcAf, irr_fac, lcA, ctx, state, Z);
                    }

                    if (success == 0)
                    {
                        success = fmpz_mpoly_factor_irred_wang(Af, A,
                                         lcAf, irr_fac, lcA, ctx, state, Z, 1);
                    }
                }
            }
        }

        if (algo & USE_ZAS)
        {
            if (success == 0)
                success = fmpz_mpoly_factor_irred_zassenhaus(Af, A, ctx, Z);
        }

        success = (success > 0);
    }

cleanup:

    fmpz_poly_clear(u);
    fmpz_poly_factor_clear(uf);
    fmpz_mpoly_clear(lcA, ctx);
    fmpz_mpoly_factor_clear(lcAf, ctx);
    zassenhaus_prune_clear(Z);

    FLINT_ASSERT(success == 0 || success == 1);
    return success;
}


/* replace f by a content factorization, f is already squarefree */
int _refine_content_factors(
    fmpz_mpolyv_t f,
    fmpz_mpolyv_t g,        /* temp */
    flint_bitcnt_t bits,
    fmpz_mpoly_univar_t u,  /* temp */
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong v, i, j;
    fmpz_mpoly_struct * c;

    for (v = 1; v < ctx->minfo->nvars; v++)
    {
        g->length = 0;
        for (j = 0; j < f->length; j++)
        {
            fmpz_mpoly_to_univar(u, f->coeffs + j, v, ctx);
            FLINT_ASSERT(u->length > 0);
            FLINT_ASSERT(fmpz_is_zero(u->exps + u->length - 1));

            fmpz_mpolyv_fit_length(g, g->length + 1, ctx);
            c = g->coeffs + g->length;
            success = _fmpz_mpoly_vec_content_mpoly(c, u->coeffs, u->length, ctx);
            if (!success)
                return 0;

            if (fmpz_mpoly_is_fmpz(c, ctx))
            {
                fmpz_mpoly_swap(c, f->coeffs + j, ctx);
                g->length++;
            }
            else
            {
                for (i = 0; i < u->length; i++)
                {
                    success = fmpz_mpoly_divides(u->coeffs + i, u->coeffs + i, c, ctx);
                    FLINT_ASSERT(success);
                }

                g->length++;

                if (u->length > 1)
                {
                    fmpz_mpolyv_fit_length(g, g->length + 1, ctx);
                    c = g->coeffs + g->length;
                    fmpz_mpoly_from_univar_bits(c, bits, u, v, ctx);
                    g->length++;
                }
                else
                {
                    FLINT_ASSERT(fmpz_mpoly_is_one(u->coeffs + 0, ctx));
                }
            }
        }

        fmpz_mpolyv_swap(f, g, ctx);
    }

    return 1;
}

/*
    A is primitive w.r.t to any variable appearing in A.
    A is square free with positive lead coeff.

    return 1 for success, 0 for failure
*/
static int _irreducible_factors(
    fmpz_mpolyv_t Af,
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t Actx,
    unsigned int algo)
{
    int success;
    slong i, j;
    flint_bitcnt_t Abits;
    fmpz_mpoly_perm_t P;
    fmpz_mpoly_t L;
    fmpz_mpolyv_t Lf;
    fmpz_mpoly_ctx_struct *  Lctx;
    flint_rand_t state;
#if WANT_ASSERT
    fmpz_mpoly_t Aorg;

    fmpz_mpoly_init(Aorg, Actx);
    fmpz_mpoly_set(Aorg, A, Actx);
#endif

    FLINT_ASSERT(A->length == 0 || fmpz_sgn(A->coeffs + 0) > 0);

    if (A->length < 2)
    {
        FLINT_ASSERT(A->length == 1);
        FLINT_ASSERT(!fmpz_mpoly_is_fmpz(A, Actx));
        fmpz_mpolyv_fit_length(Af, 1, Actx);
        Af->length = 1;
        fmpz_mpoly_swap(Af->coeffs + 0, A, Actx);
        success = 1;
        goto cleanup_less;
    }

    if (!fmpz_mpoly_degrees_fit_si(A, Actx))
    {
        success = 0;
        goto cleanup_less;
    }

    if (A->bits > FLINT_BITS &&
        !fmpz_mpoly_repack_bits_inplace(A, FLINT_BITS, Actx))
    {
        success = 0;
        goto cleanup_less;
    }

    Abits = A->bits;

    flint_randinit(state);
    fmpz_mpoly_perm_init(P, L, A, Actx);
    Lctx = P->ctx;
    fmpz_mpolyv_init(Lf, Lctx);

    if (!P->is_id)
    {
        fmpz_mpoly_univar_t U;
        fmpz_mpoly_t t;
        fmpz_mpolyv_t tf, sf;

        fmpz_mpoly_init(t, Lctx);
        fmpz_mpoly_univar_init(U, Lctx);
        fmpz_mpolyv_init(tf, Lctx);
        fmpz_mpolyv_init(sf, Lctx);

        fmpz_mpoly_to_univar(U, L, 0, Lctx);
        success = _fmpz_mpoly_vec_content_mpoly(t, U->coeffs, U->length, Lctx);
        if (!success)
            goto cleanup_more;

        if (!fmpz_mpoly_is_fmpz(t, Lctx))
        {
            success = fmpz_mpoly_divides(L, L, t, Lctx);
            FLINT_ASSERT(success);
            fmpz_mpoly_unit_normalize(L, Lctx);

            fmpz_mpolyv_fit_length(sf, 2, Lctx);
            sf->length = 2;
            fmpz_mpoly_swap(sf->coeffs + 0, L, Lctx);
            fmpz_mpoly_swap(sf->coeffs + 1, t, Lctx);

            success = _refine_content_factors(sf, tf, Abits, U, Lctx);
            if (!success)
                goto cleanup_more;

            Lf->length = 0;
            for (i = 0; i < sf->length; i++)
            {
                success = _irreducible_factors(tf, sf->coeffs + i, Lctx, algo);
                if (!success)
                    goto cleanup_more;

                fmpz_mpolyv_fit_length(Lf, Lf->length + tf->length, Lctx);
                for (j = 0; j < tf->length; j++)
                    fmpz_mpoly_swap(Lf->coeffs + Lf->length++, tf->coeffs + j, Lctx);
            }
        }
        else
        {
            success = _apply_algo(Lf, L, P->degs, Lctx, state, algo);
        }

cleanup_more:

        fmpz_mpoly_clear(t, Lctx);
        fmpz_mpoly_univar_clear(U, Lctx);
        fmpz_mpolyv_clear(tf, Lctx);
        fmpz_mpolyv_clear(sf, Lctx);
    }
    else
    {
        success = _apply_algo(Lf, L, P->degs, Lctx, state, algo);
    }

    if (success)
    {
        fmpz_mpolyv_fit_length(Af, Lf->length, Actx);
        Af->length = Lf->length;
        for (i = 0; i < Lf->length; i++)
            fmpz_mpoly_perm_expand_fmpz_mpoly(Af->coeffs + i, Abits, Actx,
                                                            Lf->coeffs + i, P);
    }

    flint_randclear(state);
    fmpz_mpolyv_clear(Lf, Lctx);
    fmpz_mpoly_perm_clear(P, L);

cleanup_less:

#if WANT_ASSERT
    if (success)
    {
        fmpz_mpoly_t prod;
        fmpz_mpoly_init(prod, Actx);
        fmpz_mpoly_one(prod, Actx);
        for (i = 0; i < Af->length; i++)
            fmpz_mpoly_mul(prod, prod, Af->coeffs + i, Actx);
        FLINT_ASSERT(fmpz_mpoly_equal(prod, Aorg, Actx));
        fmpz_mpoly_clear(prod, Actx);
        fmpz_mpoly_clear(Aorg, Actx);
    }
#endif

    FLINT_ASSERT(success == 0 || success == 1);
    return success;
}


/* assume f is a square free factorization, replace it by an irreducible one */
int fmpz_mpoly_factor_irred(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong i, j;
    fmpz_mpolyv_t t;
    fmpz_mpoly_factor_t g;

    fmpz_mpolyv_init(t, ctx);
    fmpz_mpoly_factor_init(g, ctx);

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

    return success;
}


int fmpz_mpoly_factor(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    success = fmpz_mpoly_factor_squarefree(f, A, ctx) &&
              fmpz_mpoly_factor_irred(f, ctx, USE_ZAS | USE_WANG | USE_ZIP);
    FLINT_ASSERT(!success || fmpz_mpoly_factor_matches(A, f, ctx));
    return success;
}


int fmpz_mpoly_factor_zassenhaus(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    success = fmpz_mpoly_factor_squarefree(f, A, ctx) &&
              fmpz_mpoly_factor_irred(f, ctx, USE_ZAS);
    FLINT_ASSERT(!success || fmpz_mpoly_factor_matches(A, f, ctx));
    return success;
}


int fmpz_mpoly_factor_wang(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    success = fmpz_mpoly_factor_squarefree(f, A, ctx) &&
              fmpz_mpoly_factor_irred(f, ctx, USE_WANG);
    FLINT_ASSERT(!success || fmpz_mpoly_factor_matches(A, f, ctx));
    return success;
}


int fmpz_mpoly_factor_zippel(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    success = fmpz_mpoly_factor_squarefree(f, A, ctx) &&
              fmpz_mpoly_factor_irred(f, ctx, USE_ZIP);
    FLINT_ASSERT(!success || fmpz_mpoly_factor_matches(A, f, ctx));
    return success;
}

