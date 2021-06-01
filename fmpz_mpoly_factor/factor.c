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
static int _factor_irred_compressed(
    fmpz_mpolyv_t Af,
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong i;
    slong * Adegs;
    flint_bitcnt_t Abits;
    fmpz_poly_t u;
    fmpz_poly_factor_t uf;
    fmpz_mpoly_t lcA;
    fmpz_mpoly_factor_t lcAf;
    zassenhaus_prune_t Z;
    flint_rand_t state;

    FLINT_ASSERT(A->length > 1 || fmpz_sgn(A->coeffs + 0) > 0);
    FLINT_ASSERT(fmpz_mpoly_degrees_fit_si(A, ctx));
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(fmpz_mpoly_is_canonical(A, ctx));

    if (A->bits > FLINT_BITS &&
        !fmpz_mpoly_repack_bits_inplace(A, FLINT_BITS, ctx))
    {
        return 0;
    }

    Abits = A->bits;

    flint_randinit(state);
    fmpz_poly_init(u);
    fmpz_poly_factor_init(uf);
    fmpz_mpoly_init(lcA, ctx);
    fmpz_mpoly_factor_init(lcAf, ctx);
    zassenhaus_prune_init(Z);
    Adegs = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, slong);

    fmpz_mpoly_degrees_si(Adegs, A, ctx);

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
        int zero_ok, trying_zero, image_count, sqrfree, irr_fac;
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
            zassenhaus_prune_add_factor(Z, fmpz_poly_degree(uf->p + i),
                                           uf->exp[i]);
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

        if (!(algo & (MPOLY_FACTOR_USE_WANG | MPOLY_FACTOR_USE_ZIP)))
            goto try_zassenhaus;

        _fmpz_mpoly_get_lead0(lcA, A, ctx);
        if (!fmpz_mpoly_factor_squarefree(lcAf, lcA, ctx))
            goto try_zassenhaus;

        irr_fac = 1;
        for (i = 0; i < lcAf->num; i++)
            irr_fac = irr_fac && lcAf->poly[i].length < 4;

        irr_fac = irr_fac && fmpz_mpoly_factor_irred(lcAf, ctx, algo);

        if (!(algo & MPOLY_FACTOR_USE_ZIP))
        {
            success = fmpz_mpoly_factor_irred_wang(Af, A,
                                         lcAf, irr_fac, lcA, ctx, state, Z, 1);
        }
        else if (!(algo & MPOLY_FACTOR_USE_WANG))
        {
            success = fmpz_mpoly_factor_irred_zippel(Af, A,
                                            lcAf, irr_fac, lcA, ctx, state, Z);
        }
        else
        {
            if (density > 0.002 && zero_ok)
                success = fmpz_mpoly_factor_irred_wang(Af, A,
                                         lcAf, irr_fac, lcA, ctx, state, Z, 0);

            if (success == 0 && density > 0.04)
                success = fmpz_mpoly_factor_irred_wang(Af, A,
                                         lcAf, irr_fac, lcA, ctx, state, Z, 1);

            if (success == 0)
                success = fmpz_mpoly_factor_irred_zippel(Af, A,
                                            lcAf, irr_fac, lcA, ctx, state, Z);

            if (success == 0)
                success = fmpz_mpoly_factor_irred_wang(Af, A,
                                         lcAf, irr_fac, lcA, ctx, state, Z, 1);
        }

    try_zassenhaus:

        if (algo & MPOLY_FACTOR_USE_ZAS)
        {
            if (success == 0)
                success = fmpz_mpoly_factor_irred_zassenhaus(Af, A, ctx, Z);
        }

        success = (success > 0);
    }

cleanup:

    flint_randclear(state);
    fmpz_poly_clear(u);
    fmpz_poly_factor_clear(uf);
    fmpz_mpoly_clear(lcA, ctx);
    fmpz_mpoly_factor_clear(lcAf, ctx);
    zassenhaus_prune_clear(Z);
    flint_free(Adegs);

    FLINT_ASSERT(success == 0 || success == 1);
    return success;
}


/* replace f by a content factorization, f is already squarefree */
static int _refine_content_factors(
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
                    _fmpz_mpoly_from_univar(c, bits, u, v, ctx);
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
    A is squarefree with positive lead coeff.

    return 1 for success, 0 for failure
*/
static int _factor_irred(
    fmpz_mpolyv_t Af,
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t Actx,
    unsigned int algo)
{
    int success;
    slong i, j;
    flint_bitcnt_t Abits;
    mpoly_compression_t M;
#if FLINT_WANT_ASSERT
    fmpz_mpoly_t Aorg;

    fmpz_mpoly_init(Aorg, Actx);
    fmpz_mpoly_set(Aorg, A, Actx);
#endif

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(fmpz_sgn(A->coeffs + 0) > 0);

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

    if (A->bits > FLINT_BITS &&
        !fmpz_mpoly_repack_bits_inplace(A, FLINT_BITS, Actx))
    {
        success = 0;
        goto cleanup_less;
    }

    Abits = A->bits;

    mpoly_compression_init(M);
    mpoly_compression_set(M, A->exps, A->bits, A->length, Actx->minfo);

    if (M->is_irred)
    {
        fmpz_mpolyv_fit_length(Af, 1, Actx);
        Af->length = 1;
        fmpz_mpoly_swap(Af->coeffs + 0, A, Actx);
        success = 1;
    }
    else if (M->is_trivial)
    {
        success = _factor_irred_compressed(Af, A, Actx, algo);
    }
    else
    {
        fmpz_mpoly_ctx_t Lctx;
        fmpz_mpoly_t L;
        fmpz_mpoly_univar_t U;
        fmpz_mpoly_t t;
        fmpz_mpolyv_t Lf, tf, sf;

        fmpz_mpoly_ctx_init(Lctx, M->mvars, ORD_LEX);
        fmpz_mpoly_init(L, Lctx);
        fmpz_mpolyv_init(Lf, Lctx);
        fmpz_mpoly_init(t, Lctx);
        fmpz_mpoly_univar_init(U, Lctx);
        fmpz_mpolyv_init(tf, Lctx);
        fmpz_mpolyv_init(sf, Lctx);

        fmpz_mpoly_compression_do(L, Lctx, A->coeffs, A->length, M);

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
                success = _factor_irred(tf, sf->coeffs + i, Lctx, algo);
                if (!success)
                    goto cleanup_more;

                fmpz_mpolyv_fit_length(Lf, Lf->length + tf->length, Lctx);
                for (j = 0; j < tf->length; j++)
                    fmpz_mpoly_swap(Lf->coeffs + Lf->length++, tf->coeffs + j, Lctx);
            }
        }
        else
        {
            success = _factor_irred_compressed(Lf, L, Lctx, algo);
        }

    cleanup_more:

        fmpz_mpoly_clear(t, Lctx);
        fmpz_mpoly_univar_clear(U, Lctx);
        fmpz_mpolyv_clear(tf, Lctx);
        fmpz_mpolyv_clear(sf, Lctx);

        if (success)
        {
            fmpz_mpolyv_fit_length(Af, Lf->length, Actx);
            Af->length = Lf->length;
            for (i = 0; i < Lf->length; i++)
            {
                fmpz_mpoly_compression_undo(Af->coeffs + i, Abits, Actx,
                                                      Lf->coeffs + i, Lctx, M);
            }
        }

        fmpz_mpoly_clear(L, Lctx);
        fmpz_mpolyv_clear(Lf, Lctx);
        fmpz_mpoly_ctx_clear(Lctx);
    }

    mpoly_compression_clear(M);

cleanup_less:

#if FLINT_WANT_ASSERT
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


/*
    for each factor A in f, assume A satifies _factor_irred requirements:
        A is primitive w.r.t to any variable appearing in A.
        A is squarefree with positive lead coeff.

    Replace f by and irreducible factorization.
*/
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
        success = _factor_irred(t, f->poly + j, ctx, algo);
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


/*
    append factor(f)^e to g
    assumping f is compressed and content free
*/
static int _compressed_content_to_irred(
    fmpz_mpoly_factor_t g,
    fmpz_mpoly_t f,
    const fmpz_t e,
    const fmpz_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong j, k;
    fmpz_mpoly_factor_t h;
    fmpz_mpolyv_t v;

    fmpz_mpoly_factor_init(h, ctx);
    fmpz_mpolyv_init(v, ctx);

    success = _fmpz_mpoly_factor_squarefree(h, f, e, ctx);
    if (!success)
        goto cleanup;

    for (j = 0; j < h->num; j++)
    {
        success = h->num > 1 ? _factor_irred(v, h->poly + j, ctx, algo) :
                           _factor_irred_compressed(v, h->poly + j, ctx, algo);
        if (!success)
            goto cleanup;

        fmpz_mpoly_factor_fit_length(g, g->num + v->length, ctx);
        for (k = 0; k < v->length; k++)
        {
            fmpz_set(g->exp + g->num, h->exp + j);
            fmpz_mpoly_swap(g->poly + g->num, v->coeffs + k, ctx);
            g->num++;
        }
    }

cleanup:

    fmpz_mpoly_factor_clear(h, ctx);
    fmpz_mpolyv_clear(v, ctx);

    return success;
}


/*
    no assumptions on A,
    returned factors are irreducible
*/
static int _factor(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx,
    unsigned int algo)
{
    int success;
    slong i, j;
    flint_bitcnt_t bits;
    fmpz_mpoly_factor_t g;
    mpoly_compression_t M;

    if (!fmpz_mpoly_factor_content(f, A, ctx))
        return 0;

    fmpz_mpoly_factor_init(g, ctx);
    mpoly_compression_init(M);

    /* write into g */
    fmpz_swap(g->constant, f->constant);
    g->num = 0;
    for (i = 0; i < f->num; i++)
    {
        if (f->poly[i].length < 2)
        {
            fmpz_mpoly_factor_fit_length(g, g->num + 1, ctx);
            fmpz_mpoly_swap(g->poly + g->num, f->poly + i, ctx);
            fmpz_swap(g->exp + g->num, f->exp + i);
            g->num++;
            continue;
        }

        if (f->poly[i].bits > FLINT_BITS &&
            !fmpz_mpoly_repack_bits_inplace(f->poly + i, FLINT_BITS, ctx))
        {
            success = 0;
            goto cleanup;
        }

        bits = f->poly[i].bits;

        mpoly_compression_set(M, f->poly[i].exps, bits, f->poly[i].length, ctx->minfo);
        if (M->is_irred)
        {
            fmpz_mpoly_factor_fit_length(g, g->num + 1, ctx);
            fmpz_mpoly_swap(g->poly + g->num, f->poly + i, ctx);
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
            fmpz_mpoly_ctx_t Lctx;
            fmpz_mpoly_t L;
            fmpz_mpoly_factor_t h;

            /* compression may have messed up the content factorization */
            fmpz_mpoly_ctx_init(Lctx, M->mvars, ORD_LEX);
            fmpz_mpoly_init(L, Lctx);
            fmpz_mpoly_factor_init(h, Lctx);

            fmpz_mpoly_compression_do(L, Lctx, f->poly[i].coeffs,
                                               f->poly[i].length, M);
            if (M->is_perm)
            {
                success = _compressed_content_to_irred(h, L, f->exp + i, Lctx, algo);
                fmpz_one(f->exp + i);
            }
            else
            {
                success = fmpz_mpoly_factor_squarefree(h, L, Lctx) &&
                          fmpz_mpoly_factor_irred(h, Lctx, algo);
            }

            if (success)
            {
                FLINT_ASSERT(fmpz_is_one(h->constant));
                fmpz_mpoly_factor_fit_length(g, g->num + h->num, ctx);
                for (j = 0; j < h->num; j++)
                {
                    fmpz_mul(g->exp + g->num, f->exp + i, h->exp + j);
                    fmpz_mpoly_compression_undo(g->poly + g->num, bits, ctx,
                                                         h->poly + j, Lctx, M);
                    g->num++;
                }
            }

            fmpz_mpoly_factor_clear(h, Lctx);
            fmpz_mpoly_clear(L, Lctx);
            fmpz_mpoly_ctx_clear(Lctx);

            if (!success)
                goto cleanup;
        }
    }

    fmpz_mpoly_factor_swap(f, g, ctx);

    success = 1;

cleanup:

    fmpz_mpoly_factor_clear(g, ctx);
    mpoly_compression_clear(M);

    FLINT_ASSERT(!success || fmpz_mpoly_factor_matches(A, f, ctx));
    return success;
}

int fmpz_mpoly_factor(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    return _factor(f, A, ctx, MPOLY_FACTOR_USE_ALL);
}


int fmpz_mpoly_factor_zassenhaus(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    return _factor(f, A, ctx, MPOLY_FACTOR_USE_ZAS);
}


int fmpz_mpoly_factor_wang(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    return _factor(f, A, ctx, MPOLY_FACTOR_USE_WANG);
}


int fmpz_mpoly_factor_zippel(
    fmpz_mpoly_factor_t f,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    return _factor(f, A, ctx, MPOLY_FACTOR_USE_ZIP);
}

