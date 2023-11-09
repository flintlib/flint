/*
    Copyright (C) 2018, 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"
#include "fmpz_mpoly_factor.h"

/* return an n with |gcd(A,B)|_infty < 2^n or return UWORD_MAX */
static flint_bitcnt_t fmpz_mpoly_gcd_bitbound(
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bound = UWORD_MAX;
    fmpz_t norm, M;
    slong * degs;
    TMP_INIT;

    TMP_START;

    fmpz_init(norm);
    fmpz_init(M);

    degs = TMP_ARRAY_ALLOC(ctx->minfo->nvars, slong);

    fmpz_mpoly_degrees_si(degs, A, ctx);
    _fmpz_vec_height(norm, A->coeffs, A->length);
    if (fmpz_mpoly_factor_bound_si(M, norm, degs, ctx->minfo->nvars))
        bound = FLINT_MIN(bound, fmpz_bits(M));

    fmpz_mpoly_degrees_si(degs, B, ctx);
    _fmpz_vec_height(norm, B->coeffs, B->length);
    if (fmpz_mpoly_factor_bound_si(M, norm, degs, ctx->minfo->nvars))
        bound = FLINT_MIN(bound, fmpz_bits(M));

    fmpz_clear(norm);
    fmpz_clear(M);

    TMP_END;

    return bound;
}

int fmpz_mpolyl_gcd_zippel(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    flint_bitcnt_t coeffbitbound;
    flint_bitcnt_t coeffbits;
    flint_bitcnt_t bits = G->bits;
    int success, changed;
    slong i, j, Gdegbound, Gdeg, req_zip_images;
    mp_limb_t p, t, gammap;
    fmpz_t c, gamma, modulus;
    nmod_mpoly_t Ap, Bp, Gp, Abarp, Bbarp;
    nmod_mpoly_ctx_t ctxp;
    n_poly_t Amarks, Bmarks, Gmarks;
    slong * perm = NULL;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(ctx->minfo->nvars > 1);
    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == Abar->bits);
    FLINT_ASSERT(bits == Bbar->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == G->bits);

    fmpz_init(c);
    fmpz_init(gamma);
    fmpz_init(modulus);

    nmod_mpoly_ctx_init(ctxp, ctx->minfo->nvars, ORD_LEX, 2);

    nmod_mpoly_init3(Ap, 0, bits, ctxp);
    nmod_mpoly_init3(Bp, 0, bits, ctxp);
    nmod_mpoly_init3(Gp, 0, bits, ctxp);
    nmod_mpoly_init3(Abarp, 0, bits, ctxp);
    nmod_mpoly_init3(Bbarp, 0, bits, ctxp);

    n_poly_init(Amarks);
    n_poly_init(Bmarks);
    n_poly_init(Gmarks);

    fmpz_gcd(gamma, fmpz_mpoly_leadcoeff(A), fmpz_mpoly_leadcoeff(B));

    Gdegbound = fmpz_mpoly_degree_si(A, 0, ctx);
    Gdeg = fmpz_mpoly_degree_si(B, 0, ctx);
    Gdegbound = FLINT_MIN(Gdegbound, Gdeg);

    coeffbitbound = fmpz_mpoly_gcd_bitbound(A, B, ctx);
    if (n_add_checked(&coeffbitbound, coeffbitbound, fmpz_bits(gamma)))
        coeffbitbound = UWORD_MAX;

    p = UWORD(1) << (FLINT_BITS - 1);

outer_loop:

    if (p >= UWORD_MAX_PRIME)
    {
        /* ran out of primes: absolute failure */
        success = 0;
        goto cleanup;
    }
    p = n_nextprime(p, 1);

    nmod_mpoly_ctx_change_modulus(ctxp, p);

    /* make sure mod p reduction does not kill both lc(A) and lc(B) */
    gammap = fmpz_get_nmod(gamma, ctxp->mod);
    if (gammap == 0)
        goto outer_loop;

    /* make sure mod p reduction does not kill either A or B */
    fmpz_mpoly_interp_reduce_p(Ap, ctxp, A, ctx);
    fmpz_mpoly_interp_reduce_p(Bp, ctxp, B, ctx);
    if (Ap->length == 0 || Bp->length == 0)
        goto outer_loop;

    success = nmod_mpolyl_gcdp_zippel_smprime(Gp, Abarp, Bbarp, Ap, Bp,
                                           ctx->minfo->nvars - 1, ctxp, state);
    if (!success)
        goto outer_loop;

    Gdeg = nmod_mpoly_degree_si(Gp, 0, ctxp);
    if (Gdeg > Gdegbound)
        goto outer_loop;

    if (Gdeg == 0 && nmod_mpoly_is_one(Gp, ctxp))
    {
        fmpz_mpoly_one(G, ctx);
        fmpz_mpoly_set(Abar, A, ctx);
        fmpz_mpoly_set(Bbar, B, ctx);
        success = 1;
        goto cleanup;
    }

    Gdegbound = Gdeg;

    t = nmod_div(gammap, nmod_mpoly_leadcoeff(Gp, ctxp), ctxp->mod);
    nmod_mpoly_scalar_mul_nmod_invertible(Gp, Gp, t, ctxp);

    fmpz_mpoly_interp_lift_p(G, ctx, Gp, ctxp);
    fmpz_set_ui(modulus, p);

    /* TODO: a divisibility test here if coeffs are small */

    mpoly1_fill_marks(&Gmarks->coeffs, &Gmarks->length, &Gmarks->alloc,
                                         G->exps, G->length, bits, ctx->minfo);

    perm = FLINT_ARRAY_REALLOC(perm, Gmarks->length, slong);
    for (i = 0; i < Gmarks->length; i++)
        perm[i] = i;

#define length(k) Gmarks->coeffs[(k)+1] - Gmarks->coeffs[k]

    for (i = 1; i < Gmarks->length; i++)
        for (j = i; j > 0 && length(perm[j-1]) > length(perm[j]); j--)
            FLINT_SWAP(slong, perm[j-1], perm[j]);

    req_zip_images = Gmarks->length - 2;
    j = 0;
    for (i = 0; i < Gmarks->length; i++)
    {
        req_zip_images += length(i);
        j = FLINT_MAX(j, length(i));
    }

    if (Gmarks->length > 1)
        req_zip_images = req_zip_images / (Gmarks->length - 1);

    req_zip_images = FLINT_MAX(req_zip_images, j);
    req_zip_images += 1;

inner_loop:

    if (p >= UWORD_MAX_PRIME)
    {
        /* ran out of primes: absolute failure */
        success = 0;
        goto cleanup;
    }
    p = n_nextprime(p, 1);

    nmod_mpoly_ctx_change_modulus(ctxp, p);

    /* make sure mod p reduction does not kill both lc(A) and lc(B) */
    gammap = fmpz_get_nmod(gamma, ctxp->mod);
    if (gammap == 0)
        goto inner_loop;

    /* make sure mod p reduction does not kill either A or B */
    fmpz_mpoly_interp_reduce_p(Ap, ctxp, A, ctx);
    fmpz_mpoly_interp_reduce_p(Bp, ctxp, B, ctx);
    if (Ap->length == 0 || Bp->length == 0)
        goto inner_loop;

    success = nmod_mpolyl_gcds_zippel(Gp, Gmarks->coeffs, Gmarks->length,
                       Ap, Bp, perm, req_zip_images, ctx->minfo->nvars, ctxp,
                                            state, &Gdegbound, Amarks, Bmarks);
    if (success == 0)
        goto outer_loop; /* resets modulus */

    if (success < 0 || nmod_mpoly_leadcoeff(Gp, ctxp) == 0)
        goto inner_loop;

    t = nmod_div(gammap, nmod_mpoly_leadcoeff(Gp, ctxp), ctxp->mod);
    nmod_mpoly_scalar_mul_nmod_invertible(Gp, Gp, t, ctxp);

    changed = fmpz_mpoly_interp_mcrt_p(&coeffbits, G, ctx, modulus, Gp, ctxp);
    fmpz_mul_ui(modulus, modulus, p);

    if (changed)
    {
        if (coeffbits > coeffbitbound)
            goto outer_loop; /* resets modulus */

        goto inner_loop;
    }

    _fmpz_vec_content(c, G->coeffs, G->length);
    _fmpz_vec_scalar_divexact_fmpz(G->coeffs, G->coeffs, G->length, c);

    success = fmpz_mpoly_divides(Abar, A, G, ctx) &&
              fmpz_mpoly_divides(Bbar, B, G, ctx);

    if (success)
        goto cleanup;

    /* restore interpolated state */
    _fmpz_vec_scalar_mul_fmpz(G->coeffs, G->coeffs, G->length, c);

    goto inner_loop;

cleanup:

    flint_free(perm);

    n_poly_clear(Amarks);
    n_poly_clear(Bmarks);
    n_poly_clear(Gmarks);

    nmod_mpoly_clear(Ap, ctxp);
    nmod_mpoly_clear(Bp, ctxp);
    nmod_mpoly_clear(Gp, ctxp);
    nmod_mpoly_clear(Abarp, ctxp);
    nmod_mpoly_clear(Bbarp, ctxp);

    nmod_mpoly_ctx_clear(ctxp);

    fmpz_clear(c);
    fmpz_clear(gamma);
    fmpz_clear(modulus);

    return success;
}

