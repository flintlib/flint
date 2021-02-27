/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly_factor.h"

static void fmpz_mpolyu_norm_degrees(
    fmpz_t norm,
    slong * degrees,
    const fmpz_mpolyu_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    fmpz_t M;
    slong * tdegs;
    TMP_INIT;

    fmpz_init(M);

    TMP_START;

    tdegs = TMP_ARRAY_ALLOC(ctx->minfo->nvars, slong);

    fmpz_zero(norm);

    for (j = 0; j < ctx->minfo->nvars + 1; j++)
        degrees[j] = 0;

    for (i = 0; i < A->length; i++)
    {
        _fmpz_vec_height(M, A->coeffs[i].coeffs, A->coeffs[i].length);
        if (fmpz_cmp(norm, M) < 0)
            fmpz_swap(norm, M);

        fmpz_mpoly_degrees_si(tdegs, A->coeffs + i, ctx);

        degrees[0] = FLINT_MAX(degrees[0], A->exps[i]);
        for (j = 0; j < ctx->minfo->nvars; j++)
            degrees[1 + j] = FLINT_MAX(degrees[1 + j], tdegs[j]);
    }

    fmpz_clear(M);

    TMP_END;
}

static flint_bitcnt_t fmpz_mpolyu_gcd_bitbound(
    const fmpz_mpolyu_t A,
    const fmpz_mpolyu_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bound = UWORD_MAX;
    fmpz_t Anorm, Bnorm, M;
    slong * Adegs, * Bdegs;
    TMP_INIT;

    TMP_START;

    fmpz_init(Anorm);
    fmpz_init(Bnorm);
    fmpz_init(M);

    Adegs = TMP_ARRAY_ALLOC(ctx->minfo->nvars + 1, slong);
    Bdegs = TMP_ARRAY_ALLOC(ctx->minfo->nvars + 1, slong);

    fmpz_mpolyu_norm_degrees(Anorm, Adegs, A, ctx);
    fmpz_mpolyu_norm_degrees(Bnorm, Bdegs, B, ctx);

    if (fmpz_mpoly_factor_bound_si(M, Anorm, Adegs, ctx->minfo->nvars + 1))
        bound = FLINT_MIN(bound, fmpz_bits(M));

    if (fmpz_mpoly_factor_bound_si(M, Bnorm, Bdegs, ctx->minfo->nvars + 1))
        bound = FLINT_MIN(bound, fmpz_bits(M));

    fmpz_clear(Anorm);
    fmpz_clear(Bnorm);
    fmpz_clear(M);

    TMP_END;
    return bound;
}

int fmpz_mpolyu_gcdm_zippel(
    fmpz_mpolyu_t G,
    fmpz_mpolyu_t Abar,
    fmpz_mpolyu_t Bbar,
    fmpz_mpolyu_t A,
    fmpz_mpolyu_t B,
    const fmpz_mpoly_ctx_t ctx,
    flint_rand_t randstate)
{
    flint_bitcnt_t coeffbitbound;
    flint_bitcnt_t coeffbits;
    slong degbound;
    int success, changed;
    mp_limb_t p, t, gammap;
    fmpz_t gamma, pp, gammapp, modulus;
    nmod_mpolyu_t Ap, Bp, Gp, Abarp, Bbarp, Gform;
    fmpz_mpolyu_t H;
    nmod_mpoly_ctx_t ctxp;

    fmpz_init(pp);
    fmpz_init(gammapp);
    fmpz_init_set_si(modulus, 1);
    fmpz_init(gamma);
    fmpz_gcd(gamma, fmpz_mpolyu_leadcoeff(A), fmpz_mpolyu_leadcoeff(B));

    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(A->bits == G->bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(A->exps[A->length - 1] == 0);
    FLINT_ASSERT(B->exps[B->length - 1] == 0);

    degbound = FLINT_MIN(A->exps[0], B->exps[0]);

    coeffbitbound = fmpz_mpolyu_gcd_bitbound(A, B, ctx);
    if (n_add_checked(&coeffbitbound, coeffbitbound, fmpz_bits(gamma)))
        coeffbitbound = UWORD_MAX;

    nmod_mpoly_ctx_init(ctxp, ctx->minfo->nvars, ORD_LEX, 2);

    nmod_mpolyu_init(Ap, A->bits, ctxp);
    nmod_mpolyu_init(Bp, A->bits, ctxp);
    nmod_mpolyu_init(Gp, A->bits, ctxp);
    nmod_mpolyu_init(Abarp, A->bits, ctxp);
    nmod_mpolyu_init(Bbarp, A->bits, ctxp);
    nmod_mpolyu_init(Gform, A->bits, ctxp);

    fmpz_mpolyu_init(H, A->bits, ctx);

    p = UWORD(1) << (FLINT_BITS - 1);

choose_prime_outer:

    if (p >= UWORD_MAX_PRIME)
    {
        /* ran out of machine primes - absolute failure */
        success = 0;
        goto cleanup;
    }
    p = n_nextprime(p, 1);

    /* make sure mod p reduction does not kill both lc(A) and lc(B) */
    fmpz_set_ui(pp, p);
    fmpz_mod(gammapp, gamma, pp);
    gammap = fmpz_get_ui(gammapp);
    if (gammap == UWORD(0))
        goto choose_prime_outer;

    nmod_mpoly_ctx_change_modulus(ctxp, p);

    /* make sure mod p reduction does not kill either A or B */
    fmpz_mpolyu_interp_reduce_p(Ap, ctxp, A, ctx);
    fmpz_mpolyu_interp_reduce_p(Bp, ctxp, B, ctx);
    if (Ap->length == 0 || Bp->length == 0)
        goto choose_prime_outer;

    success = nmod_mpolyu_gcdp_zippel(Gp, Abarp, Bbarp, Ap, Bp,
                                       ctx->minfo->nvars - 1, ctxp, randstate);
    if (!success || Gp->exps[0] > degbound)
        goto choose_prime_outer;
    degbound = Gp->exps[0];

    t = nmod_mpolyu_leadcoeff(Gp, ctxp);
    t = nmod_inv(t, ctxp->mod);
    t = nmod_mul(t, gammap, ctxp->mod);
    nmod_mpolyu_scalar_mul_nmod(Gp, t, ctxp);

    if (Gp->length == 1 && Gp->exps[0] == 0)
    {
        FLINT_ASSERT(nmod_mpoly_is_ui(Gp->coeffs + 0, ctxp));
        FLINT_ASSERT(!nmod_mpoly_is_zero(Gp->coeffs + 0, ctxp));
        fmpz_mpolyu_one(G, ctx);
        fmpz_mpolyu_swap(Abar, A, ctx);
        fmpz_mpolyu_swap(Bbar, B, ctx);
        success = 1;
        goto cleanup;
    }

    nmod_mpolyu_setform(Gform, Gp, ctxp);
    fmpz_mpolyu_interp_lift_p(H, ctx, Gp, ctxp);
    fmpz_set_ui(modulus, p);

choose_prime_inner:

    if (p >= UWORD_MAX_PRIME)
    {
        /* ran out of machine primes - absolute failure */
        success = 0;
        goto cleanup;
    }
    p = n_nextprime(p, 1);

    /* make sure mod p reduction does not kill both lc(A) and lc(B) */
    fmpz_set_ui(pp, p);
    fmpz_mod(gammapp, gamma, pp);
    gammap = fmpz_get_ui(gammapp);
    if (gammap == UWORD(0))
        goto choose_prime_inner;

    nmod_mpoly_ctx_change_modulus(ctxp, p);

    /* make sure mod p reduction does not kill either A or B */
    fmpz_mpolyu_interp_reduce_p(Ap, ctxp, A, ctx);
    fmpz_mpolyu_interp_reduce_p(Bp, ctxp, B, ctx);
    if (Ap->length == 0 || Bp->length == 0)
        goto choose_prime_inner;

    switch (nmod_mpolyu_gcds_zippel(Gp, Ap, Bp, Gform, ctx->minfo->nvars,
                                                   ctxp, randstate, &degbound))
    {
        default:
            FLINT_ASSERT(0);
        case nmod_gcds_form_main_degree_too_high:
        case nmod_gcds_form_wrong:
        case nmod_gcds_no_solution:
            goto choose_prime_outer;
        case nmod_gcds_scales_not_found:
        case nmod_gcds_eval_point_not_found:
        case nmod_gcds_eval_gcd_deg_too_high:
            goto choose_prime_inner;
        case nmod_gcds_success:
            break;
    }

    if (nmod_mpolyu_leadcoeff(Gp, ctxp) == UWORD(0))
        goto choose_prime_inner;

    t = nmod_mpolyu_leadcoeff(Gp, ctxp);
    t = nmod_inv(t, ctxp->mod);
    t = nmod_mul(t, gammap, ctxp->mod);
    nmod_mpolyu_scalar_mul_nmod(Gp, t, ctxp);

    changed = fmpz_mpolyu_interp_mcrt_p(&coeffbits, H, ctx, modulus, Gp, ctxp);
    fmpz_mul_ui(modulus, modulus, p);

    if (changed)
    {
        if (coeffbits > coeffbitbound)
        {
            goto choose_prime_outer;
        }
        goto choose_prime_inner;
    }

    fmpz_mpolyu_fmpz_content(pp, H, ctx);
    fmpz_mpolyu_divexact_fmpz(G, H, pp, ctx);

    if (   !fmpz_mpolyuu_divides(Abar, A, G, 1, ctx)
        || !fmpz_mpolyuu_divides(Bbar, B, G, 1, ctx))
    {
        goto choose_prime_inner;
    }

    success = 1;

cleanup:

    nmod_mpolyu_clear(Ap, ctxp);
    nmod_mpolyu_clear(Bp, ctxp);
    nmod_mpolyu_clear(Gp, ctxp);
    nmod_mpolyu_clear(Abarp, ctxp);
    nmod_mpolyu_clear(Bbarp, ctxp);
    nmod_mpolyu_clear(Gform, ctxp);

    fmpz_mpolyu_clear(H, ctx);

    nmod_mpoly_ctx_clear(ctxp);

    fmpz_clear(gammapp);
    fmpz_clear(gamma);
    fmpz_clear(modulus);
    fmpz_clear(pp);

    return success;
}

