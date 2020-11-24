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

/*
    Find a bound on the bits of the coefficients of gcd(A,B).
    If this overflows a flint_bitcnt_t, the max flint_bitcnt_t is returned.
    We will apply a Kronecker substitution and use the Landau-Mignotte bound
        for univariates.

         min(deg A, deg B)
        2                  * gcd(A[0], B[0])

                 * min( L2Norm(A) / A[0] , L2Norm(B) / B[0] )

    This number is almost certainly not going to be used,
        so it just needs to be correct and not nececessary as tight as possible.
*/
flint_bitcnt_t fmpz_mpolyu_gcd_bitbound(const fmpz_t gcdlc,
                            const fmpz_mpolyu_t A, const fmpz_mpolyu_t B,
                       const fmpz_mpoly_ctx_t ctx, const mpoly_zipinfo_t zinfo)
{
    slong i;
    fmpz_t n, an, bn;
    ulong max, len;
    flint_bitcnt_t r;

    /* find the degree of the kronecker substitution in the lesser variables */
    fmpz_init_set_ui(n, UWORD(1));
    for (i = 0; i < zinfo->nvars - 1; i++)
    {
        fmpz_mul_ui(n, n, FLINT_MAX((ulong)(zinfo->Adegs[i]),
                                    (ulong)(zinfo->Bdegs[i])) + UWORD(1));
    }

    /* n = min(deg A, deg B) after the KS */
    i = zinfo->nvars - 1;
    fmpz_addmul_ui(n, n, FLINT_MIN((ulong)(zinfo->Adegs[i]),
                                   (ulong)(zinfo->Bdegs[i])));

    /* n += log2(gcd(A[0], B[0])) */
    fmpz_add_ui(n, n, fmpz_bits(gcdlc));

    len = max = UWORD(0);
    for (i = 0; i < A->length; i++)
    {
        len += (A->coeffs + i)->length;
        max = FLINT_MAX(max, FLINT_ABS(
                    _fmpz_vec_max_bits((A->coeffs + i)->coeffs,
                                       (A->coeffs + i)->length)
              ));
    }
    fmpz_init_set_ui(an, n_clog(len, UWORD(2))/UWORD(2));
    fmpz_add_ui(an, an, max);
    fmpz_sub_ui(an, an, fmpz_bits(fmpz_mpolyu_leadcoeff(A)));
    FLINT_ASSERT(fmpz_sgn(an) >= 0);

    len = max = UWORD(0);
    for (i = 0; i < B->length; i++)
    {
        len += (B->coeffs + i)->length;
        max = FLINT_MAX(max, FLINT_ABS(
                    _fmpz_vec_max_bits((B->coeffs + i)->coeffs,
                                       (B->coeffs + i)->length)
              ));
    }
    fmpz_init_set_ui(bn, n_clog(len, UWORD(2))/UWORD(2));
    fmpz_add_ui(bn, bn, max);
    fmpz_sub_ui(bn, bn, fmpz_bits(fmpz_mpolyu_leadcoeff(B)));
    FLINT_ASSERT(fmpz_sgn(bn) >= 0);

    /* n += log2( min( L2Norm(A) / A[0] , L2Norm(B) / B[0] ) ) */
    fmpz_add(n, n, fmpz_cmp(an, bn) < 0 ? an : bn);

    FLINT_ASSERT(fmpz_sgn(n) > 0);

    r = fmpz_abs_fits_ui(n) ? fmpz_get_ui(n) : -UWORD(1);

    fmpz_clear(n);
    fmpz_clear(an);
    fmpz_clear(bn);

    return r;
}

int fmpz_mpolyu_gcdm_zippel(
    fmpz_mpolyu_t G,
    fmpz_mpolyu_t Abar,
    fmpz_mpolyu_t Bbar,
    fmpz_mpolyu_t A,
    fmpz_mpolyu_t B,
    const fmpz_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
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

    coeffbitbound = fmpz_mpolyu_gcd_bitbound(gamma, A, B, ctx, zinfo);

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
                                ctx->minfo->nvars - 1, ctxp, zinfo, randstate);
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


int fmpz_mpoly_gcd_zippel(
    fmpz_mpoly_t G,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    flint_bitcnt_t wbits;
    flint_rand_t randstate;
    int success = 0;
    mpoly_zipinfo_t zinfo;
    fmpz_mpoly_ctx_t uctx;
    fmpz_mpolyu_t Au, Bu, Gu, Abaru, Bbaru;
    fmpz_mpoly_t Ac, Bc, Gc;
    ulong * shift, * stride;
    ulong amin, bmin;

    if (fmpz_mpoly_is_zero(A, ctx))
    {
        if (fmpz_mpoly_is_zero(B, ctx))
        {
            fmpz_mpoly_zero(G, ctx);
            return 1;
        }

        if (fmpz_sgn(B->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, B, ctx);
            return 1;
        }
        else
        {
            fmpz_mpoly_set(G, B, ctx);
            return 1;
        }
    }

    if (fmpz_mpoly_is_zero(B, ctx))
    {
        if (fmpz_sgn(A->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, A, ctx);
            return 1;
        }
        else
        {
            fmpz_mpoly_set(G, A, ctx);
            return 1;
        }
    }

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
    {
        return 0;
    }

    shift = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    stride = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        shift[i] = 0;
        stride[i] = 1;
    }

    if (ctx->minfo->nvars == 1)
    {
        fmpz_poly_t a, b, g;
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(g);
        _fmpz_mpoly_to_fmpz_poly_deflate(a, A, 0, shift, stride, ctx);
        _fmpz_mpoly_to_fmpz_poly_deflate(b, B, 0, shift, stride, ctx);
        fmpz_poly_gcd(g, a, b);
        _fmpz_mpoly_from_fmpz_poly_inflate(G, A->bits, g, 0, shift, stride, ctx);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(g);

        success = 1;
        goto cleanup1;
    }

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(!fmpz_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!fmpz_mpoly_is_zero(B, ctx));
    FLINT_ASSERT(ctx->minfo->nvars > WORD(1));

    flint_randinit(randstate);

    mpoly_zipinfo_init(zinfo, ctx->minfo->nvars);
    fmpz_mpoly_degrees_si(zinfo->Adegs, A, ctx);
    fmpz_mpoly_degrees_si(zinfo->Bdegs, B, ctx);
    for (i = 0; i < ctx->minfo->nvars; i++)
        zinfo->perm[i] = i;

    wbits = FLINT_MAX(A->bits, B->bits);

    fmpz_mpoly_ctx_init(uctx, ctx->minfo->nvars - 1, ORD_LEX);

    fmpz_mpolyu_init(Au, wbits, uctx);
    fmpz_mpolyu_init(Bu, wbits, uctx);
    fmpz_mpolyu_init(Gu, wbits, uctx);
    fmpz_mpolyu_init(Abaru, wbits, uctx);
    fmpz_mpolyu_init(Bbaru, wbits, uctx);
    fmpz_mpoly_init3(Ac, 0, wbits, uctx);
    fmpz_mpoly_init3(Bc, 0, wbits, uctx);
    fmpz_mpoly_init3(Gc, 0, wbits, uctx);

    fmpz_mpoly_to_mpolyu_perm_deflate_threaded_pool(Au, uctx, A, ctx, zinfo->perm,
                                                 shift, stride, NULL, NULL, 0);
    fmpz_mpoly_to_mpolyu_perm_deflate_threaded_pool(Bu, uctx, B, ctx, zinfo->perm,
                                                 shift, stride, NULL, NULL, 0);

    FLINT_ASSERT(Au->length > 0);
    FLINT_ASSERT(Bu->length > 0);
    amin = Au->exps[Au->length - 1];
    bmin = Bu->exps[Bu->length - 1];
    fmpz_mpolyu_shift_right(Au, amin);
    fmpz_mpolyu_shift_right(Bu, bmin);

    success = fmpz_mpolyu_content_mpoly_threaded_pool(Ac, Au, uctx, NULL, 0);
    success = success &&
              fmpz_mpolyu_content_mpoly_threaded_pool(Bc, Bu, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    fmpz_mpolyu_divexact_mpoly_inplace(Au, Ac, uctx);
    fmpz_mpolyu_divexact_mpoly_inplace(Bu, Bc, uctx);

    /* after removing content, degree bounds in zinfo are still valid bounds */
    success = fmpz_mpolyu_gcdm_zippel(Gu, Abaru, Bbaru, Au, Bu,
                                                       uctx, zinfo, randstate);
    if (!success)
        goto cleanup;

    success = _fmpz_mpoly_gcd_threaded_pool(Gc, wbits, Ac, Bc, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    fmpz_mpolyu_mul_mpoly_inplace(Gu, Gc, uctx);
    fmpz_mpolyu_shift_left(Gu, FLINT_MIN(amin, bmin));

    fmpz_mpoly_from_mpolyu_perm_inflate(G, FLINT_MIN(A->bits, B->bits), ctx,
                                         Gu, uctx, zinfo->perm, shift, stride);
    success = 1;

    if (fmpz_sgn(G->coeffs + 0) < 0)
        fmpz_mpoly_neg(G, G, ctx);

cleanup:

    fmpz_mpolyu_clear(Au, uctx);
    fmpz_mpolyu_clear(Bu, uctx);
    fmpz_mpolyu_clear(Gu, uctx);
    fmpz_mpolyu_clear(Abaru, uctx);
    fmpz_mpolyu_clear(Bbaru, uctx);
    fmpz_mpoly_clear(Ac, uctx);
    fmpz_mpoly_clear(Bc, uctx);
    fmpz_mpoly_clear(Gc, uctx);

    fmpz_mpoly_ctx_clear(uctx);

    mpoly_zipinfo_clear(zinfo);

    flint_randclear(randstate);

cleanup1:

    flint_free(shift);
    flint_free(stride);

    return success;
}
