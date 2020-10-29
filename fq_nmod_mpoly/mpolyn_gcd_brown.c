/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


int fq_nmod_mpolyn_gcd_brown_smprime_bivar(
    fq_nmod_mpolyn_t G,
    fq_nmod_mpolyn_t Abar,
    fq_nmod_mpolyn_t Bbar,
    fq_nmod_mpolyn_t A,
    fq_nmod_mpolyn_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    fq_nmod_t alpha, temp, gammaeval;
    fq_nmod_poly_t Aeval, Beval, Geval, Abareval, Bbareval;
    fq_nmod_mpolyn_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    fq_nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma, trem, leadG;
    fq_nmod_poly_t modulus, tempmod;
    flint_bitcnt_t bits = A->bits;
    slong N, off, shift;
#if FLINT_WANT_ASSERT
    fq_nmod_poly_t leadA, leadB;
#endif

#if FLINT_WANT_ASSERT
    fq_nmod_poly_init(leadA, ctx->fqctx);
    fq_nmod_poly_init(leadB, ctx->fqctx);
    n_fq_poly_get_fq_nmod_poly(leadA, fq_nmod_mpolyn_leadcoeff_poly(A, ctx), ctx->fqctx);
    n_fq_poly_get_fq_nmod_poly(leadB, fq_nmod_mpolyn_leadcoeff_poly(B, ctx), ctx->fqctx);
#endif

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, A->bits, ctx->minfo);

    fq_nmod_poly_init(cA, ctx->fqctx);
    fq_nmod_poly_init(cB, ctx->fqctx);
    fq_nmod_poly_init(cG, ctx->fqctx);
    fq_nmod_poly_init(modulus, ctx->fqctx);
    fq_nmod_poly_init(tempmod, ctx->fqctx);
    fq_nmod_poly_init(trem, ctx->fqctx);
    fq_nmod_poly_init(leadG, ctx->fqctx);
    fq_nmod_poly_init(cAbar, ctx->fqctx);
    fq_nmod_poly_init(cBbar, ctx->fqctx);
    fq_nmod_poly_init(Aeval, ctx->fqctx);
    fq_nmod_poly_init(Beval, ctx->fqctx);
    fq_nmod_poly_init(Geval, ctx->fqctx);
    fq_nmod_poly_init(Abareval, ctx->fqctx);
    fq_nmod_poly_init(Bbareval, ctx->fqctx);

    fq_nmod_mpolyn_init(T, bits, ctx);

    fq_nmod_mpolyn_content_poly(cA, A, ctx);
    fq_nmod_mpolyn_content_poly(cB, B, ctx);
    fq_nmod_mpolyn_divexact_poly(A, A, cA, ctx, tempmod, trem);
    fq_nmod_mpolyn_divexact_poly(B, B, cB, ctx, tempmod, trem);

    fq_nmod_poly_gcd(cG, cA, cB, ctx->fqctx);

    fq_nmod_poly_divrem(cAbar, trem, cA, cG, ctx->fqctx);
    FLINT_ASSERT(fq_nmod_poly_is_zero(trem, ctx->fqctx));
    fq_nmod_poly_divrem(cBbar, trem, cB, cG, ctx->fqctx);
    FLINT_ASSERT(fq_nmod_poly_is_zero(trem, ctx->fqctx));

    fq_nmod_poly_init(gamma, ctx->fqctx);
    n_fq_poly_get_fq_nmod_poly(trem, fq_nmod_mpolyn_leadcoeff_poly(A, ctx), ctx->fqctx);
    n_fq_poly_get_fq_nmod_poly(tempmod, fq_nmod_mpolyn_leadcoeff_poly(B, ctx), ctx->fqctx);
    fq_nmod_poly_gcd(gamma, trem, tempmod, ctx->fqctx);

    ldegA = fq_nmod_mpolyn_lastdeg(A, ctx);
    ldegB = fq_nmod_mpolyn_lastdeg(B, ctx);
    deggamma = fq_nmod_poly_degree(gamma, ctx->fqctx);

    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    fq_nmod_poly_one(modulus, ctx->fqctx);
    fq_nmod_poly_gen(tempmod, ctx->fqctx);
    fq_nmod_poly_neg(tempmod, tempmod, ctx->fqctx);

    fq_nmod_init(gammaeval, ctx->fqctx);
    fq_nmod_init(alpha, ctx->fqctx);
    fq_nmod_init(temp, ctx->fqctx);

    fq_nmod_set_ui(alpha, 0, ctx->fqctx);

choose_prime:   /* prime is v - alpha */

    if (fq_nmod_next(alpha, ctx->fqctx) == 0)
    {
        success = 0;
        goto cleanup;
    }

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    fq_nmod_poly_evaluate_fq_nmod(gammaeval, gamma, alpha, ctx->fqctx);
    if (fq_nmod_is_zero(gammaeval, ctx->fqctx))
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    fq_nmod_mpolyn_interp_reduce_sm_poly(Aeval, A, alpha, ctx);
    fq_nmod_mpolyn_interp_reduce_sm_poly(Beval, B, alpha, ctx);
    FLINT_ASSERT(Aeval->length > 0);
    FLINT_ASSERT(Beval->length > 0);

    fq_nmod_poly_gcd(Geval, Aeval, Beval, ctx->fqctx);
    fq_nmod_poly_divrem(Abareval, trem, Aeval, Geval, ctx->fqctx); FLINT_ASSERT(fq_nmod_poly_is_zero(trem, ctx->fqctx));
    fq_nmod_poly_divrem(Bbareval, trem, Beval, Geval, ctx->fqctx); FLINT_ASSERT(fq_nmod_poly_is_zero(trem, ctx->fqctx));

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(Abareval->length);
    FLINT_ASSERT(Bbareval->length > 0);

    if (fq_nmod_poly_degree(Geval, ctx->fqctx) == 0)
    {
        fq_nmod_mpolyn_one(G, ctx);
        fq_nmod_mpolyn_swap(Abar, A);
        fq_nmod_mpolyn_swap(Bbar, B);
        goto successful_put_content;
    }

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (fq_nmod_poly_degree(Geval, ctx->fqctx) > ((G->exps + N*0)[off]>>shift))
        {
            goto choose_prime;
        }
        else if (fq_nmod_poly_degree(Geval, ctx->fqctx) < ((G->exps + N*0)[off]>>shift))
        {
            fq_nmod_poly_one(modulus, ctx->fqctx);
        }
    }

    fq_nmod_poly_scalar_mul_fq_nmod(Geval, Geval, gammaeval, ctx->fqctx);

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
    {
        fq_nmod_poly_evaluate_fq_nmod(temp, modulus, alpha, ctx->fqctx);
        fq_nmod_inv(temp, temp, ctx->fqctx);
        fq_nmod_poly_scalar_mul_fq_nmod(modulus, modulus, temp, ctx->fqctx);
        fq_nmod_mpolyn_interp_crt_sm_poly(&ldegG, G, T, Geval, modulus, alpha, ctx);
        fq_nmod_mpolyn_interp_crt_sm_poly(&ldegAbar, Abar, T, Abareval, modulus, alpha, ctx);
        fq_nmod_mpolyn_interp_crt_sm_poly(&ldegBbar, Bbar, T, Bbareval, modulus, alpha, ctx);
    }
    else
    {
        fq_nmod_mpolyn_interp_lift_sm_poly(G, Geval, ctx);
        fq_nmod_mpolyn_interp_lift_sm_poly(Abar, Abareval, ctx);
        fq_nmod_mpolyn_interp_lift_sm_poly(Bbar, Bbareval, ctx);
        ldegG = 0;
        ldegAbar = 0;
        ldegBbar = 0;
    }
    fq_nmod_poly_set_coeff(tempmod, 0, alpha, ctx->fqctx);
    fq_nmod_poly_mul(modulus, modulus, tempmod, ctx->fqctx);

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) < bound)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (   deggamma + ldegA == ldegG + ldegAbar
        && deggamma + ldegB == ldegG + ldegBbar )
    {
        goto successful;
    }

    fq_nmod_poly_one(modulus, ctx->fqctx);
    goto choose_prime;

successful:

    fq_nmod_mpolyn_content_poly(modulus, G, ctx);
    fq_nmod_mpolyn_divexact_poly(G, G, modulus, ctx, tempmod, trem);
    n_fq_poly_get_fq_nmod_poly(leadG, fq_nmod_mpolyn_leadcoeff_poly(G, ctx), ctx->fqctx);
    fq_nmod_mpolyn_divexact_poly(Abar, Abar, leadG, ctx, tempmod, trem);
    fq_nmod_mpolyn_divexact_poly(Bbar, Bbar, leadG, ctx, tempmod, trem);

successful_put_content:

    fq_nmod_mpolyn_mul_poly(G, G, cG, ctx, tempmod);
    fq_nmod_mpolyn_mul_poly(Abar, Abar, cAbar, ctx, tempmod);
    fq_nmod_mpolyn_mul_poly(Bbar, Bbar, cBbar, ctx, tempmod);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success)
    {
        n_fq_poly_get_fq_nmod_poly(leadG, fq_nmod_mpolyn_leadcoeff_poly(G, ctx), ctx->fqctx);
        n_fq_poly_get_fq_nmod_poly(trem, fq_nmod_mpolyn_leadcoeff_poly(Abar, ctx), ctx->fqctx);
        fq_nmod_poly_mul(modulus, leadG, trem, ctx->fqctx);
        FLINT_ASSERT(fq_nmod_poly_equal(modulus, leadA, ctx->fqctx));

        n_fq_poly_get_fq_nmod_poly(leadG, fq_nmod_mpolyn_leadcoeff_poly(G, ctx), ctx->fqctx);
        n_fq_poly_get_fq_nmod_poly(trem, fq_nmod_mpolyn_leadcoeff_poly(Bbar, ctx), ctx->fqctx);
        fq_nmod_poly_mul(modulus, leadG, trem, ctx->fqctx);
        FLINT_ASSERT(fq_nmod_poly_equal(modulus, leadB, ctx->fqctx));
    }
    fq_nmod_poly_clear(leadA, ctx->fqctx);
    fq_nmod_poly_clear(leadB, ctx->fqctx);
#endif

    fq_nmod_poly_clear(cA, ctx->fqctx);
    fq_nmod_poly_clear(cB, ctx->fqctx);
    fq_nmod_poly_clear(cG, ctx->fqctx);
    fq_nmod_poly_clear(cAbar, ctx->fqctx);
    fq_nmod_poly_clear(cBbar, ctx->fqctx);
    fq_nmod_poly_clear(trem, ctx->fqctx);
    fq_nmod_poly_clear(leadG, ctx->fqctx);
    fq_nmod_poly_clear(gamma, ctx->fqctx);

    fq_nmod_poly_clear(Aeval, ctx->fqctx);
    fq_nmod_poly_clear(Beval, ctx->fqctx);
    fq_nmod_poly_clear(Geval, ctx->fqctx);
    fq_nmod_poly_clear(Abareval, ctx->fqctx);
    fq_nmod_poly_clear(Bbareval, ctx->fqctx);

    fq_nmod_clear(gammaeval, ctx->fqctx);
    fq_nmod_clear(alpha, ctx->fqctx);
    fq_nmod_clear(temp, ctx->fqctx);

    fq_nmod_poly_clear(modulus, ctx->fqctx);
    fq_nmod_poly_clear(tempmod, ctx->fqctx);

    fq_nmod_mpolyn_clear(T, ctx);

    return success;
}


int fq_nmod_mpolyn_gcd_brown_smprime(
    fq_nmod_mpolyn_t G,
    fq_nmod_mpolyn_t Abar,
    fq_nmod_mpolyn_t Bbar,
    fq_nmod_mpolyn_t A,
    fq_nmod_mpolyn_t B,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    fq_nmod_t alpha, temp, gammaeval;
    fq_nmod_mpolyn_t Aeval, Beval, Geval, Abareval, Bbareval;
    fq_nmod_mpolyn_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    fq_nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma, trem, leadG;
    fq_nmod_poly_t modulus, tempmod;
    flint_bitcnt_t bits = A->bits;
    slong N, offset, shift;
#if FLINT_WANT_ASSERT
    fq_nmod_poly_t leadA, leadB;
#endif

    FLINT_ASSERT(var > 0);

    if (var == 1)
        return fq_nmod_mpolyn_gcd_brown_smprime_bivar(G, Abar, Bbar, A, B, ctx);

    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, G->bits, ctx->minfo);

#if FLINT_WANT_ASSERT
    fq_nmod_poly_init(leadA, ctx->fqctx);
    fq_nmod_poly_init(leadB, ctx->fqctx);
    n_fq_poly_get_fq_nmod_poly(leadA, fq_nmod_mpolyn_leadcoeff_poly(A, ctx), ctx->fqctx);
    n_fq_poly_get_fq_nmod_poly(leadB, fq_nmod_mpolyn_leadcoeff_poly(B, ctx), ctx->fqctx);
#endif

    fq_nmod_init(gammaeval, ctx->fqctx);
    fq_nmod_init(alpha, ctx->fqctx);
    fq_nmod_init(temp, ctx->fqctx);

    fq_nmod_poly_init(cAbar, ctx->fqctx);
    fq_nmod_poly_init(cBbar, ctx->fqctx);
    fq_nmod_poly_init(trem, ctx->fqctx);
    fq_nmod_poly_init(leadG, ctx->fqctx);
    fq_nmod_poly_init(modulus, ctx->fqctx);
    fq_nmod_poly_init(tempmod, ctx->fqctx);
    fq_nmod_poly_init(gamma, ctx->fqctx);
    fq_nmod_poly_init(cG, ctx->fqctx);
    fq_nmod_poly_init(cA, ctx->fqctx);
    fq_nmod_poly_init(cB, ctx->fqctx);

    fq_nmod_mpolyn_init(Aeval, bits, ctx);
    fq_nmod_mpolyn_init(Beval, bits, ctx);
    fq_nmod_mpolyn_init(Geval, bits, ctx);
    fq_nmod_mpolyn_init(Abareval, bits, ctx);
    fq_nmod_mpolyn_init(Bbareval, bits, ctx);
    fq_nmod_mpolyn_init(T, bits, ctx);

    fq_nmod_mpolyn_content_poly(cA, A, ctx);
    fq_nmod_mpolyn_content_poly(cB, B, ctx);
    fq_nmod_mpolyn_divexact_poly(A, A, cA, ctx, tempmod, trem);
    fq_nmod_mpolyn_divexact_poly(B, B, cB, ctx, tempmod, trem);

    fq_nmod_poly_gcd(cG, cA, cB, ctx->fqctx);

    fq_nmod_poly_divrem(cAbar, trem, cA, cG, ctx->fqctx);
    FLINT_ASSERT(fq_nmod_poly_is_zero(trem, ctx->fqctx));
    fq_nmod_poly_divrem(cBbar, trem, cB, cG, ctx->fqctx);
    FLINT_ASSERT(fq_nmod_poly_is_zero(trem, ctx->fqctx));

    n_fq_poly_get_fq_nmod_poly(trem, fq_nmod_mpolyn_leadcoeff_poly(A, ctx), ctx->fqctx);
    n_fq_poly_get_fq_nmod_poly(tempmod, fq_nmod_mpolyn_leadcoeff_poly(B, ctx), ctx->fqctx);
    fq_nmod_poly_gcd(gamma, trem, tempmod, ctx->fqctx);

    ldegA = fq_nmod_mpolyn_lastdeg(A, ctx);
    ldegB = fq_nmod_mpolyn_lastdeg(B, ctx);
    deggamma = fq_nmod_poly_degree(gamma, ctx->fqctx);

    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    fq_nmod_poly_one(modulus, ctx->fqctx);
    fq_nmod_poly_gen(tempmod, ctx->fqctx);
    fq_nmod_poly_neg(tempmod, tempmod, ctx->fqctx);

    fq_nmod_set_ui(alpha, 0, ctx->fqctx);

choose_prime:

    if (fq_nmod_next(alpha, ctx->fqctx) == 0)
    {
        success = 0;
        goto cleanup;
    }

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    fq_nmod_poly_evaluate_fq_nmod(gammaeval, gamma, alpha, ctx->fqctx);
    if (fq_nmod_is_zero(gammaeval, ctx->fqctx))
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    fq_nmod_mpolyn_interp_reduce_sm_mpolyn(Aeval, A, var, alpha, ctx);
    fq_nmod_mpolyn_interp_reduce_sm_mpolyn(Beval, B, var, alpha, ctx);
    FLINT_ASSERT(Aeval->length > 0);
    FLINT_ASSERT(Beval->length > 0);

    success = fq_nmod_mpolyn_gcd_brown_smprime(Geval, Abareval, Bbareval,
                                                   Aeval, Beval, var - 1, ctx);
    if (success == 0)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(Abareval->length > 0);
    FLINT_ASSERT(Bbareval->length > 0);

    if (fq_nmod_mpolyn_is_nonzero_fq_nmod(Geval, ctx))
    {
        fq_nmod_mpolyn_one(G, ctx);
        fq_nmod_mpolyn_swap(Abar, A);
        fq_nmod_mpolyn_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
    {
        int cmp;
        slong k;

        FLINT_ASSERT(G->length > 0);

        k = n_fq_poly_degree(Geval->coeffs + 0);
        cmp = mpoly_monomial_cmp_nomask_extra(G->exps + N*0,
                                     Geval->exps + N*0, N, offset, k << shift);
        if (cmp < 0)
        {
            goto choose_prime;
        }
        else if (cmp > 0)
        {
            fq_nmod_poly_one(modulus, ctx->fqctx);
        }
    }

    n_fq_get_fq_nmod(temp, fq_nmod_mpolyn_leadcoeff(Geval, ctx), ctx->fqctx);
    fq_nmod_inv(temp, temp, ctx->fqctx);
    fq_nmod_mul(temp, temp, gammaeval, ctx->fqctx);
    fq_nmod_mpolyn_scalar_mul_fq_nmod(Geval, temp, ctx);

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
    {
        fq_nmod_poly_evaluate_fq_nmod(temp, modulus, alpha, ctx->fqctx);
        fq_nmod_inv(temp, temp, ctx->fqctx);
        fq_nmod_poly_scalar_mul_fq_nmod(modulus, modulus, temp, ctx->fqctx);
        fq_nmod_mpolyn_interp_crt_sm_mpolyn(&ldegG, G, T, Geval, var, modulus, alpha, ctx);
        fq_nmod_mpolyn_interp_crt_sm_mpolyn(&ldegAbar, Abar, T, Abareval, var, modulus, alpha, ctx);
        fq_nmod_mpolyn_interp_crt_sm_mpolyn(&ldegBbar, Bbar, T, Bbareval, var, modulus, alpha, ctx);
    }
    else
    {
        fq_nmod_mpolyn_interp_lift_sm_mpolyn(G, Geval, var, ctx);
        fq_nmod_mpolyn_interp_lift_sm_mpolyn(Abar, Abareval, var, ctx);
        fq_nmod_mpolyn_interp_lift_sm_mpolyn(Bbar, Bbareval, var, ctx);
        ldegG = 0;
        ldegAbar = 0;
        ldegBbar = 0;
    }
    fq_nmod_poly_set_coeff(tempmod, 0, alpha, ctx->fqctx);
    fq_nmod_poly_mul(modulus, modulus, tempmod, ctx->fqctx);

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) < bound)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (   deggamma + ldegA == ldegG + ldegAbar
        && deggamma + ldegB == ldegG + ldegBbar )
    {
        goto successful;
    }

    fq_nmod_poly_one(modulus, ctx->fqctx);
    goto choose_prime;

successful:

    fq_nmod_mpolyn_content_poly(modulus, G, ctx);
    fq_nmod_mpolyn_divexact_poly(G, G, modulus, ctx, tempmod, trem);

    n_fq_poly_get_fq_nmod_poly(leadG, fq_nmod_mpolyn_leadcoeff_poly(G, ctx), ctx->fqctx);

    fq_nmod_mpolyn_divexact_poly(Abar, Abar, leadG, ctx, tempmod, trem);
    fq_nmod_mpolyn_divexact_poly(Bbar, Bbar, leadG, ctx, tempmod, trem);

successful_put_content:

    fq_nmod_mpolyn_mul_poly(G, G, cG, ctx, tempmod);
    fq_nmod_mpolyn_mul_poly(Abar, Abar, cAbar, ctx, tempmod);
    fq_nmod_mpolyn_mul_poly(Bbar, Bbar, cBbar, ctx, tempmod);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success)
    {
        n_fq_poly_get_fq_nmod_poly(leadG, fq_nmod_mpolyn_leadcoeff_poly(G, ctx), ctx->fqctx);
        n_fq_poly_get_fq_nmod_poly(trem, fq_nmod_mpolyn_leadcoeff_poly(Abar, ctx), ctx->fqctx);
        fq_nmod_poly_mul(modulus, leadG, trem, ctx->fqctx);
        FLINT_ASSERT(fq_nmod_poly_equal(modulus, leadA, ctx->fqctx));

        n_fq_poly_get_fq_nmod_poly(leadG, fq_nmod_mpolyn_leadcoeff_poly(G, ctx), ctx->fqctx);
        n_fq_poly_get_fq_nmod_poly(trem, fq_nmod_mpolyn_leadcoeff_poly(Bbar, ctx), ctx->fqctx);
        fq_nmod_poly_mul(modulus, leadG, trem, ctx->fqctx);
        FLINT_ASSERT(fq_nmod_poly_equal(modulus, leadB, ctx->fqctx));
    }
    fq_nmod_poly_clear(leadA, ctx->fqctx);
    fq_nmod_poly_clear(leadB, ctx->fqctx);
#endif

    fq_nmod_clear(gammaeval, ctx->fqctx);
    fq_nmod_clear(alpha, ctx->fqctx);
    fq_nmod_clear(temp, ctx->fqctx);

    fq_nmod_poly_clear(cA, ctx->fqctx);
    fq_nmod_poly_clear(cB, ctx->fqctx);
    fq_nmod_poly_clear(cG, ctx->fqctx);
    fq_nmod_poly_clear(cAbar, ctx->fqctx);
    fq_nmod_poly_clear(cBbar, ctx->fqctx);
    fq_nmod_poly_clear(trem, ctx->fqctx);
    fq_nmod_poly_clear(leadG, ctx->fqctx);
    fq_nmod_poly_clear(gamma, ctx->fqctx);
    fq_nmod_poly_clear(modulus, ctx->fqctx);
    fq_nmod_poly_clear(tempmod, ctx->fqctx);

    fq_nmod_mpolyn_clear(Aeval, ctx);
    fq_nmod_mpolyn_clear(Beval, ctx);
    fq_nmod_mpolyn_clear(Geval, ctx);
    fq_nmod_mpolyn_clear(Abareval, ctx);
    fq_nmod_mpolyn_clear(Bbareval, ctx);
    fq_nmod_mpolyn_clear(T, ctx);

    return success;
}


int fq_nmod_mpolyn_gcd_brown_lgprime_bivar(
    fq_nmod_mpolyn_t G,
    fq_nmod_mpolyn_t Abar,
    fq_nmod_mpolyn_t Bbar,
    fq_nmod_mpolyn_t A,
    fq_nmod_mpolyn_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    fq_nmod_t temp, gammaeval;
    fq_nmod_poly_t Aeval, Beval, Geval, Abareval, Bbareval;
    fq_nmod_mpolyn_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    fq_nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma, tquo, trem, leadG;
    fq_nmod_poly_t modulus, tempmod;
    flint_rand_t randstate;
    bad_fq_nmod_mpoly_embed_chooser_t embc;
    bad_fq_nmod_embed_struct * cur_emb;
    fq_nmod_mpoly_ctx_t ectx;
    slong N, off, shift;
#if FLINT_WANT_ASSERT
    fq_nmod_poly_t leadA, leadB;
#endif

#if FLINT_WANT_ASSERT
    fq_nmod_poly_init(leadA, ctx->fqctx);
    fq_nmod_poly_init(leadB, ctx->fqctx);
    n_fq_poly_get_fq_nmod_poly(leadA, fq_nmod_mpolyn_leadcoeff_poly(A, ctx), ctx->fqctx);
    n_fq_poly_get_fq_nmod_poly(leadB, fq_nmod_mpolyn_leadcoeff_poly(B, ctx), ctx->fqctx);
#endif

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, A->bits, ctx->minfo);

    fq_nmod_poly_init(tquo, ctx->fqctx);
    fq_nmod_poly_init(trem, ctx->fqctx);
    fq_nmod_poly_init(leadG, ctx->fqctx);
    fq_nmod_poly_init(cA, ctx->fqctx);
    fq_nmod_poly_init(cB, ctx->fqctx);
    fq_nmod_poly_init(cG, ctx->fqctx);
    fq_nmod_poly_init(cAbar, ctx->fqctx);
    fq_nmod_poly_init(cBbar, ctx->fqctx);
    fq_nmod_poly_init(gamma, ctx->fqctx);
    fq_nmod_poly_init(modulus, ctx->fqctx);
    fq_nmod_poly_init(tempmod, ctx->fqctx);

    fq_nmod_mpolyn_init(T, A->bits, ctx);

    fq_nmod_mpolyn_content_poly(cA, A, ctx);
    fq_nmod_mpolyn_content_poly(cB, B, ctx);
    fq_nmod_mpolyn_divexact_poly(A, A, cA, ctx, tquo, trem);
    fq_nmod_mpolyn_divexact_poly(B, B, cB, ctx, tquo, trem);

    fq_nmod_poly_gcd(cG, cA, cB, ctx->fqctx);

    fq_nmod_poly_divrem(cAbar, trem, cA, cG, ctx->fqctx);
    FLINT_ASSERT(fq_nmod_poly_is_zero(trem, ctx->fqctx));
    fq_nmod_poly_divrem(cBbar, trem, cB, cG, ctx->fqctx);
    FLINT_ASSERT(fq_nmod_poly_is_zero(trem, ctx->fqctx));

    n_fq_poly_get_fq_nmod_poly(trem, fq_nmod_mpolyn_leadcoeff_poly(A, ctx), ctx->fqctx);
    n_fq_poly_get_fq_nmod_poly(tempmod, fq_nmod_mpolyn_leadcoeff_poly(B, ctx), ctx->fqctx);
    fq_nmod_poly_gcd(gamma, trem, tempmod, ctx->fqctx);

    ldegA = fq_nmod_mpolyn_lastdeg(A, ctx);
    ldegB = fq_nmod_mpolyn_lastdeg(B, ctx);
    deggamma = fq_nmod_poly_degree(gamma, ctx->fqctx);

    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    fq_nmod_poly_one(modulus, ctx->fqctx);

    flint_randinit(randstate);
    cur_emb = bad_fq_nmod_mpoly_embed_chooser_init(embc, ectx, ctx, randstate);

    /*
        Once Aeval, Beval, ..., t are inited in ectx->fqctx, they do not need
        to be cleared and reinited when ectx->fqctx changes. Thanks fq_nmod!
    */
    fq_nmod_poly_init(Aeval, ectx->fqctx);
    fq_nmod_poly_init(Beval, ectx->fqctx);
    fq_nmod_poly_init(Geval, ectx->fqctx);
    fq_nmod_poly_init(Abareval, ectx->fqctx);
    fq_nmod_poly_init(Bbareval, ectx->fqctx);
    fq_nmod_init(gammaeval, ectx->fqctx);
    fq_nmod_init(temp, ectx->fqctx);

    /* initialization already picked a prime */
    goto have_prime;

choose_prime:

    cur_emb = bad_fq_nmod_mpoly_embed_chooser_next(embc, ectx, ctx, randstate);
    if (cur_emb == NULL)
    {
        /* ran out of primes */
        success = 0;
        goto cleanup;
    }

have_prime:

    /* make sure reduction does not kill both lc */
    bad_fq_nmod_embed_sm_to_lg(gammaeval, gamma, cur_emb);
    if (fq_nmod_is_zero(gammaeval, ectx->fqctx))
    {
        goto choose_prime;
    }

    /* make sure reduction does not kill either A or B */
    fq_nmod_mpolyn_interp_reduce_lg_poly(Aeval, ectx, A, ctx, cur_emb);
    fq_nmod_mpolyn_interp_reduce_lg_poly(Beval, ectx, B, ctx, cur_emb);
    if (Aeval->length == 0 || Beval->length == 0)
    {
        goto choose_prime;
    }

    fq_nmod_poly_gcd(Geval, Aeval, Beval, ectx->fqctx);
    success = fq_nmod_poly_divides(Abareval, Aeval, Geval, ectx->fqctx);
    FLINT_ASSERT(success);
    success = fq_nmod_poly_divides(Bbareval, Beval, Geval, ectx->fqctx);
    FLINT_ASSERT(success);

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(Abareval->length > 0);
    FLINT_ASSERT(Bbareval->length > 0);

    if (fq_nmod_poly_degree(Geval, ectx->fqctx) == 0)
    {
        fq_nmod_mpolyn_one(G, ctx);
        fq_nmod_mpolyn_swap(Abar, A);
        fq_nmod_mpolyn_swap(Bbar, B);
        goto successful_put_content;
    }

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
    {
        slong Gdeg;
        FLINT_ASSERT(G->length > 0);
        Gdeg = (G->exps + N*0)[off]>>shift;
        if (fq_nmod_poly_degree(Geval, ectx->fqctx) > Gdeg)
        {
            goto choose_prime;
        }
        else if (fq_nmod_poly_degree(Geval, ectx->fqctx) < Gdeg)
        {
            fq_nmod_poly_one(modulus, ctx->fqctx);
        }
    }

    FLINT_ASSERT(fq_nmod_is_one(Geval->coeffs + Geval->length - 1, ectx->fqctx));
    fq_nmod_poly_scalar_mul_fq_nmod(Geval, Geval, gammaeval, ectx->fqctx);

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
    {
        fq_nmod_mpolyn_interp_crt_lg_poly(&ldegG, G, T, modulus, ctx,
                                                         Geval, ectx, cur_emb);
        fq_nmod_mpolyn_interp_crt_lg_poly(&ldegAbar, Abar, T, modulus, ctx,
                                                      Abareval, ectx, cur_emb);
        fq_nmod_mpolyn_interp_crt_lg_poly(&ldegBbar, Bbar, T, modulus, ctx,
                                                      Bbareval, ectx, cur_emb);
    }
    else
    {
        fq_nmod_mpolyn_interp_lift_lg_poly(&ldegG, G, ctx, Geval, ectx, cur_emb);
        fq_nmod_mpolyn_interp_lift_lg_poly(&ldegAbar, Abar, ctx,
                                                      Abareval, ectx, cur_emb);
        fq_nmod_mpolyn_interp_lift_lg_poly(&ldegBbar, Bbar, ctx,
                                                      Bbareval, ectx, cur_emb);
    }

    fq_nmod_poly_mul(modulus, modulus, cur_emb->h, ctx->fqctx);

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) < bound)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (   deggamma + ldegA == ldegG + ldegAbar
        && deggamma + ldegB == ldegG + ldegBbar )
    {
        goto successful;
    }

    fq_nmod_poly_one(modulus, ctx->fqctx);
    goto choose_prime;

successful:

    fq_nmod_mpolyn_content_poly(modulus, G, ctx);
    fq_nmod_mpolyn_divexact_poly(G, G, modulus, ctx, tquo, trem);
    n_fq_poly_get_fq_nmod_poly(leadG, fq_nmod_mpolyn_leadcoeff_poly(G, ctx), ctx->fqctx);
    fq_nmod_mpolyn_divexact_poly(Abar, Abar, leadG, ctx, tquo, trem);
    fq_nmod_mpolyn_divexact_poly(Bbar, Bbar, leadG, ctx, tquo, trem);

successful_put_content:

    fq_nmod_mpolyn_mul_poly(G, G, cG, ctx, tquo);
    fq_nmod_mpolyn_mul_poly(Abar, Abar, cAbar, ctx, tquo);
    fq_nmod_mpolyn_mul_poly(Bbar, Bbar, cBbar, ctx, tquo);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success)
    {
        n_fq_poly_get_fq_nmod_poly(leadG, fq_nmod_mpolyn_leadcoeff_poly(G, ctx), ctx->fqctx);
        n_fq_poly_get_fq_nmod_poly(trem, fq_nmod_mpolyn_leadcoeff_poly(Abar, ctx), ctx->fqctx);
        fq_nmod_poly_mul(modulus, leadG, trem, ctx->fqctx);
        FLINT_ASSERT(fq_nmod_poly_equal(modulus, leadA, ctx->fqctx));

        n_fq_poly_get_fq_nmod_poly(leadG, fq_nmod_mpolyn_leadcoeff_poly(G, ctx), ctx->fqctx);
        n_fq_poly_get_fq_nmod_poly(trem, fq_nmod_mpolyn_leadcoeff_poly(Bbar, ctx), ctx->fqctx);
        fq_nmod_poly_mul(modulus, leadG, trem, ctx->fqctx);
        FLINT_ASSERT(fq_nmod_poly_equal(modulus, leadB, ctx->fqctx));
    }
    fq_nmod_poly_clear(leadA, ctx->fqctx);
    fq_nmod_poly_clear(leadB, ctx->fqctx);
#endif

    fq_nmod_poly_clear(cA, ctx->fqctx);
    fq_nmod_poly_clear(cB, ctx->fqctx);
    fq_nmod_poly_clear(cG, ctx->fqctx);
    fq_nmod_poly_clear(cAbar, ctx->fqctx);
    fq_nmod_poly_clear(cBbar, ctx->fqctx);
    fq_nmod_poly_clear(tquo, ctx->fqctx);
    fq_nmod_poly_clear(trem, ctx->fqctx);
    fq_nmod_poly_clear(leadG, ctx->fqctx);
    fq_nmod_poly_clear(gamma, ctx->fqctx);
    fq_nmod_poly_clear(modulus, ctx->fqctx);
    fq_nmod_poly_clear(tempmod, ctx->fqctx);

    fq_nmod_mpolyn_clear(T, ctx);

    fq_nmod_poly_clear(Aeval, ectx->fqctx);
    fq_nmod_poly_clear(Beval, ectx->fqctx);
    fq_nmod_poly_clear(Geval, ectx->fqctx);
    fq_nmod_poly_clear(Abareval, ectx->fqctx);
    fq_nmod_poly_clear(Bbareval, ectx->fqctx);
    fq_nmod_clear(gammaeval, ectx->fqctx);
    fq_nmod_clear(temp, ectx->fqctx);

    bad_fq_nmod_mpoly_embed_chooser_clear(embc, ectx, ctx, randstate);

    flint_randclear(randstate);

    return success;
}


int fq_nmod_mpolyn_gcd_brown_lgprime(
    fq_nmod_mpolyn_t G,
    fq_nmod_mpolyn_t Abar,
    fq_nmod_mpolyn_t Bbar,
    fq_nmod_mpolyn_t A,
    fq_nmod_mpolyn_t B,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    fq_nmod_t temp, gammaeval;
    fq_nmod_mpolyn_t Aeval, Beval, Geval, Abareval, Bbareval;
    fq_nmod_mpolyn_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    fq_nmod_poly_t cA, cB, cG, cAbar, cBbar, gamma, tquo, trem, leadG;
    fq_nmod_poly_t modulus, tempmod;
    flint_rand_t randstate;
    bad_fq_nmod_mpoly_embed_chooser_t embc;
    bad_fq_nmod_embed_struct * cur_emb;
    fq_nmod_mpoly_ctx_t ectx;
    flint_bitcnt_t bits = A->bits;
    slong N, offset, shift;
#if FLINT_WANT_ASSERT
    fq_nmod_poly_t leadA, leadB;
#endif

    FLINT_ASSERT(var > 0);

    if (var == 1)
        return fq_nmod_mpolyn_gcd_brown_lgprime_bivar(G, Abar, Bbar, A, B, ctx);

    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&offset, &shift, var - 1, bits, ctx->minfo);

#if FLINT_WANT_ASSERT
    fq_nmod_poly_init(leadA, ctx->fqctx);
    fq_nmod_poly_init(leadB, ctx->fqctx);
    n_fq_poly_get_fq_nmod_poly(leadA, fq_nmod_mpolyn_leadcoeff_poly(A, ctx), ctx->fqctx);
    n_fq_poly_get_fq_nmod_poly(leadB, fq_nmod_mpolyn_leadcoeff_poly(B, ctx), ctx->fqctx);
#endif

    fq_nmod_poly_init(cA, ctx->fqctx);
    fq_nmod_poly_init(cB, ctx->fqctx);
    fq_nmod_poly_init(cG, ctx->fqctx);
    fq_nmod_poly_init(cAbar, ctx->fqctx);
    fq_nmod_poly_init(cBbar, ctx->fqctx);
    fq_nmod_poly_init(tquo, ctx->fqctx);
    fq_nmod_poly_init(trem, ctx->fqctx);
    fq_nmod_poly_init(leadG, ctx->fqctx);
    fq_nmod_poly_init(gamma, ctx->fqctx);
    fq_nmod_poly_init(modulus, ctx->fqctx);
    fq_nmod_poly_init(tempmod, ctx->fqctx);

    fq_nmod_mpolyn_init(T, bits, ctx);

    fq_nmod_mpolyn_content_poly(cA, A, ctx);
    fq_nmod_mpolyn_content_poly(cB, B, ctx);
    fq_nmod_mpolyn_divexact_poly(A, A, cA, ctx, tquo, trem);
    fq_nmod_mpolyn_divexact_poly(B, B, cB, ctx, tquo, trem);

    fq_nmod_poly_gcd(cG, cA, cB, ctx->fqctx);

    fq_nmod_poly_divrem(cAbar, trem, cA, cG, ctx->fqctx);
    FLINT_ASSERT(fq_nmod_poly_is_zero(trem, ctx->fqctx));
    fq_nmod_poly_divrem(cBbar, trem, cB, cG, ctx->fqctx);
    FLINT_ASSERT(fq_nmod_poly_is_zero(trem, ctx->fqctx));

    n_fq_poly_get_fq_nmod_poly(trem, fq_nmod_mpolyn_leadcoeff_poly(A, ctx), ctx->fqctx);
    n_fq_poly_get_fq_nmod_poly(tempmod, fq_nmod_mpolyn_leadcoeff_poly(B, ctx), ctx->fqctx);
    fq_nmod_poly_gcd(gamma, trem, tempmod, ctx->fqctx);

    ldegA = fq_nmod_mpolyn_lastdeg(A, ctx);
    ldegB = fq_nmod_mpolyn_lastdeg(B, ctx);
    deggamma = fq_nmod_poly_degree(gamma, ctx->fqctx);

    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    fq_nmod_poly_one(modulus, ctx->fqctx);

    flint_randinit(randstate);
    cur_emb = bad_fq_nmod_mpoly_embed_chooser_init(embc, ectx, ctx, randstate);

    /*
        Once Aeval, Beval, ..., t are inited in ectx->fqctx, they do not need
        to be cleared and reinited when ectx->fqctx changes. Thanks fq_nmod!
    */
    fq_nmod_mpolyn_init(Aeval, bits, ectx);
    fq_nmod_mpolyn_init(Beval, bits, ectx);
    fq_nmod_mpolyn_init(Geval, bits, ectx);
    fq_nmod_mpolyn_init(Abareval, bits, ectx);
    fq_nmod_mpolyn_init(Bbareval, bits, ectx);
    fq_nmod_init(gammaeval, ectx->fqctx);
    fq_nmod_init(temp, ectx->fqctx);

    /* initialization already picked a prime */
    goto have_prime;

choose_prime: /* prime is irreducible element of Fq[v] (cur_emb->h) */

    cur_emb = bad_fq_nmod_mpoly_embed_chooser_next(embc, ectx, ctx, randstate);
    if (cur_emb == NULL)
    {
        /* ran out of primes */
        success = 0;
        goto cleanup;
    }

have_prime:

    /* make sure reduction does not kill both lc */
    bad_fq_nmod_embed_sm_to_lg(gammaeval, gamma, cur_emb);
    if (fq_nmod_is_zero(gammaeval, ectx->fqctx))
    {
        goto choose_prime;
    }

    /* make sure reduction does not kill either A or B */
    fq_nmod_mpolyn_interp_reduce_lg_mpolyn(Aeval, ectx, A, var, ctx, cur_emb);
    fq_nmod_mpolyn_interp_reduce_lg_mpolyn(Beval, ectx, B, var, ctx, cur_emb);
    if (Aeval->length == 0 || Beval->length == 0)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(fq_nmod_mpolyn_is_canonical(Aeval, ectx));
    FLINT_ASSERT(fq_nmod_mpolyn_is_canonical(Beval, ectx));

    success = fq_nmod_mpolyn_gcd_brown_smprime(Geval, Abareval, Bbareval,
                                                  Aeval, Beval, var - 1, ectx);
    if (success == 0)
    {
        goto choose_prime;
    }

    if (fq_nmod_mpolyn_is_nonzero_fq_nmod(Geval, ectx))
    {
        fq_nmod_mpolyn_one(G, ctx);
        fq_nmod_mpolyn_swap(Abar, A);
        fq_nmod_mpolyn_swap(Bbar, B);
        goto successful_put_content;
    }

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
    {
        /* compare leading monomials of Geval and G */
        int cmp = 0;
        slong k;

        FLINT_ASSERT(G->length > 0);

        k = n_fq_poly_degree(Geval->coeffs + 0);
        cmp = mpoly_monomial_cmp_nomask_extra(G->exps + N*0,
                                     Geval->exps + N*0, N, offset, k << shift);
        if (cmp < 0)
        {
            goto choose_prime;
        }
        else if (cmp > 0)
        {
            fq_nmod_poly_one(modulus, ctx->fqctx);
        }
    }

    n_fq_get_fq_nmod(temp, fq_nmod_mpolyn_leadcoeff(Geval, ectx), ectx->fqctx);
    fq_nmod_inv(temp, temp, ectx->fqctx);
    fq_nmod_mul(temp, temp, gammaeval, ectx->fqctx);
    fq_nmod_mpolyn_scalar_mul_fq_nmod(Geval, temp, ectx);

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) > 0)
    {
        fq_nmod_mpolyn_interp_crt_lg_mpolyn(&ldegG, G, T, modulus,
                                               var, ctx, Geval, ectx, cur_emb);
        fq_nmod_mpolyn_interp_crt_lg_mpolyn(&ldegAbar, Abar, T, modulus,
                                            var, ctx, Abareval, ectx, cur_emb);
        fq_nmod_mpolyn_interp_crt_lg_mpolyn(&ldegBbar, Bbar, T, modulus,
                                            var, ctx, Bbareval, ectx, cur_emb);
    }
    else
    {
        fq_nmod_mpolyn_interp_lift_lg_mpolyn(&ldegG, G,
                                               var, ctx, Geval, ectx, cur_emb);
        fq_nmod_mpolyn_interp_lift_lg_mpolyn(&ldegAbar, Abar,
                                            var, ctx, Abareval, ectx, cur_emb);
        fq_nmod_mpolyn_interp_lift_lg_mpolyn(&ldegBbar, Bbar,
                                            var, ctx, Bbareval, ectx, cur_emb);
    }
    fq_nmod_poly_mul(modulus, modulus, cur_emb->h, ctx->fqctx);

    if (fq_nmod_poly_degree(modulus, ctx->fqctx) < bound)
    {
        goto choose_prime;
    }

    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (   deggamma + ldegA == ldegG + ldegAbar
        && deggamma + ldegB == ldegG + ldegBbar )
    {
        goto successful;
    }

    fq_nmod_poly_one(modulus, ctx->fqctx);
    goto choose_prime;

successful:

    fq_nmod_mpolyn_content_poly(modulus, G, ctx);
    fq_nmod_mpolyn_divexact_poly(G, G, modulus, ctx, tquo, trem);
    n_fq_poly_get_fq_nmod_poly(leadG, fq_nmod_mpolyn_leadcoeff_poly(G, ctx), ctx->fqctx);
    fq_nmod_mpolyn_divexact_poly(Abar, Abar, leadG, ctx, tquo, trem);
    fq_nmod_mpolyn_divexact_poly(Bbar, Bbar, leadG, ctx, tquo, trem);

successful_put_content:

    fq_nmod_mpolyn_mul_poly(G, G, cG, ctx, tquo);
    fq_nmod_mpolyn_mul_poly(Abar, Abar, cAbar, ctx, tquo);
    fq_nmod_mpolyn_mul_poly(Bbar, Bbar, cBbar, ctx, tquo);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success)
    {
        n_fq_poly_get_fq_nmod_poly(leadG, fq_nmod_mpolyn_leadcoeff_poly(G, ctx), ctx->fqctx);
        n_fq_poly_get_fq_nmod_poly(trem, fq_nmod_mpolyn_leadcoeff_poly(Abar, ctx), ctx->fqctx);
        fq_nmod_poly_mul(modulus, leadG, trem, ctx->fqctx);
        FLINT_ASSERT(fq_nmod_poly_equal(modulus, leadA, ctx->fqctx));

        n_fq_poly_get_fq_nmod_poly(leadG, fq_nmod_mpolyn_leadcoeff_poly(G, ctx), ctx->fqctx);
        n_fq_poly_get_fq_nmod_poly(trem, fq_nmod_mpolyn_leadcoeff_poly(Bbar, ctx), ctx->fqctx);
        fq_nmod_poly_mul(modulus, leadG, trem, ctx->fqctx);
        FLINT_ASSERT(fq_nmod_poly_equal(modulus, leadB, ctx->fqctx));
    }
    fq_nmod_poly_clear(leadA, ctx->fqctx);
    fq_nmod_poly_clear(leadB, ctx->fqctx);
#endif

    fq_nmod_poly_clear(cA, ctx->fqctx);
    fq_nmod_poly_clear(cB, ctx->fqctx);
    fq_nmod_poly_clear(cG, ctx->fqctx);
    fq_nmod_poly_clear(cAbar, ctx->fqctx);
    fq_nmod_poly_clear(cBbar, ctx->fqctx);
    fq_nmod_poly_clear(tquo, ctx->fqctx);
    fq_nmod_poly_clear(trem, ctx->fqctx);
    fq_nmod_poly_clear(leadG, ctx->fqctx);
    fq_nmod_poly_clear(gamma, ctx->fqctx);
    fq_nmod_poly_clear(modulus, ctx->fqctx);
    fq_nmod_poly_clear(tempmod, ctx->fqctx);

    fq_nmod_mpolyn_clear(T, ctx);

    fq_nmod_mpolyn_clear(Aeval, ectx);
    fq_nmod_mpolyn_clear(Beval, ectx);
    fq_nmod_mpolyn_clear(Geval, ectx);
    fq_nmod_mpolyn_clear(Abareval, ectx);
    fq_nmod_mpolyn_clear(Bbareval, ectx);
    fq_nmod_clear(gammaeval, ectx->fqctx);
    fq_nmod_clear(temp, ectx->fqctx);

    bad_fq_nmod_mpoly_embed_chooser_clear(embc, ectx, ctx, randstate);

    flint_randclear(randstate);

    return success;
}

