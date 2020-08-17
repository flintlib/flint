/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"


int fmpz_mod_mpolyn_gcd_brown_bivar(
    fmpz_mod_mpolyn_t G,
    fmpz_mod_mpolyn_t Abar,
    fmpz_mod_mpolyn_t Bbar,
    fmpz_mod_mpolyn_t A,
    fmpz_mod_mpolyn_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;
    slong bound;
    fmpz_t alpha, temp, gammaeval;
    fmpz_mod_poly_t Aeval, Beval, Geval, Abareval, Bbareval;
    fmpz_mod_mpolyn_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    fmpz_mod_poly_t cA, cB, cG, cAbar, cBbar, gamma, r;
    fmpz_mod_poly_t modulus, modulus2;
    slong N, off, shift;
    flint_bitcnt_t bits = A->bits;
#if WANT_ASSERT
    fmpz_mod_poly_t leadA, leadB;
#endif

    FLINT_ASSERT(A->bits == bits);
    FLINT_ASSERT(B->bits == bits);
    FLINT_ASSERT(G->bits == bits);
    FLINT_ASSERT(Abar->bits == bits);
    FLINT_ASSERT(Bbar->bits == bits);

    fmpz_init(alpha);
    fmpz_init(temp);
    fmpz_init(gammaeval);

#if WANT_ASSERT
    fmpz_mod_poly_init(leadA, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(leadB, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_set(leadA, fmpz_mod_mpolyn_leadcoeff_poly(A, ctx));
    fmpz_mod_poly_set(leadB, fmpz_mod_mpolyn_leadcoeff_poly(B, ctx));
#endif

    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, bits, ctx->minfo);

    fmpz_mod_mpolyn_init(T, bits, ctx);
    fmpz_mod_poly_init(r, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(cA, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(cB, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(cG, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(cAbar, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(cBbar, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(gamma, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(Aeval, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(Beval, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(Geval, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(Abareval, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(Bbareval, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(modulus, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_poly_init(modulus2, fmpz_mod_ctx_modulus(ctx->ffinfo));

    fmpz_mod_mpolyn_content_poly(cA, A, ctx);
    fmpz_mod_mpolyn_content_poly(cB, B, ctx);
    fmpz_mod_mpolyn_divexact_poly(A, cA, ctx);
    fmpz_mod_mpolyn_divexact_poly(B, cB, ctx);

    fmpz_mod_poly_gcd(cG, cA, cB);

    fmpz_mod_poly_divrem(cAbar, r, cA, cG);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r));
    fmpz_mod_poly_divrem(cBbar, r, cB, cG);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r));

    fmpz_mod_poly_gcd(gamma, fmpz_mod_mpolyn_leadcoeff_poly(A, ctx),
                             fmpz_mod_mpolyn_leadcoeff_poly(B, ctx));

    ldegA = fmpz_mod_mpolyn_lastdeg(A, ctx);
    ldegB = fmpz_mod_mpolyn_lastdeg(B, ctx);
    deggamma = fmpz_mod_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    fmpz_mod_poly_set_ui(modulus, 1);

    fmpz_sub_ui(alpha, fmpz_mod_ctx_modulus(ctx->ffinfo), 1);

choose_prime:

    fmpz_sub_ui(alpha, alpha, 1);
    if (fmpz_sgn(alpha) <= 0)
    {
        success = 0;
        goto cleanup;
    }

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    fmpz_mod_poly_evaluate_fmpz(gammaeval, gamma, alpha);
    if (fmpz_is_zero(gammaeval))
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    fmpz_mod_mpolyn_intp_reduce_sm_poly(Aeval, A, alpha, ctx);
    fmpz_mod_mpolyn_intp_reduce_sm_poly(Beval, B, alpha, ctx);
    FLINT_ASSERT(Aeval->length > 0);
    FLINT_ASSERT(Beval->length > 0);

    fmpz_mod_poly_gcd(Geval, Aeval, Beval);
    fmpz_mod_poly_divrem(Abareval, r, Aeval, Geval);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r));
    fmpz_mod_poly_divrem(Bbareval, r, Beval, Geval);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r));

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(Abareval->length > 0);
    FLINT_ASSERT(Bbareval->length > 0);

    if (fmpz_mod_poly_degree(Geval) == 0)
    {
        fmpz_mod_mpolyn_one(G, ctx);
        fmpz_mod_mpolyn_swap(Abar, A, ctx);
        fmpz_mod_mpolyn_swap(Bbar, B, ctx);
        goto successful_put_content;    
    }

    if (fmpz_mod_poly_degree(modulus) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (fmpz_mod_poly_degree(Geval) > ((G->exps + N*0)[off]>>shift))
        {
            goto choose_prime;
        }
        else if (fmpz_mod_poly_degree(Geval) < ((G->exps + N*0)[off]>>shift))
        {
            fmpz_mod_poly_set_ui(modulus, 1);
        }
    }

    fmpz_mod_poly_scalar_mul_fmpz(Geval, Geval, gammaeval);

    if (fmpz_mod_poly_degree(modulus) > 0)
    {
        fmpz_mod_poly_evaluate_fmpz(temp, modulus, alpha);
        fmpz_mod_inv(temp, temp, ctx->ffinfo);
        fmpz_mod_poly_scalar_mul_fmpz(modulus, modulus, temp);
        fmpz_mod_mpolyn_intp_crt_sm_poly(&ldegG, G, T, Geval, modulus, alpha, ctx);
        fmpz_mod_mpolyn_intp_crt_sm_poly(&ldegAbar, Abar, T, Abareval, modulus, alpha, ctx);
        fmpz_mod_mpolyn_intp_crt_sm_poly(&ldegBbar, Bbar, T, Bbareval, modulus, alpha, ctx);
    }
    else
    {
        fmpz_mod_mpolyn_intp_lift_sm_poly(G, Geval, ctx);
        fmpz_mod_mpolyn_intp_lift_sm_poly(Abar, Abareval, ctx);
        fmpz_mod_mpolyn_intp_lift_sm_poly(Bbar, Bbareval, ctx);
        ldegG = 0;
        ldegAbar = 0;
        ldegBbar = 0;
    }

    fmpz_mod_poly_scalar_mul_fmpz(modulus2, modulus, alpha);
    fmpz_mod_poly_shift_left(modulus, modulus, 1);
    fmpz_mod_poly_sub(modulus, modulus, modulus2);

    if (fmpz_mod_poly_degree(modulus) < bound)
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

    fmpz_mod_poly_set_ui(modulus, 1);
    goto choose_prime;

successful:

    fmpz_mod_mpolyn_content_poly(modulus, G, ctx);
    fmpz_mod_mpolyn_divexact_poly(G, modulus, ctx);
    fmpz_mod_mpolyn_divexact_poly(Abar, fmpz_mod_mpolyn_leadcoeff_poly(G, ctx), ctx);
    fmpz_mod_mpolyn_divexact_poly(Bbar, fmpz_mod_mpolyn_leadcoeff_poly(G, ctx), ctx);

successful_put_content:

    fmpz_mod_mpolyn_mul_poly(G, cG, ctx);
    fmpz_mod_mpolyn_mul_poly(Abar, cAbar, ctx);
    fmpz_mod_mpolyn_mul_poly(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(fmpz_is_one(fmpz_mod_mpolyn_leadcoeff_last_ref(G, ctx)));
        fmpz_mod_poly_mul(modulus, fmpz_mod_mpolyn_leadcoeff_poly(G, ctx),
                                   fmpz_mod_mpolyn_leadcoeff_poly(Abar, ctx));
        FLINT_ASSERT(fmpz_mod_poly_equal(modulus, leadA));
        fmpz_mod_poly_mul(modulus, fmpz_mod_mpolyn_leadcoeff_poly(G, ctx),
                                   fmpz_mod_mpolyn_leadcoeff_poly(Bbar, ctx));
        FLINT_ASSERT(fmpz_mod_poly_equal(modulus, leadB));
    }
    fmpz_mod_poly_clear(leadA);
    fmpz_mod_poly_clear(leadB);
#endif

    fmpz_mod_poly_clear(r);
    fmpz_mod_poly_clear(cA);
    fmpz_mod_poly_clear(cB);
    fmpz_mod_poly_clear(cG);
    fmpz_mod_poly_clear(cAbar);
    fmpz_mod_poly_clear(cBbar);

    fmpz_mod_poly_clear(gamma);

    fmpz_mod_poly_clear(Aeval);
    fmpz_mod_poly_clear(Beval);
    fmpz_mod_poly_clear(Geval);
    fmpz_mod_poly_clear(Abareval);
    fmpz_mod_poly_clear(Bbareval);

    fmpz_mod_mpolyn_clear(T, ctx);

    fmpz_mod_poly_clear(modulus);
    fmpz_mod_poly_clear(modulus2);

    fmpz_clear(alpha);
    fmpz_clear(temp);
    fmpz_clear(gammaeval);

    return success;
}

