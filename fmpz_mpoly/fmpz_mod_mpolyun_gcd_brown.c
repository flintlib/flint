/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"
#include "fmpz_mpoly.h"

/*
    E = A(v = alpha)
    A is in R[X][v]
    E is in R[X]
*/
void fmpz_mod_mpolyun_eval_last_bivar(
    fmpz_mod_poly_t E,
    const fmpz_mod_mpolyun_t A,
    const fmpz_t alpha,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    fmpz_t v;
    slong Ai, Alen;
    fmpz_mod_mpolyn_struct * Acoeff;
    ulong * Aexp;

    fmpz_init(v);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = A->length;
    Ai = 0;
    fmpz_mod_poly_zero(E);
    for (Ai = 0; Ai < Alen; Ai++)
    {
        fmpz_mod_poly_evaluate_fmpz(v, (Acoeff + Ai)->coeffs + 0, alpha);
        fmpz_mod_poly_set_coeff_fmpz(E, Aexp[Ai], v);
    }

    fmpz_clear(v);
}

/*
    A = B
    A, B are in R[X]
*/
void fmpz_mod_mpolyun_startinterp_bivar(
    fmpz_mod_mpolyun_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    slong Bexp;
    slong Blen = fmpz_mod_poly_length(B);
    fmpz * Bcoeff = B->coeffs;
    fmpz_mod_mpolyn_struct * Acoeff;
    ulong * Aexp;
    slong Ai;

    fmpz_mod_mpolyun_fit_length(A, Blen, ctx, fpctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    Ai = 0;
    for (Bexp = Blen - 1; Bexp >= 0; Bexp--)
    {
        if (!fmpz_is_zero(Bcoeff + Bexp))
        {
            FLINT_ASSERT(Ai < A->alloc);

            fmpz_mod_mpolyn_fit_length(Acoeff + Ai, 1, ctx, fpctx);
            mpoly_monomial_zero((Acoeff + Ai)->exps + N*0, N);
            fmpz_mod_poly_zero((Acoeff + Ai)->coeffs + 0);
            fmpz_mod_poly_set_coeff_fmpz((Acoeff + Ai)->coeffs + 0, 0, Bcoeff + Bexp);
            Aexp[Ai] = Bexp;
            (Acoeff + Ai)->length = 1;
            Ai++;
        }
    }
    A->length = Ai;
}


/*
    F = F + modulus*(A - F(v = alpha))
    no assumptions about matching monomials
    F is in Fp[X][v]
    A is in Fp[X]
    it is expected that modulus(alpha) == 1
*/
int fmpz_mod_mpolyun_addinterp_bivar(
    slong * lastdeg_,
    fmpz_mod_mpolyun_t F,
    fmpz_mod_mpolyun_t T,
    const fmpz_mod_poly_t A,
    const fmpz_mod_poly_t modulus,
    const fmpz_t alpha,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    int changed = 0;
    slong lastdeg = -WORD(1);
    slong N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    fmpz_t v;
    slong Fi, Toff, Aexp;
    fmpz * Acoeff = A->coeffs;
    slong Flen = F->length;
    fmpz_mod_mpolyn_struct * Fcoeff = F->coeffs;
    ulong * Fexp = F->exps;
    fmpz_mod_mpolyn_struct * Tcoeff;
    ulong * Texp;
    fmpz_mod_poly_t tp;
    
    Fi = 0;

    fmpz_init(v);

    Aexp = fmpz_mod_poly_degree(A);

    fmpz_mod_poly_init(tp, fmpz_mod_ctx_modulus(fpctx));

    fmpz_mod_mpolyun_fit_length(T, Flen + Aexp + 1, ctx, fpctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;
    Toff = 0;

    while (Fi < Flen || Aexp >= 0)
    {
        FLINT_ASSERT(Toff < T->alloc);

        if (Fi < Flen)
        {
            FLINT_ASSERT(!fmpz_mod_poly_is_zero((Fcoeff + Fi)->coeffs + 0));
            FLINT_ASSERT(fmpz_mod_poly_degree((Fcoeff + Fi)->coeffs + 0) < fmpz_mod_poly_degree(modulus));
        }

        if (Aexp >= 0)
        {
            FLINT_ASSERT(Acoeff[Aexp] != UWORD(0));
        }

        fmpz_mod_mpolyn_fit_length(Tcoeff + Toff, 1, ctx, fpctx);

        if (Fi < Flen && Aexp >= 0 && Fexp[Fi] == Aexp)
        {
            /* F term ok, A term ok */
            fmpz_mod_poly_evaluate_fmpz(v, (Fcoeff + Fi)->coeffs + 0, alpha);
            fmpz_mod_sub(v, Acoeff + Aexp, v, fpctx);
            if (!fmpz_is_zero(v))
            {
                changed = 1;
                fmpz_mod_poly_scalar_mul_fmpz(tp, modulus, v);
                fmpz_mod_poly_add((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0, tp);
            }
            else
            {
                fmpz_mod_poly_set((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0);
            }
            Texp[Toff] = Aexp;
            Fi++;
            do {
                Aexp--;
            } while (Aexp >= 0 && fmpz_is_zero(Acoeff + Aexp));
        }
        else if (Fi < Flen && (Aexp < 0 || Fexp[Fi] > Aexp))
        {
            /* F term ok, A term missing */
            fmpz_mod_poly_evaluate_fmpz(v, (Fcoeff + Fi)->coeffs + 0, alpha);
            if (!fmpz_is_zero(v))
            {
                changed = 1;
                fmpz_mod_poly_scalar_mul_fmpz(tp, modulus, v);
                fmpz_mod_poly_sub((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0, tp);
            }
            else
            {
                fmpz_mod_poly_set((Tcoeff + Toff)->coeffs + 0, (Fcoeff + Fi)->coeffs + 0);                
            }
            Texp[Toff] = Fexp[Fi];
            Fi++;
        }
        else if (Aexp >= 0 && (Fi >= Flen || Fexp[Fi] < Aexp))
        {
            /* F term missing, A term ok */
            changed = 1;
            fmpz_mod_poly_scalar_mul_fmpz((Tcoeff + Toff)->coeffs + 0, modulus, Acoeff + Aexp);
            Texp[Toff] = Aexp;
            do {
                Aexp--;
            } while (Aexp >= 0 && fmpz_is_zero(Acoeff + Aexp));
        }
        else
        {
            FLINT_ASSERT(0);
        }

        lastdeg = FLINT_MAX(lastdeg, fmpz_mod_poly_degree((Tcoeff + Toff)->coeffs + 0));
        mpoly_monomial_zero((Tcoeff + Toff)->exps + N*0, N);
        FLINT_ASSERT(!fmpz_mod_poly_is_zero((Tcoeff + Toff)->coeffs + 0));
        (Tcoeff + Toff)->length = 1;
        Toff++;
    }
    T->length = Toff;

    fmpz_mod_poly_clear(tp);

    if (changed)
    {
        fmpz_mod_mpolyun_swap(T, F);
    }

    fmpz_clear(v);

    *lastdeg_ = lastdeg;
    return changed;
}


int fmpz_mod_mpolyun_gcd_brown_bivar(
    fmpz_mod_mpolyun_t G,
    fmpz_mod_mpolyun_t Abar,
    fmpz_mod_mpolyun_t Bbar,
    fmpz_mod_mpolyun_t A,
    fmpz_mod_mpolyun_t B,
    const fmpz_mpoly_ctx_t ctx,
    const fmpz_mod_ctx_t fpctx)
{
    int success;
    slong bound;
    fmpz_t alpha, temp, gammaeval;
    fmpz_mod_poly_t Aeval, Beval, Geval, Abareval, Bbareval;
    fmpz_mod_mpolyun_t T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    fmpz_mod_poly_t cA, cB, cG, cAbar, cBbar, gamma, r;
    fmpz_mod_poly_t modulus, modulus2;
#if WANT_ASSERT
    fmpz_mod_poly_t leadA, leadB;
#endif

    fmpz_init(alpha);
    fmpz_init(temp);
    fmpz_init(gammaeval);
/*
printf("fmpz_mod_mpolyun_gcd_bivar called\n");
printf("A: "); fmpz_mod_mpolyun_print_pretty(A, NULL, ctx, fpctx); printf("\n");
printf("B: "); fmpz_mod_mpolyun_print_pretty(B, NULL, ctx, fpctx); printf("\n");
*/

#if WANT_ASSERT
    fmpz_mod_poly_init(leadA, fmpz_mod_ctx_modulus(fpctx));
    fmpz_mod_poly_init(leadB, fmpz_mod_ctx_modulus(fpctx));
    fmpz_mod_poly_set(leadA, fmpz_mod_mpolyun_leadcoeff_ref(A, ctx, fpctx));
    fmpz_mod_poly_set(leadB, fmpz_mod_mpolyun_leadcoeff_ref(B, ctx, fpctx));
#endif

    fmpz_mod_poly_init(r, fmpz_mod_ctx_modulus(fpctx));
    fmpz_mod_poly_init(cA, fmpz_mod_ctx_modulus(fpctx));
    fmpz_mod_poly_init(cB, fmpz_mod_ctx_modulus(fpctx));
    fmpz_mod_mpolyun_content_last(cA, A, ctx, fpctx);
    fmpz_mod_mpolyun_content_last(cB, B, ctx, fpctx);
    fmpz_mod_mpolyun_divexact_last(A, cA, ctx, fpctx);
    fmpz_mod_mpolyun_divexact_last(B, cB, ctx, fpctx);

    fmpz_mod_poly_init(cG, fmpz_mod_ctx_modulus(fpctx));
    fmpz_mod_poly_gcd_euclidean(cG, cA, cB);

    fmpz_mod_poly_init(cAbar, fmpz_mod_ctx_modulus(fpctx));
    fmpz_mod_poly_init(cBbar, fmpz_mod_ctx_modulus(fpctx));
    fmpz_mod_poly_divrem(cAbar, r, cA, cG);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r));
    fmpz_mod_poly_divrem(cBbar, r, cB, cG);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(r));

    fmpz_mod_poly_init(gamma, fmpz_mod_ctx_modulus(fpctx));
    fmpz_mod_poly_gcd(gamma, fmpz_mod_mpolyun_leadcoeff_ref(A, ctx, fpctx),
                             fmpz_mod_mpolyun_leadcoeff_ref(B, ctx, fpctx));

    ldegA = fmpz_mod_mpolyun_lastdeg(A, ctx, fpctx);
    ldegB = fmpz_mod_mpolyun_lastdeg(B, ctx, fpctx);
    deggamma = fmpz_mod_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);
/*
flint_printf("bound: %wd\n", bound);
*/
    fmpz_mod_poly_init(Aeval, fmpz_mod_ctx_modulus(fpctx));
    fmpz_mod_poly_init(Beval, fmpz_mod_ctx_modulus(fpctx));
    fmpz_mod_poly_init(Geval, fmpz_mod_ctx_modulus(fpctx));
    fmpz_mod_poly_init(Abareval, fmpz_mod_ctx_modulus(fpctx));
    fmpz_mod_poly_init(Bbareval, fmpz_mod_ctx_modulus(fpctx));

    fmpz_mod_mpolyun_init(T, A->bits, ctx, fpctx);

    fmpz_mod_poly_init(modulus, fmpz_mod_ctx_modulus(fpctx));
    fmpz_mod_poly_init(modulus2, fmpz_mod_ctx_modulus(fpctx));
    fmpz_mod_poly_set_ui(modulus, 1);

    fmpz_sub_ui(alpha, fmpz_mod_ctx_modulus(fpctx), 2);

choose_prime: /* prime is v - alpha */

    if (fmpz_cmp_ui(alpha, 2) < 1)
    {
        success = 0;
        goto cleanup;
    }

    fmpz_sub_ui(alpha, alpha, 1);

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    fmpz_mod_poly_evaluate_fmpz(gammaeval, gamma, alpha);
    if (fmpz_is_zero(gammaeval))
    {
        goto choose_prime;
    }

    /* evaluation point should kill neither A nor B */
    fmpz_mod_mpolyun_eval_last_bivar(Aeval, A, alpha, ctx, fpctx);
    fmpz_mod_mpolyun_eval_last_bivar(Beval, B, alpha, ctx, fpctx);
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
        fmpz_mod_mpolyun_one(G, ctx, fpctx);
        fmpz_mod_mpolyun_swap(Abar, A);
        fmpz_mod_mpolyun_swap(Bbar, B);
        goto successful_put_content;    
    }

    if (fmpz_mod_poly_degree(modulus) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (fmpz_mod_poly_degree(Geval) > G->exps[0])
        {
            goto choose_prime;
        }
        else if (fmpz_mod_poly_degree(Geval) < G->exps[0])
        {
            fmpz_mod_poly_set_ui(modulus, 1);
        }
    }

    fmpz_mod_poly_scalar_mul_fmpz(Geval, Geval, gammaeval);

    if (fmpz_mod_poly_degree(modulus) > 0)
    {
        fmpz_mod_poly_evaluate_fmpz(temp, modulus, alpha);
        fmpz_mod_inv(temp, temp, fpctx);
        fmpz_mod_poly_scalar_mul_fmpz(modulus, modulus, temp);
        fmpz_mod_mpolyun_addinterp_bivar(&ldegG, G, T, Geval, modulus, alpha, ctx, fpctx);
        fmpz_mod_mpolyun_addinterp_bivar(&ldegAbar, Abar, T, Abareval, modulus, alpha, ctx, fpctx);
        fmpz_mod_mpolyun_addinterp_bivar(&ldegBbar, Bbar, T, Bbareval, modulus, alpha, ctx, fpctx);
    }
    else
    {
        fmpz_mod_mpolyun_startinterp_bivar(G, Geval, ctx, fpctx);
        fmpz_mod_mpolyun_startinterp_bivar(Abar, Abareval, ctx, fpctx);
        fmpz_mod_mpolyun_startinterp_bivar(Bbar, Bbareval, ctx, fpctx);
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

    fmpz_mod_mpolyun_content_last(modulus, G, ctx, fpctx);
    fmpz_mod_mpolyun_divexact_last(G, modulus, ctx, fpctx);
    fmpz_mod_mpolyun_divexact_last(Abar, fmpz_mod_mpolyun_leadcoeff_ref(G, ctx, fpctx), ctx, fpctx);
    fmpz_mod_mpolyun_divexact_last(Bbar, fmpz_mod_mpolyun_leadcoeff_ref(G, ctx, fpctx), ctx, fpctx);

successful_put_content:

    fmpz_mod_mpolyun_mul_last(G, cG, ctx, fpctx);
    fmpz_mod_mpolyun_mul_last(Abar, cAbar, ctx, fpctx);
    fmpz_mod_mpolyun_mul_last(Bbar, cBbar, ctx, fpctx);

    success = 1;

cleanup:

#if WANT_ASSERT
    if (success)
    {
        FLINT_ASSERT(fmpz_is_one(fmpz_mod_mpolyun_leadcoeff_last_ref(G, ctx, fpctx)));
        fmpz_mod_poly_mul(modulus, fmpz_mod_mpolyun_leadcoeff_ref(G, ctx, fpctx),
                                   fmpz_mod_mpolyun_leadcoeff_ref(Abar, ctx, fpctx));
        FLINT_ASSERT(fmpz_mod_poly_equal(modulus, leadA));
        fmpz_mod_poly_mul(modulus, fmpz_mod_mpolyun_leadcoeff_ref(G, ctx, fpctx),
                                   fmpz_mod_mpolyun_leadcoeff_ref(Bbar, ctx, fpctx));
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

    fmpz_mod_mpolyun_clear(T, ctx, fpctx);

    fmpz_mod_poly_clear(modulus);
    fmpz_mod_poly_clear(modulus2);

    fmpz_clear(alpha);
    fmpz_clear(temp);
    fmpz_clear(gammaeval);

    return success;
}
