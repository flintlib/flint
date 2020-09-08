/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

/*
    E = A(v = alpha)
    A is in R[X][v]
    E is in R[X]
*/
void fmpz_mod_mpolyn_intp_reduce_sm_poly(
    fmpz_mod_poly_t E,
    const fmpz_mod_mpolyn_t A,
    const fmpz_t alpha,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t v;
    slong Ai, Alen, k;
    fmpz_mod_poly_struct * Acoeff;
    ulong * Aexp;
    slong N, off, shift;

    fmpz_init(v);

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, A->bits, ctx->minfo);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = A->length;
    Ai = 0;
    fmpz_mod_poly_zero(E, ctx->ffinfo);
    for (Ai = 0; Ai < Alen; Ai++)
    {
        fmpz_mod_poly_evaluate_fmpz(v, Acoeff + Ai, alpha, ctx->ffinfo);
        k = (Aexp + N*Ai)[off] >> shift;
        fmpz_mod_poly_set_coeff_fmpz(E, k, v, ctx->ffinfo);
    }

    fmpz_clear(v);
}

/*
    A = B
    A, B are in R[X]
*/
void fmpz_mod_mpolyn_intp_lift_sm_poly(
    fmpz_mod_mpolyn_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong Bi;
    slong Blen = fmpz_mod_poly_length(B, ctx->ffinfo);
    fmpz * Bcoeff = B->coeffs;
    fmpz_mod_poly_struct * Acoeff;
    ulong * Aexp;
    slong Ai;
    slong N, off, shift;

    N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, A->bits, ctx->minfo);

    fmpz_mod_mpolyn_fit_length(A, Blen, ctx);
    Acoeff = A->coeffs;
    Aexp = A->exps;

    Ai = 0;
    for (Bi = Blen - 1; Bi >= 0; Bi--)
    {
        if (!fmpz_is_zero(Bcoeff + Bi))
        {
            FLINT_ASSERT(Ai < A->alloc);

            fmpz_mod_poly_set_fmpz(Acoeff + Ai, Bcoeff + Bi, ctx->ffinfo);
            mpoly_monomial_zero(Aexp + N*Ai, N);
            (Aexp + N*Ai)[off] = Bi << shift;
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
int fmpz_mod_mpolyn_intp_crt_sm_poly(
    slong * lastdeg_,
    fmpz_mod_mpolyn_t F,
    fmpz_mod_mpolyn_t T,
    const fmpz_mod_poly_t A,
    const fmpz_mod_poly_t modulus,
    const fmpz_t alpha,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong lastdeg = -WORD(1);
    fmpz_t u, v;
    slong Fi, Ti, Ai;
    fmpz * Acoeff = A->coeffs;
    slong Flen = F->length;
    fmpz_mod_poly_struct * Fcoeff = F->coeffs;
    ulong * Fexp = F->exps;
    fmpz_mod_poly_struct * Tcoeff;
    ulong * Texp;
    fmpz_mod_poly_t tp;
    slong N, off, shift;

    N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, 0, F->bits, ctx->minfo);

    Fi = 0;
    Ai = fmpz_mod_poly_degree(A, ctx->ffinfo);

    fmpz_init(u);
    fmpz_init(v);
    fmpz_mod_poly_init(tp, ctx->ffinfo);

    fmpz_mod_mpolyn_fit_length(T, Flen + Ai + 1, ctx);
    Tcoeff = T->coeffs;
    Texp = T->exps;
    Ti = 0;

    while (Fi < Flen || Ai >= 0)
    {
        FLINT_ASSERT(Ti < T->alloc);

        if (Fi < Flen)
        {
            FLINT_ASSERT(!fmpz_mod_poly_is_zero(Fcoeff + Fi, ctx->ffinfo));
            FLINT_ASSERT(fmpz_mod_poly_degree(Fcoeff + Fi, ctx->ffinfo) <
                                   fmpz_mod_poly_degree(modulus, ctx->ffinfo));
        }

        if (Ai >= 0)
        {
            FLINT_ASSERT(!fmpz_is_zero(Acoeff + Ai));
        }

        mpoly_monomial_zero(Texp + N*Ti, N);

        if (Fi < Flen && Ai >= 0 && ((Fexp + N*Fi)[off]>>shift) == Ai)
        {
            /* F term ok, A term ok */
            fmpz_mod_poly_evaluate_fmpz(u, Fcoeff + Fi, alpha, ctx->ffinfo);
            fmpz_mod_sub(v, Acoeff + Ai, u, ctx->ffinfo);
            if (!fmpz_is_zero(v))
            {
                changed = 1;
                fmpz_mod_poly_scalar_mul_fmpz(tp, modulus, v, ctx->ffinfo);
                fmpz_mod_poly_add(Tcoeff + Ti, Fcoeff + Fi, tp, ctx->ffinfo);
            }
            else
            {
                fmpz_mod_poly_set(Tcoeff + Ti, Fcoeff + Fi, ctx->ffinfo);
            }

            (Texp + N*Ti)[off] = Ai << shift;
            Fi++;
            do {
                Ai--;
            } while (Ai >= 0 && fmpz_is_zero(Acoeff + Ai));
        }
        else if (Fi < Flen && (Ai < 0 || ((Fexp + N*Fi)[off]>>shift) > Ai))
        {
            /* F term ok, A term missing */
            fmpz_mod_poly_evaluate_fmpz(v, Fcoeff + Fi, alpha, ctx->ffinfo);
            if (!fmpz_is_zero(v))
            {
                changed = 1;
                fmpz_mod_poly_scalar_mul_fmpz(tp, modulus, v, ctx->ffinfo);
                fmpz_mod_poly_sub(Tcoeff + Ti, Fcoeff + Fi, tp, ctx->ffinfo);
            }
            else
            {
                fmpz_mod_poly_set(Tcoeff + Ti, Fcoeff + Fi, ctx->ffinfo);                
            }

            (Texp + N*Ti)[off] = (Fexp + N*Fi)[off];
            Fi++;
        }
        else if (Ai >= 0 && (Fi >= Flen || ((Fexp + N*Fi)[off]>>shift) < Ai))
        {
            /* F term missing, A term ok */
            changed = 1;
            fmpz_mod_poly_scalar_mul_fmpz(Tcoeff + Ti, modulus, Acoeff + Ai, ctx->ffinfo);

            (Texp + N*Ti)[off] = Ai << shift;
            do {
                Ai--;
            } while (Ai >= 0 && fmpz_is_zero(Acoeff + Ai));
        }
        else
        {
            FLINT_ASSERT(0);
        }

        FLINT_ASSERT(!fmpz_mod_poly_is_zero(Tcoeff + Ti, ctx->ffinfo));
        lastdeg = FLINT_MAX(lastdeg, fmpz_mod_poly_degree(Tcoeff + Ti, ctx->ffinfo));

        Ti++;
    }
    T->length = Ti;

    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_mod_poly_clear(tp, ctx->ffinfo);

    if (changed)
    {
        fmpz_mod_mpolyn_swap(T, F, ctx);
    }

    *lastdeg_ = lastdeg;
    return changed;
}
