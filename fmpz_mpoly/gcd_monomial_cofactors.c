/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"


/*
    scan to check if the cofactors G/A and G/B are monomials
    If so, compute the GCD and return 1. Otherwise, return 0.
    This is a special internal version that takes info filled in as part of the 
    general gcd routine. It also assumes that A and B have the same length
    and are packed into <= FLINT_BITS.

    If Gbits = 0, the function has no restriction on the bits into which
    it can pack its answer.
    If Gbits != 0, the function was called from an internal context
    that expects all of G,A,B to be packed with bits = Gbits <= FLINT_BITS.

*/
int _fmpz_mpoly_gcd_monomial_cofactors_sp(fmpz_mpoly_t G, mp_bitcnt_t Gbits,
                 const fmpz_mpoly_t A, const mpoly_gcd_var_info_struct * Ainfo,
                 const fmpz_mpoly_t B, const mpoly_gcd_var_info_struct * Binfo,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i, j;
    slong NA, NB;
    slong nvars = ctx->minfo->nvars;
    ulong * Aexp, * Bexp;
    fmpz_t t1, t2;
    fmpz_t gA, gB;
    fmpz_mpoly_t T;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->length == B->length);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);

    for (j = 0; j < nvars; j++)
    {
        if (   Ainfo[j].max_exp - Ainfo[j].min_exp
            != Binfo[j].max_exp - Binfo[j].min_exp)
        {
            return 0;
        }
    }

    TMP_START;

    Aexp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
    Bexp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));

    NA = mpoly_words_per_exp(A->bits, ctx->minfo);
    NB = mpoly_words_per_exp(A->bits, ctx->minfo);
    
    fmpz_init(t1);
    fmpz_init(t2);
    fmpz_init(gA);
    fmpz_init(gB);
    success = 0;
    for (i = 0; i < A->length; i++)
    {
        fmpz_mul(t1, A->coeffs + 0, B->coeffs + i);
        fmpz_mul(t2, B->coeffs + 0, A->coeffs + i);
        if (!fmpz_equal(t1, t2))
        {
            success = 0;
            goto cleanup;
        }

        fmpz_gcd(gA, gA, A->coeffs + i);
        fmpz_gcd(gB, gB, B->coeffs + i);

        mpoly_get_monomial_ui(Aexp, A->exps + NA*i, A->bits, ctx->minfo);
        mpoly_get_monomial_ui(Bexp, B->exps + NB*i, B->bits, ctx->minfo);
        for (j = 0; j < nvars; j++)
        {
            if (   Aexp[j] - Ainfo[j].min_exp
                != Bexp[j] - Binfo[j].min_exp)
            {
                success = 0;
                goto cleanup;
            }
        }
    }

    success = 1;

    fmpz_mpoly_init3(T, A->length, FLINT_MIN(A->bits, B->bits), ctx);

    /* put A's cofactor monomial in Bexp */
    for (j = 0; j < nvars; j++)
    {
        Bexp[j] = FLINT_MAX(Ainfo[j].min_exp, Binfo[j].min_exp)
                    - Binfo[j].min_exp;
    }

    /* put A's cofactor coefficient in t1 */
    fmpz_gcd(t2, gA, gB);
    fmpz_divexact(t1, gA, t2);
    if (fmpz_sgn(A->coeffs + 0) < 0)
        fmpz_neg(t1, t1);

    for (i = 0; i < A->length; i++)
    {
        mpoly_get_monomial_ui(Aexp, A->exps + NA*i, A->bits, ctx->minfo);
        for (j = 0; j < nvars; j++)
        {
            Aexp[j] -= Bexp[j];
        }
        fmpz_divexact(t2, A->coeffs, t1);
        _fmpz_mpoly_emplacebackterm_fmpz_ui(T, t2, Aexp, ctx);
        fmpz_init(t2);
    }

    fmpz_mpoly_swap(T, G, ctx);
    fmpz_mpoly_clear(T, ctx);

cleanup:

    TMP_END;

    fmpz_clear(t1);
    fmpz_clear(t2);
    fmpz_clear(gA);
    fmpz_clear(gB);

    return success;
}


/*
    General slow version for arbitrary A and B.
*/
int _fmpz_mpoly_gcd_monomial_cofactors(fmpz_mpoly_t G,
                                const fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    return 0;
}


