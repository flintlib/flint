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
    Scan A and fill in the min and max exponents of each variable along
    with the count of terms attached to each.
*/
static void _init_info(mpoly_gcd_var_info_struct * info, const fmpz_mpoly_t A,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    ulong * Aexp;
    slong i, j, N;
    slong nvars = ctx->minfo->nvars;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->bits <= FLINT_BITS);

    TMP_START;

    Aexp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    i = 0;
    mpoly_get_monomial_ui(Aexp, A->exps + N*i, A->bits, ctx->minfo);
    for (j = 0; j < nvars; j++)
    {
        info[j].min_exp = Aexp[j];
        info[j].max_exp = Aexp[j];
        info[j].min_exp_term_count = 1;
        info[j].max_exp_term_count = 1;
    }
    for (i = 1; i < A->length; i++)
    {
        mpoly_get_monomial_ui(Aexp, A->exps + N*i, A->bits, ctx->minfo);

        for (j = 0; j < nvars; j++)
        {
            if (info[j].min_exp > Aexp[j])
            {
                info[j].min_exp = Aexp[j];
                info[j].min_exp_term_count = 1;            
            }
            else if (info[j].min_exp == Aexp[j])
            {
                info[j].min_exp_term_count += 1;
            }

            if (info[j].max_exp < Aexp[j])
            {
                info[j].max_exp = Aexp[j];
                info[j].max_exp_term_count = 1;            
            }
            else if (info[j].max_exp == Aexp[j])
            {
                info[j].max_exp_term_count += 1;
            }
        }
    }

    TMP_END;
}



/*
    If Gbits = 0, the function has no restriction on the bits into which
    it can pack its answer.
    If Gbits != 0, the function was called from an internal context
    that expects all of G,A,B to be packed with bits = Gbits <= FLINT_BITS.

    Both A and B have to be packed into bits <= FLINT_BITS

    return is 1 for success, 0 for failure.
*/
int _fmpz_mpoly_gcd(fmpz_mpoly_t G, mp_bitcnt_t Gbits,
                               const fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong nvars = ctx->minfo->nvars;
    mpoly_gcd_var_info_struct * Ainfo, * Binfo;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(Gbits <= FLINT_BITS);
    FLINT_ASSERT(Gbits == 0 || (Gbits == A->bits && Gbits == B->bits));

    if (A->length == 1)
    {
        return _fmpz_mpoly_gcd_monomial(G, Gbits, B, A, ctx);
    }
    else if (B->length == 1)
    {
        return _fmpz_mpoly_gcd_monomial(G, Gbits, A, B, ctx);
    }

    TMP_START;

    Ainfo = (mpoly_gcd_var_info_struct *) TMP_ALLOC(nvars
                                           *sizeof(mpoly_gcd_var_info_struct));
    Binfo = (mpoly_gcd_var_info_struct *) TMP_ALLOC(nvars
                                           *sizeof(mpoly_gcd_var_info_struct));
    _init_info(Ainfo, A, ctx);
    _init_info(Binfo, B, ctx);

    /* check if the cofactors could be monomials */
    if (A->length == B->length)
    {
        success = _fmpz_mpoly_gcd_monomial_cofactors_sp(G, Gbits,
                                                      A, Ainfo, B, Binfo, ctx);
        if (success)
        {
            goto cleanup;
        }
    }

    success = fmpz_mpoly_gcd_brown(G, A, B, ctx);

cleanup:

    TMP_END;

    return success;
}


int fmpz_mpoly_gcd(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_zero(A, ctx))
    {
        if (B->length == 0)
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

    if (A->bits <= FLINT_BITS && B->bits <= FLINT_BITS)
    {
        return _fmpz_mpoly_gcd(G, 0, A, B, ctx);
    }

    if (A->length == 1)
    {
        return _fmpz_mpoly_gcd_monomial(G, 0, B, A, ctx);
    }
    else if (B->length == 1)
    {
        return _fmpz_mpoly_gcd_monomial(G, 0, A, B, ctx);
    }
    else if (_fmpz_mpoly_gcd_monomial_cofactors(G, A, B, ctx))
    {
        return 1;
    }
    else
    {
    /*
        Since we do not have truly sparse GCD algorithms,
        if either of A or B cannot be packed into FLINT_BITS and neither is a
        monomial, assume that the gcd computation is hopeless and fail.
    */
        int success;
        int useAnew = 0;
        int useBnew = 0;
        fmpz_mpoly_t Anew;
        fmpz_mpoly_t Bnew;

        if (A->bits > FLINT_BITS)
        {
            fmpz_mpoly_init(Anew, ctx);
            useAnew = fmpz_mpoly_repack_bits(Anew, A, FLINT_BITS, ctx);
            if (!useAnew)
            {
                useAnew = 1;
                success = 0;
                goto cleanup;
            }
        }

        if (B->bits > FLINT_BITS)
        {
            fmpz_mpoly_init(Bnew, ctx);
            useBnew = fmpz_mpoly_repack_bits(Bnew, B, FLINT_BITS, ctx);
            if (!useBnew)
            {
                useBnew = 1;
                success = 0;
                goto cleanup;
            }
        }

        success = _fmpz_mpoly_gcd(G, 0, useAnew ? Anew : A,
                                        useBnew ? Bnew : B, ctx);

cleanup:
        if (useAnew)
            fmpz_mpoly_clear(Anew, ctx);
        if (useBnew)
            fmpz_mpoly_clear(Bnew, ctx);

        return success;
    }
}
