/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void _fmpz_mpoly_set_coeff_fmpz_fmpz(fmpz_mpoly_t poly,
                  const fmpz_t c, const fmpz * exp, const fmpz_mpoly_ctx_t ctx)
{
    flint_bitcnt_t exp_bits;
    slong i, N, index;
    ulong * cmpmask;
    ulong * packed_exp;
    int exists;
    TMP_INIT;

    TMP_START;

    exp_bits = mpoly_exp_bits_required_ffmpz(exp, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    fmpz_mpoly_fit_bits(poly, exp_bits, ctx);

    N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, poly->bits, ctx->minfo);

    packed_exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    mpoly_set_monomial_ffmpz(packed_exp, exp, poly->bits, ctx->minfo);
    exists = mpoly_monomial_exists(&index, poly->exps,
                                  packed_exp, poly->length, N, cmpmask);

    if (!exists)
    {
        if (!fmpz_is_zero(c)) /* make new term only if coeff is nonzero*/
        {       

            fmpz_mpoly_fit_length(poly, poly->length + 1, ctx);

            for (i = poly->length; i >= index + 1; i--)
            {
                fmpz_set(poly->coeffs + i, poly->coeffs + i - 1);
                mpoly_monomial_set(poly->exps + N*i, poly->exps + N*(i - 1), N);
            }

            fmpz_set(poly->coeffs + index, c);
            mpoly_monomial_set(poly->exps + N*index, packed_exp, N);

            poly->length++; /* safe because length is increasing */
        }
    } else if (fmpz_is_zero(c)) /* zero coeff, remove term */
    {
        for (i = index; i < poly->length - 1; i++)
        {
            fmpz_set(poly->coeffs + i, poly->coeffs + i + 1);
            mpoly_monomial_set(poly->exps + N*i, poly->exps + N*(i + 1), N);
        }

        _fmpz_mpoly_set_length(poly, poly->length - 1, ctx);

    } else /* term with that monomial exists, coeff is nonzero */
    {
        fmpz_set(poly->coeffs + index, c);  
    }

   TMP_END; 
}


void fmpz_mpoly_set_coeff_fmpz_fmpz(fmpz_mpoly_t poly,
                const fmpz_t c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong i, nvars = ctx->minfo->nvars;
    fmpz * newexp;
    TMP_INIT;

    TMP_START;
    newexp = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (i = 0; i < nvars; i++)
    {
        fmpz_init(newexp + i);
        fmpz_set(newexp + i, exp[i]);
    }

    _fmpz_mpoly_set_coeff_fmpz_fmpz(poly, c, newexp, ctx);

    for (i = 0; i < nvars; i++)
        fmpz_clear(newexp + i);

    TMP_END;
}
