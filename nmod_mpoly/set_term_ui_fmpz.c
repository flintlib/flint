/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"


void _nmod_mpoly_set_term_ui_fmpz(nmod_mpoly_t poly,
                         ulong c, const fmpz * exp, const nmod_mpoly_ctx_t ctx)
{
    mp_bitcnt_t exp_bits;
    slong i, N, index;
    ulong * cmpmask;
    ulong * packed_exp;
    mp_limb_t cr;
    int exists;
    TMP_INIT;

    TMP_START;

    exp_bits = mpoly_exp_bits_required_fmpz(exp, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    nmod_mpoly_fit_bits(poly, exp_bits, ctx);

    N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, poly->bits, ctx->minfo);

    packed_exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    mpoly_set_monomial_fmpz(packed_exp, exp, poly->bits, ctx->minfo);
    exists = mpoly_monomial_exists(&index, poly->exps,
                                  packed_exp, poly->length, N, cmpmask);

    NMOD_RED(cr, c, ctx->ffinfo->mod);

    if (!exists)
    {
        if (cr != 0) /* make new term only if coeff is nonzero*/
        {
            nmod_mpoly_fit_length(poly, poly->length + 1, ctx);

            for (i = poly->length; i >= index + 1; i--)
            {
                poly->coeffs[i] = poly->coeffs[i - 1];
                mpoly_monomial_set(poly->exps + N*i, poly->exps + N*(i - 1), N);
            }

            poly->coeffs[index] = cr;
            mpoly_monomial_set(poly->exps + N*index, packed_exp, N);

            poly->length++; /* safe because length is increasing */
        }
    } else if (cr == 0) /* zero coeff, remove term */
    {
        for (i = index; i < poly->length - 1; i++)
        {
            poly->coeffs[i] = poly->coeffs[i + 1];
            mpoly_monomial_set(poly->exps + N*i, poly->exps + N*(i + 1), N);
        }

        _nmod_mpoly_set_length(poly, poly->length - 1, ctx);

    } else /* term with that monomial exists, coeff is nonzero */
    {
        poly->coeffs[index] = cr;
    }

   TMP_END; 
}

void nmod_mpoly_set_term_ui_fmpz(nmod_mpoly_t poly,
                        ulong c, const fmpz ** exp, const nmod_mpoly_ctx_t ctx)
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

    _nmod_mpoly_set_term_ui_fmpz(poly, c, newexp, ctx);

    for (i = 0; i < nvars; i++)
        fmpz_clear(newexp + i);

    TMP_END;
}
