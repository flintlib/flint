/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "nmod_mpoly.h"

void nmod_mpoly_set_term_ui(nmod_mpoly_t poly,
                        ulong const * exp, ulong c, const nmod_mpoly_ctx_t ctx)
{
    slong i, N, index, exp_bits;
    ulong * cmpmask;
    ulong * packed_exp;
    mp_limb_t cr;
    int exists;
    TMP_INIT;

    TMP_START;

    exp_bits = mpoly_exp_bits_required_ui(exp, ctx->minfo);
    if (exp_bits > FLINT_BITS)
        flint_throw(FLINT_EXPOF, "Exponent overflow in fmpz_mpoly_set_term_fmpz");

    /* reallocate the number of bits of the exponents of the polynomial */
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    nmod_mpoly_fit_bits(poly, exp_bits, ctx);

    N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, poly->bits, ctx->minfo);

    packed_exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    /* pack exponent vector */
    mpoly_set_monomial_ui(packed_exp, exp, poly->bits, ctx->minfo);

    /* work out at what index term should be placed */
    exists = mpoly_monomial_exists(&index, poly->exps,
                                  packed_exp, poly->length, N, cmpmask);

    NMOD_RED(cr, c, ctx->ffinfo->mod);
    if (!exists) /* term with that exponent doesn't exist */
    {
        if (cr != 0) /* only set if coeff is nonzero */
        {       
            nmod_mpoly_fit_length(poly, poly->length + 1, ctx);

            /* shift coeffs and exps by one to make space */
            for (i = poly->length; i >= index + 1; i--)
            {
                poly->coeffs[i] = poly->coeffs[i - 1];
                mpoly_monomial_set(poly->exps + N*i, poly->exps + N*(i - 1), N);
            }

            mpoly_monomial_set(poly->exps + N*index, packed_exp, N);
            poly->coeffs[index] = cr;

            poly->length++;
        }
    } else if (cr == 0) /* zero coeff, remove term */
    {
        /* shift coeffs and exps by one to make space */
        for (i = index; i < poly->length - 1; i++)
        {
            poly->coeffs[i] = poly->coeffs[i + 1];
            mpoly_monomial_set(poly->exps + N*i, poly->exps + N*(i + 1), N);
        }
        _nmod_mpoly_set_length(poly, poly->length - 1, ctx);

    } else
    {
        /* term with that monomial exists, coeff is nonzero */
        poly->coeffs[index] = cr;
    }

   TMP_END; 
}
