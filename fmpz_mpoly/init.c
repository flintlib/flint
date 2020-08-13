/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_init(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
    /* default to MPOLY_MIN_BITS bits per exponent */
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = MPOLY_MIN_BITS;
}

void fmpz_mpoly_init3(fmpz_mpoly_t A, slong alloc, flint_bitcnt_t bits,
                                                    const fmpz_mpoly_ctx_t ctx)
{
   slong N = mpoly_words_per_exp(bits, ctx->minfo);

    /* sanitize alloc input */
    alloc = FLINT_MAX(alloc, WORD(0));

   if (alloc != 0)
   {
      A->coeffs = (fmpz *) flint_calloc(alloc, sizeof(fmpz));
      A->exps   = (ulong *) flint_malloc(alloc*N*sizeof(ulong));
   } else
   {
      A->coeffs = NULL;
      A->exps = NULL;
   }
   A->alloc = alloc;
   A->length = 0;
   A->bits = bits;
}

void fmpz_mpoly_init2(fmpz_mpoly_t A, slong alloc, const fmpz_mpoly_ctx_t ctx)
{
    /* default to MPOLY_MIN_BITS bits per exponent */
    fmpz_mpoly_init3(A, alloc, MPOLY_MIN_BITS, ctx);
}
